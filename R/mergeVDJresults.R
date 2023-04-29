#' @title Merge Immcantation VDJ results
#' @description Takes input single cell immcantation results dataframes and merges them so that each cell has two chain assemblies.
#' Each BCR should have one heavy and one light chain. Each TCR should have one alpha and one beta chain.
#' @param df1 Immcantation results data frame either BCR heavy chains or TCR Beta chains <default: NULL>
#' @param df2 Immcantation results data frame either BCR light chains or TCR Alpha chains <default: NULL>
#' @param umi.thresh UMI count threshold for all contigs <default: 2>
#' @param assay what assay type is being used among only two options (bcr,tcr) <default: 'bcr'>
#' @return merge BCR or TCR dataframe.
#' @examples bcr_merged_df <- mergeVDJresults(df1 = bcr_heavy_df, df2 = bcr_light_df,umi.thresh = 3,assay = "bcr")
#' @import stringr
#' @import Matrix
#' @import BiocGenerics
#' @import data.table
#' @export
mergeVDJresults <- function(df1,
                            df2,
                            umi.thresh = 2,
                            assay = "bcr") {
  # create unique cell id
  df1$unique_cell_id <- tstrsplit(df1$sequence_id,"-")[[1]]
  df2$unique_cell_id <- tstrsplit(df2$sequence_id,"-")[[1]]
  
  # apply umi count filter
  df1 <- df1[df1$umi_count >= umi.thresh,]
  df2 <- df2[df2$umi_count >= umi.thresh,]
  
  # Filter for specific vdj loci
  if(assay == "bcr"){
    df1 <- df1[df1$locus %in% "IGH",]
    df2 <- df2[df2$locus %in% c("IGK","IGL"),]
  }
  if(assay == "tcr"){
    df1 <- df1[df1$locus %in% "TRB",]
    df2 <- df2[df2$locus %in% "TRA",]
  }
  
  # Filter contigs so each unique ID has two chains
  df1 <- filterVDJcontigs(df1,unique.id.col = "unique_cell_id")
  df2 <- filterVDJcontigs(df2,unique.id.col = "unique_cell_id")
  
  # one keep cells that overlap between both data frames
  overlapping_ids <- df1$unique_cell_id[df1$unique_cell_id %in% df2$unique_cell_id]
  df1             <- df1[df1$unique_cell_id %in% overlapping_ids,]
  df2             <- df2[df2$unique_cell_id %in% overlapping_ids,]
  
  # merge 
  df <- rbind(df1,df2)
  
  # QC filter to make sure there are only two chains per cell
  df <- qcfilterVDJcontigs(df,unique.id.col = "unique_cell_id",locus.col = "locus",assay=assay)
  
  # Adjust clonotype ID so it is uniform between chains (BCR or TCR)
  uni_cell_list <- as.list(unique(df$unique_cell_id))
  df$clone_id   <- unlist(lapply(uni_cell_list, formatCellCloneID,df=df,
                                 assay=assay,
                                 unique.id.col = "unique_cell_id",
                                 locus.col = "locus",
                                 clone.col = "clone_id"))

  return(df)
}


# Get max contig - used in filterVDJcontigs.cell.2chains
getMaxContig <- function(cell_df, chain,locus.col = "locus") {
  cell_df_sub <- cell_df[cell_df[,locus.col] %in% chain,]
  max_umi_contig  <- cell_df_sub[cell_df_sub$umi_count == max(cell_df_sub$umi_count),]
  max_umi_contig  <- max_umi_contig[max_umi_contig$consensus_count == max(max_umi_contig$consensus_count),]
  max_umi_contig  <- max_umi_contig$sequence_id
  # if UMI and reads counts are the same then arbitrality pick the first contig 
  if(length(max_umi_contig) > 1){max_umi_contig <- max_umi_contig[1]}
  
  return(max_umi_contig)
}


# Filter contigs for one cell all chains
filterVDJcontigs.cell.2chains <- function(cell,df,unique.id.col,locus.col,assay) {
  cell_df <- df[df[,unique.id.col] %in% cell,]
  chains  <- unique(cell_df[,locus.col])
  
  if(nrow(cell_df) > 1 && length(chains) > 1){
    seq_ids_keep <- c()
    
    if(assay == "bcr"){
      chain_heavy = "IGH"
      chain_light = c("IGK","IGL")
      heavy_contig <- getMaxContig(cell_df,chain = chain_heavy,locus.col=locus.col)
      light_contig <- getMaxContig(cell_df,chain = chain_light,locus.col=locus.col)
      
      seq_ids_keep <- c(heavy_contig,light_contig)
    }
    
    if(assay == "tcr"){
      chain_beta  = "TRB"
      chain_alpha = "TRA"
      beta_contig <- getMaxContig(cell_df,chain = chain_beta,locus.col=locus.col)
      alpha_contig <- getMaxContig(cell_df,chain = chain_alpha,locus.col=locus.col)
      
      seq_ids_keep <- c(beta_contig,alpha_contig)
    }
    
    cell_df <- cell_df[cell_df$sequence_id %in% seq_ids_keep,]
  }else{
    cell_df <- NULL
  }
  
  # return contigs to keep
  
  return(cell_df)
}


# Filter contigs for one cell one chain
filterVDJcontigs.cell <- function(cell,df,unique.id.col) {
  cell_df <- df[df[,unique.id.col] %in% cell,]

  if(nrow(cell_df) > 1){
    max_umi_contig  <- cell_df[cell_df$umi_count == max(cell_df$umi_count),]
    max_umi_contig  <- max_umi_contig[max_umi_contig$consensus_count == max(max_umi_contig$consensus_count),]
    max_umi_contig  <- max_umi_contig$sequence_id
    # if UMI and reads counts are the same then arbitrality pick the first contig 
    if(length(max_umi_contig) > 1){max_umi_contig <- max_umi_contig[1]}
    
    cell_df <- cell_df[cell_df$sequence_id %in% max_umi_contig,]
  }
  
  return(cell_df)
}

# Filter contigs for immcatation results one chain only (i.e. TCR beta chains)
filterVDJcontigs <- function(df,unique.id.col = "unique_cell_id") {
  uni_id_list <- as.list(unique(df[,unique.id.col]))
  
  df_filt <- lapply(uni_id_list, filterVDJcontigs.cell,
                    df=df,unique.id.col=unique.id.col)
  
  df_filt <- do.call("rbind",df_filt)
  
  return(df_filt)
}


# Filter contigs for immcatation results all chains (i.e. TCR alpha & beta chains)
qcfilterVDJcontigs <- function(df,unique.id.col = "unique_cell_id",locus.col = "locus",assay = "bcr") {
  uni_id_list <- as.list(unique(df[,unique.id.col]))
  
  df_filt <- lapply(uni_id_list, filterVDJcontigs.cell.2chains,
                    df=df,unique.id.col=unique.id.col,locus.col=locus.col,assay=assay)
  df_filt <- Filter(function(x) !is.null(x), df_filt)
  
  df_filt <- do.call("rbind",df_filt)
  
  return(df_filt)
}

# Format clonotype IDs so they are aligned across all cells
#  - function supports BCR and TCR data
#  - BCR - Heavy chain clonotype is mapped to light chains
#  - TCR - Alpha and beta chain clone IDs are merged
formatCellCloneID <- function(cell,df,assay = "bcr",
                              unique.id.col = "unique_cell_id",
                              locus.col = "locus",
                              clone.col = "clone_id") {
    cell_df <- df[df[,unique.id.col] %in% cell,]
    cell_df <- cell_df[,c(unique.id.col,locus.col,clone.col)]
    
    if(assay == "bcr"){
      clone_heavy         <- cell_df[cell_df[,locus.col] %in% "IGH",clone.col]
      cell_df[,clone.col] <- clone_heavy
    }
    
    if(assay == "tcr"){
      clone_alpha <- cell_df[cell_df[,locus.col] %in% "TRA",clone.col]
      clone_beta  <- cell_df[cell_df[,locus.col] %in% "TRB",clone.col]
      clone_ab    <- paste(clone_alpha,clone_beta,sep = ":")
      cell_df[,clone.col] <- clone_ab
    }
    
    return(cell_df[,clone.col])
}
