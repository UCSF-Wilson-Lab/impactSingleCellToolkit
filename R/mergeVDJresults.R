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
  
  return(df)
}


# Filter contigs for one cell
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

