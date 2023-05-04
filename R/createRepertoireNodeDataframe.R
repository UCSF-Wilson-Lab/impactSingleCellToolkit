#' @title Create Immune Repertoire Network Node Dataframe
#' @description Takes an input dataframe of VDJ results (from Immcantation pipeline) and creates a node dataframe which can be used to construct an Immune Repertoire network
#' @param imm_results results dataframe outputted from immcantation and merged using mergeVDJresults function
#' @param patient.col Column name containing the patient ID <default: 'patient'>
#' @param seqid.col Column name containing the sequence ID <default: 'unique_cell_id'>
#' @param locus.col Column name which contains the type of BCR or TCR chain <default: 'locus'>
#' @param clone.col Column name which contains the clone ID <default: 'clone_id'>
#' @param cdr3.aa.col Column name which contains the AA sequences of the junction including the CDR3 region <default: 'junction_aa'>
#' @param umi.count.col Column name which contains the UMI count per contig <default: 'umi_count'>
#' @param read.count.col Column name which contains the read count per contig <default: 'consensus_count'>
#' @param prefix.clone.id String prefix added to values in the clone ID column <default: 'clone'> 
#' @param assay what assay type is being used among only two options (bcr,tcr) <default: 'bcr'>
#' @param select.patients  A vector of select patient ID used to subset the input dataframe <default: NULL>
#' @param expansion.thresh Minimum cell count threshold required for a clonotype to be defined as expanded <default: 2>
#' @return A dataframe with node IDs that can be used to create a network.
#' @examples node_df <- createRepertoireNodeDataframe(imm_results = bcr_df,expansion.thresh = 2,assay = "bcr")
#' @import data.table
#' @import stringr
#' @import dplyr
#' @export

createRepertoireNodeDataframe <- function(imm_results,
                                          patient.col = "patient",
                                          seqid.col = "unique_cell_id",
                                          locus.col = "locus",
                                          clone.col = "clone_id",
                                          cdr3.aa.col = "junction_aa",
                                          umi.count.col = "umi_count",
                                          read.count.col = "consensus_count",
                                          prefix.clone.id = "clone",
                                          assay = "bcr",
                                          select.patients = NULL,
                                          expansion.thresh = 2) {
  # Immcantation outputted columns to keep
  COLS_TO_KEEP = c(seqid.col,patient.col,clone.col,cdr3.aa.col,umi.count.col,read.count.col)
  
  ### 1. subset immcantation results
  if(!is.null(select.patients)){
    imm_results <- imm_results[imm_results[,patient.col] %in% select.patients,]
  }
  # Filter to one chain
  if(assay == "bcr"){imm_results <- imm_results[imm_results[,locus.col] %in% "IGH",]}
  if(assay == "tcr"){imm_results <- imm_results[imm_results[,locus.col] %in% "TRB",]}

  imm_results <- imm_results[,COLS_TO_KEEP]
  
  ### 2. Create clone count table
  clone_count_table             <- as.data.frame(table(imm_results[,clone.col]))
  names(clone_count_table)      <- c(clone.col,"count")
  clone_count_table$count       <- as.numeric(clone_count_table$count)
  clone_count_table[,clone.col] <- as.character(clone_count_table[,clone.col])
  clone_count_table$exp_status  <- clone_count_table$count
  clone_count_table$exp_status[clone_count_table$exp_status >= expansion.thresh] <- "Expanded"
  clone_count_table$exp_status[! clone_count_table$exp_status %in% "Expanded"] <- "Non_Expanded"
  
  clone_count_table$clone_id_fmt <- paste(prefix.clone.id,clone_count_table[,clone.col],sep = "")
  
  ### 3. initialize Node df
  nseqs <- nrow(imm_results)
  
  pt_nodes <- data.frame(matrix(nrow = nseqs,ncol = 5))
  names(pt_nodes) <- c("id","clone_id","clone_type","clone_type_label","clone_size")
  pt_nodes$id <- imm_results[,seqid.col]
  pt_nodes$clone_id <- imm_results[,clone.col]
  pt_nodes$clone_id <- paste(prefix.clone.id,pt_nodes$clone_id,sep = "")
  
  
  ### 4. Add extra info to node df (expansion status, clone count, etc.)
  
  # Conversion lists from clone table
  clone_status_list        <- clone_count_table$exp_status
  names(clone_status_list) <- clone_count_table$clone_id_fmt
  
  clone_status_num_list    <- clone_status_list
  clone_status_num_list[clone_status_num_list %in% "Expanded"]     <- 2
  clone_status_num_list[clone_status_num_list %in% "Non_Expanded"] <- 1
  
  clone_count_list        <- as.character(clone_count_table$count)
  names(clone_count_list) <- clone_count_table$clone_id_fmt
  
  clone_status_list        <- as.list(clone_status_list)
  clone_status_num_list    <- as.list(clone_status_num_list)
  clone_count_list         <- as.list(clone_count_list)
  
  # Clonotype labels and size columns
  pt_nodes$clone_type_label <- createNodeCol(pt_nodes$clone_id,clone_status_list)
  pt_nodes$clone_type       <- createNodeCol(pt_nodes$clone_id,clone_status_num_list)
  pt_nodes$clone_size       <- createNodeCol(pt_nodes$clone_id,clone_count_list)
  
  return(pt_nodes)
}


# Function for formatting columns
createNodeCol <- function(clone_col, conv_list) {
  uni_clones <- unique(clone_col)
  for (clone in uni_clones) {
    clone <- as.character(clone)
    categ <- conv_list[[clone]]
    clone_col[clone_col %in% clone] <- categ
  }
  
  return(clone_col)
}