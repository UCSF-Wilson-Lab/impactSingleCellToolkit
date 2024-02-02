#' @title Create Immune Repertoire Network Edge Dataframe
#' @description Takes an input node dataframe of VDJ results (from createRepertoireNodeDataframe) and creates an edge dataframe which can be used to construct an Immune Repertoire network
#' @param node_df node dataframe outputted from createRepertoireNodeDataframe
#' @param interaction.type Column name containing the interaction type <default: 'hyperlink'>
#' @return A dataframe with all pairwise edge interactions used to create a network.
#' @examples edge_df <- createRepertoireEdgeDataframe(node_df)
#' @import data.table
#' @import stringr
#' @import dplyr
#' @export

createRepertoireEdgeDataframe <- function(node_df, 
                                          interaction.type = "hyperlink") {
  EDGE_COLS = c("from","to","type","weight")
  
  # Initialize edge df
  pt_links <- data.frame(matrix(nrow = 0,ncol = 4))
  names(pt_links) <- EDGE_COLS
  
  # create clone count list
  clone_count_list <- list()
  clone_count_vec <- unique(paste(node_df$clone_id,node_df$clone_size,sep = ":sep:"))
  for (id_count in clone_count_vec) {
    info_vec <- unlist(str_split(id_count,":sep:"))
    clone_count_list[[info_vec[1]]] <- as.numeric(info_vec[2])
  }
  
  # Only make connections for Expanded clones
  exp_clones <- node_df[node_df$clone_type_label %in% "Expanded",]
  nonexp_clones <- node_df[node_df$clone_type_label %in% "Non_Expanded",]
  
  if(nrow(exp_clones) > 0){
    for (clone in unique(exp_clones$clone_id)) {
      clone_interacs <- getAllUniPairwiseInteractions(node_df,clone)
      
      # Create entry for all clone interactions
      clone_links <- data.frame(matrix(nrow = length(clone_interacs),ncol = 4))
      names(clone_links) <- EDGE_COLS
      clone_links$from <- tstrsplit(clone_interacs,":interac:")[[1]]
      clone_links$to <- tstrsplit(clone_interacs,":interac:")[[2]]
      
      clone_links$type   <- interaction.type
      clone_links$weight <- clone_count_list[[clone]]
      
      # rbind
      pt_links <- rbind(pt_links,clone_links)
    }
    pt_links$edge_arrows <- 0
  }
  
  return(pt_links)
}


# getAllUniPairwiseInteractions
#  - list all pairwise clonotype interactions for edge dataframe
getAllUniPairwiseInteractions <- function(pt_nodes,clone) {
  clone_df  <- pt_nodes[pt_nodes$clone_id %in% clone,]
  node_list <- clone_df$id
  list_len  <- length(node_list)
  
  all_uni_pairs <- c()
  for (i in 1:(list_len-1)) {
    curr_node <- node_list[i]
    rem_nodes <- node_list[(i+1):list_len]
    pairs     <- paste(curr_node,rem_nodes,sep = ":interac:")
    all_uni_pairs <- c(all_uni_pairs,pairs)
  }
  
  all_uni_pairs <- unique(all_uni_pairs)
  
  return(all_uni_pairs)
}
