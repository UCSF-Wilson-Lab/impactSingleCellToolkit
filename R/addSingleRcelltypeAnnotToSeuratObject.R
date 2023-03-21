#' @title Add SingleR cell-type annotations into Seurat Object
#' @description Takes an input Seurat Object and adds inputted SingleR cell-type annotations into the metadata of the object
#' @param seurat Seurat Object
#' @param singler_results_list list of SingleR results objects
#' @param cluster_col Column name with the Seurat Object metadata that contains the cluster IDs <default: 'seurat_clusters'>
#' @return Seurat Object with SingleR annotations.
#' @examples seurat <- addSingleRcelltypeAnnotToSeuratObject(seurat,singler_results_list)
#' @import Seurat
#' @import SingleR
#' @export

addSingleRcelltypeAnnotToSeuratObject <- function(seurat,singler_results_list,
                                                  cluster_col = "seurat_clusters")
{
  annot_categories <- names(singler_results_list)
  
  for (annot in annot_categories) {
    clus_df   <- seurat[[cluster_col]]
    clust_col <- as.character(clus_df[,cluster_col])
    names(clust_col) <- row.names(clus_df)
    
    results_df          <- singler_results_list[[annot]]
    annot_labels        <- results_df$labels
    names(annot_labels) <- row.names(results_df)
    for (c in names(annot_labels)) {
      celltype <- as.character(annot_labels[names(annot_labels) %in% c])
      clust_col[clust_col %in% c] <- celltype
    }
    
    seurat[[annot]]     <- clust_col
  }
  
  return(seurat)
}
