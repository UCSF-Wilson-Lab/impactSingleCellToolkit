#' @title Visualize the optimal PCs for Seurat Object
#' @description Takes an input Seurat Object and generates a plot of variation captured across different numbers of principal components (PC)
#' @param seurat_obj Seurat Object
#' @return A plot showing increasing numbers of PCs where each point represents PCs which cumulatively account for greater than 90% of variation and less than 5% std. dev. of variation between PCs.
#' @references This code is from https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_clustering_analysis.html
#' @examples plot <- plotOptimalPCsforSeuratObject(seurat_obj)
#' @import Seurat
#' @export

plotOptimalPCsforSeuratObject <- function(seurat_obj) {
  # Source: https://hbctraining.github.io/scRNA-seq/lessons/sc_exercises_clustering_analysis.html
  pct <- seurat_obj[["pca"]]@stdev / sum(seurat_obj[["pca"]]@stdev) * 100
  # Calculate cumulative percents for each PC
  cumu <- cumsum(pct)
  # Determine which PC exhibits cumulative percent greater than 90% and % variation associated with the PC as less than 5
  co1 <- which(cumu > 90 & pct < 5)[1]
  co2 <- sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1), decreasing = T)[1] + 1
  
  # Minimum of the two calculation
  pcs <- min(co1, co2)
  
  plot_df <- data.frame(pct = pct, 
                        cumu = cumu, 
                        rank = 1:length(pct))
  # Elbow plot to visualize 
  plot <- ggplot(plot_df, aes(cumu, pct, label = rank, color = rank > pcs)) + 
    geom_text() + geom_vline(xintercept = 90, color = "grey") + 
    geom_hline(yintercept = min(pct[pct > 5]), color = "grey") + theme_bw()
  
  return(plot)
}
