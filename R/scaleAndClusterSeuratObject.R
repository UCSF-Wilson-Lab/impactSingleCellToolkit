#' @title Scale and Cluster a Seurat Object
#' @description Takes an input Seurat Object and applies normalization (regression on UMI count) and dimension reduction
#' @param seurat_obj Seurat Object
#' @param dims Dimensions of reduction used as input <default: 1:30>
#' @param npca  Maximum dimensions of reduction for tSNE if number of cells are too small <default: 10>
#' @param normalize Bool of whether to normalize using the SCT transform <default: TRUE>
#' @param dim.reduction Bool of whether to perform dimension reduction with PCA, tSNE and UMAP <default: TRUE>
#' @param tsne Bool of whether to perform dimension reduction with tSNE <default: TRUE>
#' @param umap Bool of whether to perform dimension reduction with UMAP <default: TRUE>
#' @return A re-clustered Seurat object.
#' @references 
#' Hao Y, Hao S, et al (2021). 
#' Integrated analysis of multimodal single-cell data.
#' \emph{Cell}. doi:10.1016/j.cell.2021.04.048 
#' @examples seurat.obj <- scaleAndClusterSeuratObject(seurat.obj,npca = 20,normalize = TRUE,dim.reduction = TRUE)
#' @import Seurat
#' @import parallel
#' @export

scaleAndClusterSeuratObject <- function(seurat_obj,dims = 1:30,npca = 10,
                                        normalize = TRUE,dim.reduction = TRUE,tsne = TRUE, umap = TRUE) 
{
  if(normalize){
    seurat_obj <- SCTransform(object = seurat_obj, verbose = FALSE, do.correct.umi = T, vars.to.regress = "nCount_RNA")
  }
  
  ## Dimension reduction with PCA, tSNE and UMAP
  if(dim.reduction){
    seurat_obj <- RunPCA(object = seurat_obj, verbose = FALSE)
    seurat_obj <- FindNeighbors(object = seurat_obj, dims = dims)
    seurat_obj <- FindClusters(object = seurat_obj)
    
    if(tsne){
      if (ncol(seurat_obj@assays$RNA@data) < 100) {
        seurat_obj <- RunTSNE(seurat_obj, perplexity = 10, dims = 1:npca)
      }
      else {
        seurat_obj <- RunTSNE(seurat_obj, dims = dims)
      }
    }
    
    if(umap){
      seurat_obj <- RunUMAP(seurat_obj, dims = dims, verbose = FALSE)
    }
  }
  
  return(seurat_obj)
}
