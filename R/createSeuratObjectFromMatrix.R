#' @title Create Seurat Object from Count Matrix
#' @description Takes input single cell sparse matrix of gene counts and creates an analysis object using the Seurat R package
#' @param sc.data Input sparse matrix of gene counts per cell
#' @param project.name Character string of the project name that will be added into the Seurat object
#' @param min.genes Minimum number of genes with expression greater than zero per cell <default: 200>
#' @param min.cells Minimum number of cells required where features are detected <default: 2>
#' @param npca  Maximum dimensions of reduction for tSNE if number of cells are too small <default: 10>
#' @param dims Dimensions of reduction used as input <default: 1:30>
#' @param normalize Bool of whether to normalize using the SCT transform <default: FALSE>
#' @param dim.reduction Bool of whether to perform dimension reduction with PCA, tSNE and UMAP <default: FALSE>
#' @return Seurat object generated from input sparse matrix.
#' @references 
#' Hao Y, Hao S, et al (2021). 
#' Integrated analysis of multimodal single-cell data.
#' \emph{Cell}. doi:10.1016/j.cell.2021.04.048 
#' @examples seurat.obj <- createSeuratObjectFromMatrix(sc.data,project.name = "Analysis1",npca = 20, min.genes = 700)
#' @import Seurat
#' @import parallel
#' @export

createSeuratObjectFromMatrix <-  function (sc.data, project.name,min.genes = 200, min.cells = 2, 
                                           npca = 10, dims= 1:30, normalize = FALSE,dim.reduction=FALSE) 
{
  mtgenes = "^mt-"
  ## Create Object
  sc = CreateSeuratObject(sc.data, min.cells = min.cells, 
                          min.features = min.genes, project = project.name)
  ## Add % MT genes
  percent.mito <- PercentageFeatureSet(object = sc, pattern = "^(?i)mt-")
  sc <- AddMetaData(object = sc, metadata = percent.mito, 
                    col.name = "percent.mito")
  
  ## Normalize
  if(normalize){
    sc <- SCTransform(object = sc, verbose = FALSE, do.correct.umi = T, vars.to.regress = "nCount_RNA")
  }
  
  
  ## Dimension reduction with PCA, tSNE and UMAP
  if(dim.reduction){
    sc <- RunPCA(object = sc, verbose = FALSE)
    sc <- FindNeighbors(object = sc, dims = dims)
    sc <- FindClusters(object = sc)
    if (ncol(sc@assays$RNA$counts) < 100) {
      sc <- RunTSNE(sc, perplexity = 10, dims = 1:npca)
    }
    else {
      sc <- RunTSNE(sc, dims = dims)
    }
    sc <- RunUMAP(sc, dims = dims, verbose = FALSE)
  }
  
  return(sc)
}
