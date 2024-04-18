#' @title Generate cell-type annotations using SingleR
#' @description Takes an input Seurat Object and generates cell-type annotations for human data single cell RNA-Seq data using the SingleR R package
#' @param seurat_obj Seurat Object
#' @param blueprintencode.annot Bool of whether to generate celltype annotations using the Blueprint ENCODE dataset as a reference <default: TRUE>
#' @param monaco.annot Bool of whether to generate celltype annotations using the Monaco dataset as a reference <default: TRUE>
#' @return SingleR results object.
#' @references 
#' Aran D, Looney AP, Liu L et al. (2019).
#' Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage.
#' \emph{Nat. Immunology} 20, 163–172. 
#' @examples singler_results <- runSingleRcelltypeAnnot(seurat_obj)
#' @import Seurat
#' @import celldex
#' @import SingleCellExperiment
#' @import SingleR
#' @export

runSingleRcelltypeAnnot <- function(seurat_obj,
                                    blueprintencode.annot = TRUE,
                                    monaco.annot = TRUE) 
{
  seurat_obj[["RNA3"]] <- as(object = seurat_obj[["RNA"]], Class = "Assay") # Convert to Seurat v3 assay class
  sc                   <- as.SingleCellExperiment(seurat_obj)               # Supports Seurat v3 assay class over v5
  results_list         <- list()
  
  if (blueprintencode.annot) {
    cat(paste0(">> Running SingleR: BlueprintENCODE\n\n"))
    immun_ref <- celldex::BlueprintEncodeData()
    singler.obj.bencodeRef.fine <- SingleR(sc,ref = immun_ref,labels = immun_ref$label.fine,clusters = seurat_obj$seurat_clusters)
    singler.obj.bencodeRef.main <- SingleR(sc,ref = immun_ref,labels = immun_ref$label.main,clusters = seurat_obj$seurat_clusters)
    results_list[["main_bencode"]] <- singler.obj.bencodeRef.main
    results_list[["fine_bencode"]] <- singler.obj.bencodeRef.fine
  }
  
  if (monaco.annot) {
    cat(paste0(">> Running SingleR: Monaco\n\n"))
    immun_ref <- celldex::MonacoImmuneData()
    singler.obj.monacoRef.fine <- SingleR(sc,ref = immun_ref,labels = immun_ref$label.fine,clusters = seurat_obj$seurat_clusters)
    singler.obj.monacoRef.main <- SingleR(sc,ref = immun_ref,labels = immun_ref$label.main,clusters = seurat_obj$seurat_clusters)
    results_list[["main_monaco"]] <- singler.obj.monacoRef.main
    results_list[["fine_monaco"]] <- singler.obj.monacoRef.fine
  }
  
  return(results_list)
}
