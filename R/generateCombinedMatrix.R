#' @title Create single cell count matrix
#' @description Takes input directory of cellranger (10x) and combines all samples into one sparse matrix.
#' This function only works for scRNA-Seq and Cite-Seq data. Each sample sub-directory must contain:
#' - barcodes.tsv.gz
#' - features.tsv.gz
#' - matrix.mtx.gz
#' @param dataset_loc directory path for all cellranger results
#' @param samples.vec character vector of all sample names
#' @param THREADS number of threads <default: 4>
#' @param multi.results bool value for where these are results outputted by cellranger multi <default: TRUE>
#' @param assay what assay type is being used among only two options (gex, csp) <default: 'gex'>
#' @param min.genes.per.cell minimum number of genes per cell with counts greater than zero <default: 400>
#' @param max.genes.per.cell maximum number of genes per cell with counts greater than zero <default: NULL>
#' @return sparse matrix where cells are the columns and genes are the rows.
#' @examples generateCombinedMatrix(dataset_loc, samples.vec,THREADS = 15, multi.results=T,
#' assay = "gex",min.genes.per.cell = 700,max.genes.per.cell = 2500)
#' @import Seurat
#' @import parallel
#' @export

generateCombinedMatrix <- function(dataset_loc, samples.vec, THREADS = 4, multi.results=T,assay = "gex",
                                   min.genes.per.cell = 400,max.genes.per.cell = NULL) {
  # Access More clusters
  cl <- makeCluster(getOption("cl.cores", THREADS))
  clusterExport(cl,c('Read10X', "dataset_loc", "samples.vec"))
  start.time <- Sys.time()
  
  # Load sample Matricies
  d10x.data <- sapply(samples.vec, function(i){
    data.dir <- file.path(dataset_loc,i,"/")
    d10x <- Read10X(data.dir = data.dir)
    if(assay == "gex")
      if(multi.results == TRUE){d10x <- d10x$`Gene Expression`}
    colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
    if(assay == "csp"){
      if(multi.results == TRUE){d10x <- d10x$`Antibody Capture`}
      colnames(d10x) <- paste(sapply(strsplit(colnames(d10x),split="-"),'[[',1L),i,sep="-")
    }
    d10x
  })
  
  # Stop running mult-threads
  end.time <- Sys.time()
  time.taken <- end.time - start.time
  time.taken
  
  stopCluster(cl)
  rm(cl)
  
  # Combine all patient matricies
  experiment.data <- do.call("cbind", d10x.data)
  
  # Filter cells based on number genes expressed per cell
  num_genes_per_cell <- (colSums(experiment.data !=0))
  target_bool <- (num_genes_per_cell > min.genes.per.cell)
  if(!is.null(max.genes.per.cell)){
    target_bool <- (num_genes_per_cell > min.genes.per.cell & num_genes_per_cell < max.genes.per.cell)
  }
  cells_to_keep <- target_bool[target_bool == T]
  experiment.data <- experiment.data[,colnames(experiment.data) %in% names(cells_to_keep)]
  
  return(experiment.data)
}
