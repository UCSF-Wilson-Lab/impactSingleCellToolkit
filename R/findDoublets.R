#' @title Annotate Doublets inside a Seurat Object
#' @description Takes an input Seurat Object and annotates all doublets using the DoubletFinder R package
#' @param seurat.obj Seurat Object
#' @param sample.col Metadata column name containing sample names within the input Seurat Object <default: 'sample'>
#' @param threads  number of threads <default: 5>
#' @return A Seurat object with a new column 'doublet_finder' which contains doublet annotations.
#' @examples seurat.obj <- findDoublets(seurat.obj, sample.col = "sample", threads = 10)
#' @import Seurat
#' @import DoubletFinder
#' @import parallel
#' @export

findDoublets <- function(seurat.obj,sample.col = "sample", threads = 5) {
  scobject.split <- SplitObject(seurat.obj, split.by = sample.col) 
  
  # Doublet finder column
  doublet_finder_column <- row.names(seurat.obj@meta.data)
  
  # loop through samples to find doublets
  for (i in 1:length(scobject.split)) {
    # Pre-process seurat object with standard seurat workflow
    sc.sample <- NormalizeData(scobject.split[[i]])
    sc.sample <- FindVariableFeatures(sc.sample)
    sc.sample <- ScaleData(sc.sample)
    sc.sample <- RunPCA(sc.sample, nfeatures.print = 10)
    
    # Find significant PCs
    stdv <- sc.sample[["pca"]]@stdev
    sum.stdv <- sum(sc.sample[["pca"]]@stdev)
    percent.stdv <- (stdv / sum.stdv) * 100
    cumulative <- cumsum(percent.stdv)
    co1 <- which(cumulative > 90 & percent.stdv < 5)[1]
    co2 <- sort(which((percent.stdv[1:length(percent.stdv) - 1] - 
                         percent.stdv[2:length(percent.stdv)]) > 0.1), 
                decreasing = T)[1] + 1
    min.pc <- min(co1, co2)
    
    # finish pre-processing
    sc.sample <- RunUMAP(sc.sample, dims = 1:min.pc)
    sc.sample <- FindNeighbors(object = sc.sample, dims = 1:min.pc)              
    sc.sample <- FindClusters(object = sc.sample, resolution = 0.1)
    
    # pK identification (no ground-truth)
    sweep.list <- paramSweep_v3(sc.sample, PCs = 1:min.pc,num.cores = threads,sct = F)
    sweep.stats <- summarizeSweep(sweep.list,GT = F)
    bcmvn <- find.pK(sweep.stats)
    
    # Optimal pK is the max of the bomodality coefficent (BCmvn) distribution
    bcmvn.max <- bcmvn[which.max(bcmvn$BCmetric),]
    optimal.pk <- bcmvn.max$pK
    optimal.pk <- as.numeric(levels(optimal.pk))[optimal.pk]
    
    # Get 10x doublet rate
    sample_cell_count <- nrow(sc.sample@meta.data)
    doublet.rate      <- get10xDoubletRate(sample_cell_count)
    
    ## Homotypic doublet proportion estimate
    annotations <- sc.sample@meta.data$seurat_clusters
    homotypic.prop <- modelHomotypic(annotations) 
    nExp.poi <- round(doublet.rate * nrow(sc.sample@meta.data))
    nExp.poi.adj <- round(nExp.poi * (1 - homotypic.prop))
    
    # run DoubletFinder
    sc.sample <- doubletFinder_v3(seu = sc.sample, 
                                  PCs = 1:min.pc, 
                                  pK = optimal.pk,
                                  nExp = nExp.poi.adj)
    metadata <- sc.sample@meta.data
    colnames(metadata)[ncol(metadata)] <- "doublet_finder"
    sc.sample@meta.data <- metadata 
    
    # subset and store in doublet finder column
    sc.singlets <- subset(sc.sample, doublet_finder == "Singlet")
    scobject.split[[i]] <- sc.singlets
    
    # Add to doublet annotation column
    doublet_finder_column[doublet_finder_column %in% colnames(sc.singlets)] <- "Singlet"
    
    remove(sc.singlets)
  }
  rm(scobject.split)
  
  # Add DoubletFinder results back into Seurat Object
  doublet_finder_column[! doublet_finder_column %in% "Singlet"] <- "Doublet"
  seurat.obj@meta.data$doublet_finder <- doublet_finder_column
  
  return(seurat.obj)
}


#' @title Calculate doublet rate based on sample cell count
#' @description Takes and input sample cell count and calculates a doublet rate using published doublet rates from the 3' scRNA-Seq assay (10x Genomics)
#' @param sample_cell_count Seurat Object
#' @return A sample doublet rate which will get used in the 'findDoublets' function.
#' @examples sample.doublet.rate <- get10xDoubletRate(sample_cell_count = 1000)
#' @import Seurat
#' @export

get10xDoubletRate <- function(sample_cell_count) {
  min.doublet.rate = 0.004
  max.doublet.rate = 0.080
  
  # doublet rates for 3 prime kit 
  multiplet_rates <- list()
  multiplet_rates[[1000]]  <- 0.008
  multiplet_rates[[2000]]  <- 0.016
  multiplet_rates[[3000]]  <- 0.024
  multiplet_rates[[4000]]  <- 0.032
  multiplet_rates[[5000]]  <- 0.040
  multiplet_rates[[6000]]  <- 0.048
  multiplet_rates[[7000]]  <- 0.056
  multiplet_rates[[8000]]  <- 0.064
  multiplet_rates[[9000]]  <- 0.072
  multiplet_rates[[10000]] <- 0.080
  
  
  rounded_count <- round(sample_cell_count)
  if(sample_cell_count >= 1000){
    rounded_count <- round(sample_cell_count / 1000) * 1000
  }
  
  if (rounded_count < 1000){
    doublet.rate <- min.doublet.rate
  } else{
    if (rounded_count <=10000){
      doublet.rate <- multiplet_rates[[rounded_count]]
    } else{
      doublet.rate <- max.doublet.rate
    }
  }
  
  return(doublet.rate)
}

