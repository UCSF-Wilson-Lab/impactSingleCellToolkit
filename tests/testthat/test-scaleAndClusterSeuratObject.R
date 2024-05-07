context("scaleAndClusterSeuratObject")

library(Seurat)
library(SeuratData)
pbmc       <- LoadData("pbmc3k", type = "default")
pbmc_small <- pbmc[,colnames(pbmc)[1:200]]

testthat::test_that("scaleAndClusterSeuratObject works", {
  set.seed(0)
  
  pbmc_small <- scaleAndClusterSeuratObject(pbmc_small,normalize = FALSE,dim.reduction = FALSE,tsne = FALSE, umap = FALSE)
  expect_equal(ncol(pbmc_small), 200)
})
