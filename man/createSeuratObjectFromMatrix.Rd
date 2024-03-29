% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/createSeuratObjectFromMatrix.R
\name{createSeuratObjectFromMatrix}
\alias{createSeuratObjectFromMatrix}
\title{Create Seurat Object from Count Matrix}
\usage{
createSeuratObjectFromMatrix(
  sc.data,
  project.name,
  min.genes = 200,
  min.cells = 2,
  npca = 10,
  dims = 1:30,
  normalize = FALSE,
  dim.reduction = FALSE
)
}
\arguments{
\item{sc.data}{Input sparse matrix of gene counts per cell}

\item{project.name}{Character string of the project name that will be added into the Seurat object}

\item{min.genes}{Minimum number of genes with expression greater than zero per cell <default: 200>}

\item{min.cells}{Minimum number of cells required where features are detected <default: 2>}

\item{npca}{Maximum dimensions of reduction for tSNE if number of cells are too small <default: 10>}

\item{dims}{Dimensions of reduction used as input <default: 1:30>}

\item{normalize}{Bool of whether to normalize using the SCT transform <default: FALSE>}

\item{dim.reduction}{Bool of whether to perform dimension reduction with PCA, tSNE and UMAP <default: FALSE>}
}
\value{
Seurat object generated from input sparse matrix.
}
\description{
Takes input single cell sparse matrix of gene counts and creates an analysis object using the Seurat R package
}
\examples{
seurat.obj <- createSeuratObjectFromMatrix(sc.data,project.name = "Analysis1",npca = 20, min.genes = 700)
}
\references{
Hao Y, Hao S, et al (2021). 
Integrated analysis of multimodal single-cell data.
\emph{Cell}. doi:10.1016/j.cell.2021.04.048
}
