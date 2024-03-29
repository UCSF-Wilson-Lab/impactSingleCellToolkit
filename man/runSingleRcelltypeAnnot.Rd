% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/runSingleRcelltypeAnnot.R
\name{runSingleRcelltypeAnnot}
\alias{runSingleRcelltypeAnnot}
\title{Generate cell-type annotations using SingleR}
\usage{
runSingleRcelltypeAnnot(
  seurat_obj,
  blueprintencode.annot = TRUE,
  monaco.annot = TRUE
)
}
\arguments{
\item{seurat_obj}{Seurat Object}

\item{blueprintencode.annot}{Bool of whether to generate celltype annotations using the Blueprint ENCODE dataset as a reference <default: TRUE>}

\item{monaco.annot}{Bool of whether to generate celltype annotations using the Monaco dataset as a reference <default: TRUE>}
}
\value{
SingleR results object.
}
\description{
Takes an input Seurat Object and generates cell-type annotations for human data single cell RNA-Seq data using the SingleR R package
}
\examples{
singler_results <- runSingleRcelltypeAnnot(seurat_obj)
}
\references{
Aran D, Looney AP, Liu L et al. (2019).
Reference-based analysis of lung single-cell sequencing reveals a transitional profibrotic macrophage.
\emph{Nat. Immunology} 20, 163–172.
}
