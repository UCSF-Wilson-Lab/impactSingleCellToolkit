% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/addSingleRcelltypeAnnotToSeuratObject.R
\name{addSingleRcelltypeAnnotToSeuratObject}
\alias{addSingleRcelltypeAnnotToSeuratObject}
\title{Add SingleR cell-type annotations into Seurat Object}
\usage{
addSingleRcelltypeAnnotToSeuratObject(
  seurat,
  singler_results_list,
  cluster_col = "seurat_clusters"
)
}
\arguments{
\item{seurat}{Seurat Object}

\item{singler_results_list}{list of SingleR results objects}

\item{cluster_col}{Column name with the Seurat Object metadata that contains the cluster IDs <default: 'seurat_clusters'>}
}
\value{
Seurat Object with SingleR annotations.
}
\description{
Takes an input Seurat Object and adds inputted SingleR cell-type annotations into the metadata of the object
}
\examples{
seurat <- addSingleRcelltypeAnnotToSeuratObject(seurat,singler_results_list)
}
