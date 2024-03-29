% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mergeVDJresults.R
\name{mergeVDJresults}
\alias{mergeVDJresults}
\title{Merge Immcantation VDJ results}
\usage{
mergeVDJresults(df1, df2, umi.thresh = 2, assay = "bcr")
}
\arguments{
\item{df1}{Immcantation results data frame either BCR heavy chains or TCR Beta chains <default: NULL>}

\item{df2}{Immcantation results data frame either BCR light chains or TCR Alpha chains <default: NULL>}

\item{umi.thresh}{UMI count threshold for all contigs <default: 2>}

\item{assay}{what assay type is being used among only two options (bcr,tcr) <default: 'bcr'>}
}
\value{
merge BCR or TCR dataframe.
}
\description{
Takes input single cell immcantation results dataframes and merges them so that each cell has two chain assemblies.
Each BCR should have one heavy and one light chain. Each TCR should have one alpha and one beta chain.
}
\examples{
bcr_merged_df <- mergeVDJresults(df1 = bcr_heavy_df, df2 = bcr_light_df,umi.thresh = 3,assay = "bcr")
}
