# impactSingleCellToolkit

This R package is designed to make use of the ever expanding pool of single cell analysis tools. This package is particularly focused on the analysis of human scRNA-Seq, VDJ and CITE-Seq data sets. 

## Installation

Use [devtools](https://github.com/hadley/devtools "devtools") to install this package. You will first need to install all dependencies.

```{r}
install.packages("devtools")
```

```r
library(devtools)
install_github("UCSF-Wilson-Lab/impactSingleCellToolkit")
```

## Dependencies

This package requires pre-installation of the following R packages:

* [Seurat](https://github.com/satijalab/seurat "Seurat")
* [SingleCellExperiment](https://bioconductor.org/packages/release/bioc/html/SingleCellExperiment.html "SingleCellExperiment")
* [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder "DoubletFinder")
* [SingleR](https://github.com/LTLA/SingleR "SingleR")
* [celldex](http://bioconductor.org/packages/release/data/experiment/html/celldex.html "celldex")
* [Biostrings](https://bioconductor.org/packages/release/bioc/html/Biostrings.html "Biostrings")

