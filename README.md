# omicsR2

Method to investigate how much variance of a given trait is explained through DNA methylation as a whole, centered around the `greml` function described in the `qgg` R package

It is possible to define covariates (or genetic similarity) that are taking into account before estimating the variance explained through DNA methylation.


## Installation

         library(devtools)
         install_github("aneumann-science/omicsR2")

## Step-by-step guide

You will need at least:

* `meth`: your DNA methylation data (e.g. Illumina 450k or EPIC array)
* `pheno`: your phenotype file with your trait of interest (e.g. BMI) and covariates (e.g. cell type, batch), for all participants (and in the same order) present in your DNA methylation file

optional data:

* genome-wide data

In the first step, you will create similarity matrixes based on your DNA methylation (and, if required, genetic) data. This step is centered around the `getG` function in the `BGData` package to compute a genomic relationship matrix.

         Gmt <- getG(meth)

In the second step, you will fit your fit your greml model. You first need to define your model:

         X <- model.matrix(trait.of.interest ~ covariates, data=pheno)

You then fit greml:

         fitM <- greml(y=trait.of.interest, X=X, GRM=list(M=Gmt))

You can include several relationship matrices (e.g. genetic or batch), using `GRM=list(M=Gmt, G=..., Batch=...)`

Variance explained (R^2) is defined as the `theta$M / (theta$M [+ theta$G + ...] + theta$E)` with theta being the covariance estimates from REML.
Crossvalidated estimates are computed using 20% of data for testing and a resampling of n=100. 

 
 
 
