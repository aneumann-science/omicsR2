# Load BGData for similarity matrix computation
library(BGData)

# Load omicsR2 package for GREML analysis
library(omicsR2)

# Load methylation and phenotype data
data("cpg_simulated") # Methylation data
data("phenotype") # Simulated phenotypes data
data("Batch")     # Batch similarity matrix

# The methylation data contains methylation levels
# for 1000 CpG probes (columns)
# and 500 participants (rows)
cpg_simulated[1:5,1:5]

# The phenotype data contains simulated outcome data and two covariates
head(phenotype, 5)

# Also a batch matrix is included, defining which participants are on the same
# batch (1)
Batch[c(1:5,100:105),c(1:5,100:105)]

# Compute similarity matrix for DNA methylation
# getG does the following:
# 1. Z-score standardizes the methylation values
# 2. Takes the product of the transpose
# 3. Scales the similarity matrix so mean diagonal is 1
Gmt <- getG(cpg_simulated)
Gmt[1:5, 1:5]

# Match phenotype and methylation similarity matrix by ID variable,
# so that participants appear in the same order
matched.list <- match_pheno_similarity(phenotype, Gmt, "ID")

# Save the data as separate objects again
phenotype_matched <- matched.list[[1]]
Gmt_matched <- matched.list[[2]]

# Regress outcome on methylation similarity matrix, batch similarity matrix,
# two fixed effect covariates with GREML and Monte-Carlo Cross-Validation.
# Returns variance explained by DNA methylation minus the variance explained by
# fixed effect covariates and similarity matrices specified in random_baseline.
# Test:validation ratio is 80:20 and randomly sampled 100 times.
Gmt_variance_explained <- omicsR2(outcome = "outcome",
                                  fixed_covar = "covariate1 + covariate2",
 random_full = list(Methylation=Gmt_matched, Batch=Batch),
  random_baseline = list(Batch=Batch),
 data = phenotype_matched,
 validation_proportion = 0.2, repetitions = 100, seed = 20190405)

# Examine the variance explained distribution obtained from cross-validation
head(Gmt_variance_explained)
mean(Gmt_variance_explained)
quantile(Gmt_variance_explained)
hist(Gmt_variance_explained, breaks = 10)
