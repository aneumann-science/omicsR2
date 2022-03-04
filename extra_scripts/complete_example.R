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
# Training:validation ratio is 80:20 and randomly sampled 100 times.
Gmt_variance_explained <- omicsR2(outcome = "outcome",
                                  fixed_covar = "covariate1 + covariate2",
 random_full = list(Methylation=Gmt_matched, Batch=Batch),
  random_baseline = list(Batch=Batch),
 data = phenotype_matched,
 validation_proportion = 0.2, repetitions = 100, seed = 20190405)

# Examine the variance explained distribution obtained from cross-validation
head(Gmt_variance_explained$r2_diff)
mean(Gmt_variance_explained$r2_diff)
quantile(Gmt_variance_explained$r2_diff)
hist(Gmt_variance_explained$r2_diff, breaks = 10)

# Multiple Imputation
# Code adapted from semTools runMI documentation

# Randomly set covariates to missing
set.seed(20190405)
phenotype_matched_na <- phenotype_matched
phenotype_matched_na$covariate1 <- ifelse(phenotype_matched_na$covariate1 <= quantile(phenotype_matched_na$covariate1, .3), NA, phenotype_matched_na$covariate1)
phenotype_matched_na$covariate2 <- ifelse(phenotype_matched_na$covariate2 <= quantile(phenotype_matched_na$covariate2, .3), NA, phenotype_matched_na$covariate2)

# Impute data five times
library(mice)
miceImps <- mice(phenotype_matched_na)

# Create a list of imputed datasets
impList <- list()
for (i in 1:miceImps$m) impList[[i]] <- complete(miceImps, action = i)

# Repeat GREML analysis, but this time on multiply imputed data
Gmt_variance_explained_imp <- omicsR2_imp(outcome = "outcome",
                                      fixed_covar = "covariate1 + covariate2",
                                      random_full = list(Methylation=Gmt_matched, Batch=Batch),
                                      random_baseline = list(Batch=Batch),
                                      imp_list = impList,
                                      validation_proportion = 0.2, repetitions = 100, seed = 20190405)

# Examine the variance explained distribution across cross-validations and
# imputed datasets
all_estimates.data <- do.call(rbind, Gmt_variance_explained_imp[[1]])
hist(all_estimates.data$r2_diff, breaks = 10)
Gmt_variance_explained_imp[[3]]
