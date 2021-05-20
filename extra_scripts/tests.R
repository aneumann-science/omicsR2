# Testing whether the simulated data works
# Load all simulated data
data("phenotype")
data("Gmt")
data("Batch")

library(qgg) #GREML

# Match phenotype data and methylation similatry matrix
# First randomly rearrange the datasets
phenotype <- phenotype[sample(1:nrow(phenotype)), ]
Gmt <- Gmt[sample(1:nrow(phenotype)),]

# Match phenotype and methylation similarity matrix by ID variable
matched.list <- match_pheno_similarity(phenotype, Gmt, "ID")

# Save the data as seperate objects again
phenotype_matched <- matched.list[[1]]
Gmt_matched <- matched.list[[2]]

###  Predict outcome with covariates and mehtylation matrix
# Set up covariates
X <- model.matrix(outcome ~ covariate1 + covariate2, data=phenotype)

# GREML
fitM_full <- greml(y=phenotype$outcome, X=X, GRM=list(batch=Gmt_matched))

# Calculate ICC
# (variance explained by methylation matrix after accounting for covariates)
fitM_full$theta["methylation"]/(fitM_full$theta["methylation"] + fitM_full$theta["batch"] +
                                  fitM_full$theta["E"])

### Predict outcome with mehtylation matrix only
# Set up intercept only for covariates
X_no_covariates <- model.matrix(outcome ~ 1, data=phenotype)

# GREML
fitM_no_covariates <- greml(y=phenotype$outcome, X=X_no_covariates,
                            GRM=list(methylation=Gmt))

# Calculate ICC
# (variance explained by methylation matrix)
fitM_no_covariates$theta["methylation"]/
  (fitM_no_covariates$theta["methylation"] + fitM_no_covariates$theta["E"])
