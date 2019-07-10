# Testing whether the simulated data works
# Load all simulated data
data("phenotype")
data("Gmt")

library(qgg) #GREML

###  Predict outcome with covariates and mehtylation matrix
# Set up covariates
X <- model.matrix(outcome ~ covariate1 + covariate2, data=phenotype)

# GREML
fitM_full <- greml(y=phenotype$outcome, X=X, GRM=list(methylation=Gmt))

# Calculate ICC
# (variance explained by methylation matrix after accounting for covariates)
fitM_full$theta["methylation"]/(fitM_full$theta["methylation"] +
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
