#' Estimate variance explained by a similarity matrix above
#' variance explained by covariates
#'
#' @param outcome A string specifying the outcome variable name.
#' @param fixed_covar A string specyfying the fixed effects.
#' @param random_full A list of similarity matrices to test.
#' @param random_baseline A list of similarity matrices to adjust for.
#' @param validation_proportion The size of the validation set in percentage.
#' @param repetitions How often should cross-validation be repeated.
#' @param seed Seed to ensure reproducibility of cross-validation.
#' @return Vector of variance explained estimates from cross-validation.
#' @examples
#' # Load datasets
#' data("phenotype") # Phenotype data
#' data("Gmt")       # Methylation similarity matrix
#' data("Batch")     # Batch similarity matrix
#'
#' # Regress outcome on methylation similarity matrix, batch similarity matrix,
#' # two fixed effect covariates with GREML and Monte-Carlo Cross-Validation.
#' # Returns variance explained minus the variance explained by fixed effect
#' # covariates and similarity matrices specified in random_baseline.
#' Gmt_variance_explained <- omicsR2(outcome = outcome, fixed_covar = "covariate1 + covariate2",
#'  random_full = list(Methylation=Gmt, Batch=Batch),
#'   random_baseline = list(Batch=Batch), validation_proportion = 0.2,
#'   repetitions = 100, seed = 20190405)
#'
#' # Examine the variance explained distribution from cross-validation
#' mean(Gmt_variance_explained)
#' quantile(Gmt_variance_explained)
#' hist(Gmt_variance_explained)


omicsR2 <- function(outcome, fixed_covar, random_full, random_baseline, validation_proportion = 0.2, repetitions = 100, seed = 0) {
  # Set up fixed effects models
  fixed_covar.model <- as.formula(paste(outcome, fixed_covar, sep = " ~ "))

  # Set up fixed effects design matrix
  fixed_covar.model.matrix <- model.matrix(fixed_covar.model, data=data)

  # Set up monte-carlo cross-validation
  n <- length(data[,outcome])
  set.seed(seed)
  validate <- replicate(repetitions, sample(1:n, validation_proportion*n))

  # Fit baseline model
  baseline.fit <- qgg::greml(y=data[,outcome], X=fixed_covar.model.matrix, GRM=random_baseline, validate = validate)
  # Save cross-validated R2 for baseline model
  r2_covariates <- baseline.fit$accuracy$R2

  # Fit full model
  full.fit <- qgg::greml(y=data[,outcome], X=fixed_covar.model.matrix, GRM=random_full, validate = validate)

  # Save cross-validated R2 and substract R2 from baseline model
  r2_full <- full.fit$accuracy$R2
  r2_diff <- r2_full - r2_covariates

  # Return R2 of full model minus baseline model
  return(r2_diff)
}
