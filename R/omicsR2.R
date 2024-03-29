#' Estimate variance explained by a similarity matrix above
#' variance explained by covariates
#'
#' @param outcome A string specifying the outcome variable name.
#' @param fixed_covar A string specifying the fixed effects.
#' @param random_full A list of similarity matrices to test.
#' @param random_baseline A list of similarity matrices to adjust for.
#' @param data A data.frame containing outcome and covariate data.
#' @param validation_proportion The size of the validation set in percentage.
#' @param repetitions How often should cross-validation be repeated.
#' @param ncores Number of cores.
#' @param seed Seed to ensure reproducibility of cross-validation.
#' @return Data.frame with variance explained estimates from cross-validation
#' for baseline (r2_covariates) and full model (r2_full), plus the difference
#' (r2_diff).
#' @examples
#' # Load datasets
#' data("phenotype") # Phenotype data
#' data("Gmt")       # Methylation similarity matrix
#' data("Batch")     # Batch similarity matrix
#'
#' # Regress outcome on methylation similarity matrix, batch similarity matrix,
#' # two fixed effect covariates with GREML and Monte-Carlo Cross-Validation.
#' # Returns variance explained by DNA methylation minus the variance explained by
#' # fixed effect covariates and similarity matrices specified in random_baseline.
#' # Training:validation ratio is 80:20 and randomly sampled 100 times.
#' Gmt_variance_explained <- omicsR2(outcome = "outcome",
#'                                   fixed_covar = "covariate1 + covariate2",
#'                                   random_full = list(Methylation=Gmt, Batch=Batch),
#'                                   random_baseline = list(Batch=Batch),
#'                                   data = phenotype,
#'                                   validation_proportion = 0.2, repetitions = 100, seed = 20190405)
#'
#' # Examine the variance explained distribution obtained from cross-validation
#' head(Gmt_variance_explained$r2_diff)
#' mean(Gmt_variance_explained$r2_diff)
#' quantile(Gmt_variance_explained$r2_diff)
#' hist(Gmt_variance_explained$r2_diff, breaks = 10)

omicsR2 <- function(outcome, fixed_covar, random_full, random_baseline, data, validation_proportion = 0.2, repetitions = 100, ncores = 1, seed = 0) {
  # Set up fixed effects models
  fixed_covar.model <- as.formula(paste(outcome, fixed_covar, sep = " ~ "))

  # Set up fixed effects design matrix
  fixed_covar.model.matrix <- model.matrix(fixed_covar.model, data=data)

  # Set up monte-carlo cross-validation
  n <- length(data[,outcome])
  set.seed(seed)
  validate <- replicate(repetitions, sample(1:n, validation_proportion*n))

  # Fit baseline model
  baseline.fit <- qgg::greml(y=data[,outcome], X=fixed_covar.model.matrix, GRM=random_baseline, validate = validate, ncores = ncores)
  # Save cross-validated R2 for baseline model
  r2_covariates <- baseline.fit$accuracy$R2

  # Fit full model
  full.fit <- qgg::greml(y=data[,outcome], X=fixed_covar.model.matrix, GRM=random_full, validate = validate, ncores = ncores)

  # Save cross-validated R2 and substract R2 from baseline model
  r2_full <- full.fit$accuracy$R2
  r2_diff <- r2_full - r2_covariates

  # Collect all R2 into one data.frame
  r2.data <- data.frame(r2_covariates, r2_full, r2_diff)

  # Return R2 of baseline and full model, plus the change in R2
  return(r2.data)
}
