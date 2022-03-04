#' Estimate variance explained by a similarity matrix above
#' variance explained by covariates for multiply imputed data. This function
#' runs \code{\link{omicsR2}} for each imputed datasets. Results
#' are then pooled using Rubin's rule.
#'
#' @param outcome A string specifying the outcome variable name.
#' @param fixed_covar A string specifying the fixed effects.
#' @param random_full A list of similarity matrices to test.
#' @param random_baseline A list of similarity matrices to adjust for.
#' @param imp_list A list of imputed datasets.
#' @param validation_proportion The size of the validation set in percentage.
#' @param repetitions How often should cross-validation be repeated.
#' @param ncores_imp number of cores used to run imputed datasets in parallel.
#' @param ncores Number of cores to run cross-validations in parallel.
#' @param seed Seed to ensure reproducibility of cross-validation.
#' @return List with three data.frames:
#' 1. Per imputed dataset: data.frame with variance explained estimates
#' per cross-validation for baseline (r2_covariates) and full model (r2_full),
#' plus the difference (r2_diff).
#' 2. Per imputed dataset, the mean and sd of the r2 estimates across
#' cross-validations.
#' 3. The pooled estimates across imputations using Rubin's rule.
#' @examples
#' # Load datasets
#' data("phenotype") # Phenotype data
#' data("Gmt")       # Methylation similarity matrix
#' data("Batch")     # Batch similarity matrix
#'
#' # Randomly set covariates to missing
#' set.seed(20190405)
#' phenotype_matched_na <- phenotype
#' phenotype_matched_na$covariate1 <- ifelse(phenotype_matched_na$covariate1 <= quantile(phenotype_matched_na$covariate1, .3), NA, phenotype_matched_na$covariate1)
#' phenotype_matched_na$covariate2 <- ifelse(phenotype_matched_na$covariate2 <= quantile(phenotype_matched_na$covariate2, .3), NA, phenotype_matched_na$covariate2)
#'
#' # Impute data five times
#' miceImps <- mice(phenotype_matched_na)
#'
#' # Create a list of imputed datasets
#' impList <- list()
#' for (i in 1:miceImps$m) impList[[i]] <- complete(miceImps, action = i)
#'
#' # Regress outcome on methylation similarity matrix, batch similarity matrix,
#' # two fixed effect covariates with GREML and Monte-Carlo Cross-Validation.
#' # Returns variance explained by DNA methylation minus the variance explained by
#' # fixed effect covariates and similarity matrices specified in random_baseline.
#' # Training:validation ratio is 80:20 and randomly sampled 100 times.
#' # Results are pooled across the list of imputed datasets.
#' Gmt_variance_explained <- omicsR2_imp(outcome = "outcome",
#'                                   fixed_covar = "covariate1 + covariate2",
#'                                   random_full = list(Methylation=Gmt_matched, Batch=Batch),
#'                                   random_baseline = list(Batch=Batch),
#'                                   imp_list = impList,
#'                                   validation_proportion = 0.2, repetitions = 100, seed = 20190405)
#'
#' # Examine the variance explained distribution across cross-validations and
#' # imputed datasets
#' all_estimates.data <- do.call(rbind, Gmt_variance_explained_imp[[1]])
#' hist(all_estimates.data$r2_diff, breaks = 10)
#' Gmt_variance_explained_imp[[3]]
omicsR2_imp <- function(outcome, fixed_covar, random_full, random_baseline,
                    imp_list, validation_proportion = 0.2, repetitions = 100,
                    ncores_imp = 1, ncores = 1, seed = 0) {
# Iterate over each imputed dataset and apply cross-validation, to estimate R2
# for each imputed dataset.
  R2.list <- pbmcapply::pbmclapply(imp_list, function(imp.data) {
    omicsR2(outcome=outcome, fixed_covar=fixed_covar, random_full=random_full,
            random_baseline=random_baseline, data=imp.data,
            validation_proportion = validation_proportion,
            repetitions = repetitions, ncores = ncores, seed = seed)
  }, mc.cores = ncores_imp)

  # Effect estimates per imputed dataset
  results_per_imp.list <- lapply(R2.list, function(cv) {
    r2_covariates_mean <- mean(cv$r2_covariates)
    r2_covariates_sd <- sd(cv$r2_covariates)
    r2_full_mean <- mean(cv$r2_full)
    r2_full_sd <- sd(cv$r2_full)
    r2_diff_mean <- mean(cv$r2_diff)
    r2_diff_sd <- sd(cv$r2_diff)
    r2_mean.data <- data.frame(r2_covariates_mean,r2_covariates_sd,r2_full_mean,r2_full_sd,r2_diff_mean,r2_diff_sd)
    return(r2_mean.data)
  })
  results_per_imp.data <- do.call(rbind, results_per_imp.list)

  # Pooled results
  # Pooled estimate
  r2_covariates_mean_pooled <- mean(results_per_imp.data$r2_covariates_mean)
  r2_full_mean_pooled <- mean(results_per_imp.data$r2_full_mean)
  r2_diff_mean_pooled <- mean(results_per_imp.data$r2_diff_mean)

  # Within imputation variance
  r2_covariates_within_variance <- mean(results_per_imp.data$r2_covariates_sd^2)
  r2_full_within_variance <- mean(results_per_imp.data$r2_full_sd^2)
  r2_diff_within_variance <- mean(results_per_imp.data$r2_diff_sd^2)

  # Between imputation variance
  # Number of imputations
  m <- length(results_per_imp.list)
  # Vb
  r2_covariates_between_variance <- sum((results_per_imp.data$r2_covariates_mean - r2_covariates_mean_pooled)^2)/(m-1)
  r2_full_between_variance <- sum((results_per_imp.data$r2_full_mean - r2_full_mean_pooled)^2)/(m-1)
  r2_diff_between_variance <- sum((results_per_imp.data$r2_diff_mean - r2_diff_mean_pooled)^2)/(m-1)

  # Total variance
  r2_covariates_total_variance <- r2_covariates_within_variance + r2_covariates_between_variance + r2_covariates_between_variance/3
  r2_full_total_variance <- r2_full_within_variance + r2_full_between_variance + r2_full_between_variance/3
  r2_diff_total_variance <- r2_diff_within_variance + r2_diff_between_variance + r2_diff_between_variance/3

  # Total SD
  r2_covariates_total_sd <- sqrt(r2_covariates_total_variance)
  r2_full_total_sd <- sqrt(r2_full_total_variance)
  r2_diff_total_sd <- sqrt(r2_diff_total_variance)

  # Combine all pooled results
  pooled_results.data <- data.frame(r2_covariates_mean_pooled,r2_covariates_total_sd,r2_full_mean_pooled,r2_full_total_sd,r2_diff_mean_pooled,r2_diff_total_sd)

  # Aggregate and return all results as list
  results.list <- list(R2.list,results_per_imp.data,pooled_results.data)
  return(results.list)
}
