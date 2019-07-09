#' @title A matrix of DNA methylation at 1000 simulated CpG sites
#'
#' @description  This simulated dataset of DNA methylation is based on 1000 randomly chosen
#' CpG sites of the Generation R methylation dataset at birth. Mean and SD were
#' extracted from the real dataset and used to simulate normal distributions
#' representing DNA methylation with some added noise.
#'
#' See extra_scripts/simulate_cpg.R for script used to generate the dataset
#'
#' @format A matrix of 1000 CpG simulated datasets for 500 participants
#' \describe{
#'   \item{row.names}{participant ID}
#'   \item{cgXXXXXXXX}{CpG ID}
#' }
#' @source \url{https://generationr.nl/}
"cpg_simulated"

#' @title A data.frame with simulated confounders and outcome
#'
#' @description This simulated dataset contains two confounders and ane outcome. These
#' variables were specified to be associated with 5 CpG sites of the
#' cpg_simulated dataset.
#'
#' See extra_scripts/simulate_covariates_outcomes.R for the script used to
#' generate the dataset
#'
#' @format A data.frame of 3 variables for 500 participants
#' \describe{
#'   \item{ID}{participant ID}
#'   \item{covariate1}{Simulated covariate 1}
#'   \item{covariate1}{Simulated covariate 2}
#'   \item{outcome}{Simulated outcome}
#' }
"phenotype"

#' @title Similarity DNA methylation matrix
#'
#' @description This similarity matrix is based on the simulated DNA
#' methylation dataset cpg_simulated
#'
#' See extra_scripts/make_similarity_matrix.R for the script used to
#' generate the similarity matrix
#'
#' @format A matrix of 500x500 participants
#' }
"Gmt"
