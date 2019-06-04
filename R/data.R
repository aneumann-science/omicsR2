#' A matrix of DNA methylation at 1000 simulated CpG sites
#'
#' This simulated dataset of DNA methylation is based on 1000 randomly chosen
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
