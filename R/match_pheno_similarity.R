#' Match phenotype and similarity matrix participant order by ID variables
#'
#' @param phenotype A phenotype data.frame.
#' @param similarity A similarity matrix.
#' @param id A string for ID variable in phenotype data \code{x},
#'  which corresponds to \code{rownames} in similarity matrix \code{y}
#' @return A list of matched phenotype and similarity matrix data.
#' @examples
#' # Load datasets
#' data("phenotype") # Phenotype data
#' data("Gmt")       # Methylation similarity matrix
#'
#' # Match phenotype and methylation similarity matrix by ID variable
#' matched.list <- match_pheno_similarity(phenotype, Gmt, "ID")
#'
#' # Save the data as separate objects again
#' phenotype_matched <- matched.list[[1]]
#' Gmt_matched <- matched.list[[2]]
match_pheno_similarity <- function(phenotype, similarity, id) {
  # First match order of similarity matrix to the phenotype ID and remove
  # non-matching participants
  similarity_matched <- similarity[na.omit(match(phenotype[,id],
                                                 rownames(similarity))),
                                   na.omit(match(phenotype[,id],
                                                 colnames(similarity)))]
  # Then match phenotype data to the already matched similarity matrix
  phenotype_matched <- phenotype[na.omit(match(rownames(similarity_matched),
                                               phenotype[,id])),]
  # Perform  a check whether matching was succesfull
  print(ifelse(all(phenotype_matched[,id]==rownames(similarity_matched)),
         "phenotype and similarity matrix successfully matched :)",
         "Data not matched :("))
  # Return a list of matched phenotype and similarity matrix data
  return(list(phenotype_matched, similarity_matched))
}


