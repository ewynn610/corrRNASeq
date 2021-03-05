#' Example RNA-Seq data from paired observations
#'
#' A list that has 2 components: A matrix with  RNA-Seq counts
#' (simulated from a negative binomial distribution) and
#' a data frame that has the sample metadata.
#'
#' @format A list with 2 named objects:
#' \describe{
#'   \item{counts}{gene expression values}
#'   \item{metadata}{data frame with variables for group, time, and subject ID.  Rows are in the same order as columns of counts}
#'   ...
#' }
"simdata"