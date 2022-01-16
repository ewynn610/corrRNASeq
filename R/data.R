#' Example RNA-Seq data from paired observations
#'
#' A list containing 4 components: a data frame that has the sample metadata,
#' a matrix with  RNA-Seq counts (simulated from a negative binomial distribution), a matrix with transformed gene expression values,
#' and a vector denoting if each gene in the count/transformed expression matrices is differentially expressed
#' (changes in expression from baseline to follow up in the treatment group).
#'
#'
#' .
#'
#' @format A list with 4 named objects:
#' \describe{
#'   \item{metadata}{Data frame with variables for group, time, and subject ID. The \pkg{DESeq2} size factors (log scale) are also included for use as an offset term. Rows are in the same order as columns of counts.}
#'   \item{counts}{Gene expression values.}
#'   \item{vst_expr}{Transformed gene expression values using the \code{vst} function from the \pkg{DESeq2}-package}
#'   \item{de}{Vector of logical values denoting if each gene (row) in counts and vst_expr is differentially expressed.}
#' }
"simdata"
