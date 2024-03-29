#' Function to fit various types of models to longitudinal RNA-Seq data
#'
#' Wrapper function that fits one of five types of models to RNA-Seq data. Available model fitting methods are linear mixed models (lmm) (using transformed data), generalized estimating equations with an optional small sample adjustment (gee), and negative binomial models using
#' a pseudo-likelhood approach (nbmm_pl), a maximum likelihood approach with Laplace approximation (nbmm_lp), or a maximum likelihood approach with adaptive Gaussian quadrature (nbmm_agq).
#'
#' @param formula A one-sided linear formula describing the model variables. For models that include them, random effects should be included in the formula using the syntax of the lme4 package.
#' @param expr_mat A (G x N) numeric matrix RNA-seq expression data with genes in rows and samples in columns. For \code{method="gee"}, \code{method="nbmm_pl"}, \code{method="nbmm_lp"} and \code{method="nbmm_agq"}, the matrix should contain raw counts and for \code{method="lmm"} the matrix should contain transformed counts (e.g. using VST from DESeq2). G = number of genes.  N = number of samples.
#' @param gene_names An optional character vector of gene names (length G).
#' @param sample_data Data frame with N rows containing the fixed- and random-effects terms included in the formula.  The rows of the data frame must correspond (and be in the same order as) the columns of the expression matrix.
#' @param method Method to use to fit the models. Possible options are \code{"lmm"}, \code{"gee"}, \code{"nbmm_pl"}, \code{"nbmm_lp"} and \code{"nbmm_agq"}.
#' @param id Only applicable for models fit using \code{method="gee"}. A vector or data column name which identifies the clusters. The length of
#' ‘id’ should be the same as the number of observations. Data are
#' assumed to be sorted so that observations on each cluster appear
#' as contiguous rows in data. If data is not sorted this way, the
#' function will not identify the clusters correctly. If \code{sort=TRUE} (default),
#' the dataframe from the \code{data} argument is sorted by the id column to avoid
#' this issue.
#' @param small.samp.method Only applicable for models fit using \code{method="gee"}. A character string specifying the
#' small sample method. The following are permitted: "pan" for the
#' Pan (2001) method, "md" for the Mancl and Derouen (2001) method, and "wl" for the Wang and Long (2011) method.
#' If \code{small.samp.method} is null, small sample variance estimates are not computed.
#' @param parallel A logical variable indicating whether forking (via mclapply) should be used to parallelize fits. Only available on Mac or linux machines.
#' @param cores Number of cores to use if parallelizing (default is 2)
#' @param ... additional arguments passed on to \code{\link[lmerTest]{lmer}} (\code{method="lmm"}), \code{\link{gee_small_sample}} (\code{method="gee"}),
#' \code{\link{glmm_nb_lmer}} (\code{method="nbmm_pl"}), \code{\link[glmmADMB]{glmmadmb}} (\code{method="nbmm_lp"}), or \code{\link[GLMMadaptive]{mixed_model}} (method="nbmm_agq").
#'
#' @return A list of length G of model objects from the following functions: \code{\link[lmerTest]{lmer}} (\code{method="lmm"}), \code{\link{gee_small_sample}}(\code{method="lmm"}),
#' \code{\link{glmm_nb_lmer}} (\code{method="nbmm_pl"}), \code{\link[glmmADMB]{glmmadmb}} (\code{method="nbmm_lp"}), or \code{\link[GLMMadaptive]{mixed_model}} (method="nbmm_agq").
#'
#'
#' @author Elizabeth Wynn
#'
#' @seealso \code{\link{corrSeq_summary}}, \code{\link{glmm_nb_lmer}}, \code{\link[lmerTest]{lmer}}, \code{\link{gee_small_sample}}, \code{\link[glmmADMB]{glmmadmb}}, and \code{\link[GLMMadaptive]{mixed_model}}
#'
#' @examples
#' data("simdata")
#' sample_meta_data <- simdata$metadata
#'
#' ## Subset down to 10 observation (i.e. gene)
#' counts=simdata$counts[1:10,]
#'
#' ## Fit GEE models using Wang-Long small sample size estimator
#' ## log_offset is the log size factors from the DESeq2 package
#' gee_fit <- corrSeq_fit(formula = ~ group * time+offset(log_offset),
#'                            expr_mat = counts,
#'                            sample_data = sample_meta_data,
#'                            method="gee",
#'                            id=ids,
#'                            small.samp.method="wl")
#'
#' ## Fit NBMM-PL models
#' nbmm_pl_fit <- corrSeq_fit(formula = ~ group * time+(1|ids)+offset(log_offset),
#'                            expr_mat = counts,
#'                            sample_data = sample_meta_data,
#'                            method="nbmm_pl")
#'
#'
#' ## Fit NBMM-AGQ models
#' nbmm_agq_fit <- corrSeq_fit(formula = ~ group * time+(1|ids)+offset(log_offset),
#'                            expr_mat = counts,
#'                            sample_data = sample_meta_data,
#'                            method="nbmm_agq")
#'
#' ## Fit NBMM-LP models
#' ## Random effects must be factors
#' sample_meta_data$ids<-factor(sample_meta_data$ids)
#' nbmm_lp_fit <- corrSeq_fit(formula = ~ group * time+(1|ids)+offset(log_offset),
#'                            expr_mat = counts,
#'                            sample_data = sample_meta_data,
#'                            method="nbmm_lp")
#'
#' ## Fit LMM models to transformed data
#' ## Use the variance transformed counts in the simdata object
#' ## Subset down to 10 genes
#' vst_expr<-simdata$vst_expr[1:10,]
#'
#' ##  Fit the Models
#' lmm_fit<- corrSeq_fit(formula = ~ group * time + (1|ids),
#'                            expr_mat = vst_expr,
#'                            sample_data = sample_meta_data,
#'                            method="lmm")
#' @export
#'

corrSeq_fit <- function(formula = NULL, # Formula for fixed effects
                        expr_mat = NULL, # Matrix of transformed RNA-Seq counts where rows are genes and columns are samples
                        gene_names = NULL, # A vector of gene names (the length of the number of rows in the expression matrix).  If unspecified, rownames from the expression matrix will be used.
                        sample_data = NULL, # A data frame with sample meta data
                        method,
                        id,
                        small.samp.method=NULL,
                        parallel = F,
                        cores = 2,
                        ...
){
  if(method=="lmm"){
    # Gene Names
    if(is.null(gene_names)){
      if(is.null(rownames(expr_mat))==T){rownames(expr_mat)<-paste("gene_", seq(1,nrow(expr_mat),1))}
      gene_names =  rownames(expr_mat)
    }
    ret=lmerSeq::lmerSeq.fit(form = formula, expr_mat = expr_mat,
                             gene_names = gene_names,
                             sample_data = sample_data, parallel = parallel,
                             cores = cores,...)
  }else{
    ############################################################################################################
    #Error Messages for insufficient or inconsistent information
    ############################################################################################################

    ### Insufficient Information ###

    if(is.null(expr_mat)==T ) {
      stop("An expression matrix must be provided.")}

    if(is.null(formula)==T ) {
      stop("A formula must be provided.")}

    if(is.null(sample_data)==T ) {
      stop("sample_data is missing.")}

    ### Inconsistent information ###

    method_choices=c("gee", "nbmm_pl", "nbmm_lp", "nbmm_agq", "lmm")

    if(!(method%in%method_choices) ) {
      stop("Invalid method")}

    if((ncol(expr_mat)==nrow(sample_data))==F ) {
      stop("The expression matrix and sample data include differing numbers of samples.")}

    if(is.null(gene_names)==F & (nrow(expr_mat)==length(gene_names))==F ) {
      print("The expression matrix and gene_names indicate differing numbers of genes.
          Row names of the expression matrix will be used as gene names.")
      gene_names = rownames(expr_mat)
    }

    ################################################################################################
    # Calculate Default Values if none supplied
    ################################################################################################

    # Gene Names
    if(is.null(gene_names)){
      if(is.null(rownames(expr_mat))==T){rownames(expr_mat)<-paste0("gene_", seq(1,nrow(expr_mat),1))}
      gene_names =  rownames(expr_mat)
    }

    ############################################################################################################
    # Begin Analysis
    ############################################################################################################
    # Make sure expr_mat is a matrix
    expr_mat <- as.matrix(expr_mat)

    # Ensure that contrast, gene and fixed effect names are supplied as characters
    gene_names <- as.character(gene_names)
    form_sub <- update(formula, expr ~ .)

    if(method=="gee"){

      args<-list(formula=form_sub, data=quote(dat_sub), id=substitute(id), small.samp.method=small.samp.method, ...)

      method_call=geeglm_small_samp
    }else if(method=="nbmm_lp"){
      ## Make sure random effects are factor
      random=paste("(", lme4::findbars(form_sub), ")")
      random<-sapply(random, function(x) gsub(" ", "", stringr::str_split(stringr::str_split(x, pattern = "\\)")[[1]], "\\|")[[1]][2]))
      if(sum(sapply(random, function(x) class(sample_data[,x])!="factor"))>0) stop("All grouping variables in random effects must be factors")

      args<-list(formula=form_sub, data=quote(dat_sub), #random=random,
                 family = "nbinom",...)
      method_call=glmmADMB::glmmadmb

    }else if(method=="nbmm_pl"){
      args<-list(formula=form_sub, data=quote(dat_sub),...)
      method_call=glmm_nb_lmer
    }else if(method=="nbmm_agq"){
      ## Separate random and fixed effects
      random=as.formula(paste0("~",paste0(lme4::findbars(form_sub), collapse = "+")))
      fixed_eff=lme4::nobars(form_sub)
      args=args2=list(fixed=fixed_eff, random=random, data=quote(dat_sub),
                      family=GLMMadaptive::negative.binomial(), ...)
      args2$family=poisson()
      method_call=GLMMadaptive::mixed_model
    }

    if(parallel == F){
      ret <- pbapply::pblapply(X = 1:nrow(expr_mat),
                               FUN = function(i){
                                 dat_sub <- cbind(sample_data, data.frame(expr = as.numeric(expr_mat[i, ])))
                                 ret_sub <- tryCatch({
                                   tmp1 <- suppressWarnings(suppressMessages(do.call(method_call, args)))
                                 }, error = function(e) {
                                   if(method=="nbmm_agq"){
                                     #If error, try fitting with a poisson model
                                     ret_sub2  <- tryCatch({
                                       tmp1 <- suppressWarnings(suppressMessages(do.call(method_call, args2)))
                                     }, error=function(e){NULL})
                                   }else ret_sub2 <- NULL
                                   ret_sub2
                                 })
                                 ret2 <- ret_sub
                                 if(method=="nbmm_lp"&!is.null(ret2)) ret2$data=dat_sub
                                 ret2
                               })
    }
    else{
      ret <- parallel::mclapply(X = 1:nrow(expr_mat),
                                mc.silent = F,
                                mc.cores = cores,
                                FUN = function(i){
                                  dat_sub <- cbind(sample_data, data.frame(expr = as.numeric(expr_mat[i, ])))
                                  ret_sub <- tryCatch({
                                    tmp1 <- suppressWarnings(suppressMessages(do.call(method_call, args)))
                                  }, error = function(e) {
                                    if(method=="nbmm_agq"){
                                      #If error, try fitting with a poisson model
                                      ret_sub2  <- tryCatch({
                                        tmp1 <- suppressWarnings(suppressMessages(do.call(method_call, args2)))
                                      }, error=function(e){NULL})
                                    }else ret_sub2 <- NULL
                                    ret_sub2
                                  })
                                  ret2 <- ret_sub
                                  if(method=="nbmm_lp"&!is.null(ret2)) ret2$data=dat_sub
                                  ret2
                                })
    }
  }
  names(ret)<-gene_names
  return(ret)
}
