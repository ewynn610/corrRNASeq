#' Fit Generalize Estimating Equations with Small Sample Variance Estimators
#'
#'This function fits a GEE model using the \code{geeglm} function
#'from the \pkg{geepack}-package. Additionally, small sample variance
#'estimates are calculated using one of three methods proposed in the following
#'papers: Pan (2001), Mancl and Derouen (2001), and Wang and Long (2011).
#'
#'This function borrows heavily from the corresponding small sample variance estimating
#'functions in the \pkg{geesmv}-package
#'(\code{\link[geesmv]{GEE.var.pan}}, \code{\link[geesmv]{GEE.var.md}} and \code{\link[geesmv]{GEE.var.wl}}).
#'In addition to combining these functions into a single function and using the \code{geeglm} function
#'in model fitting, this function also varys from these functions in that it
#'ensures that model offsets are properly accounted for in the calculation of
#'small sample estimators.
#'
#' @param formula A two-sided linear formula object with the response
#'  on the left of a ~ operator and the terms, separated by + operators, on the right.
#' @param family a description of the error distribution and link
#' function to be used in the model. Only \code{guassian} and \code{poisson}
#' families are supported.
#' @param data 	an optional data frame, list or environment (or object coercible by as.data.frame to a data frame)
#' containing the variables in the model and the id variable. If not
#' environment(formula), typically the environment from which the function is called.
#' @param id a vector or data column name which identifies the clusters. The length of
#' id should be the same as the number of observations. Data are
#' assumed to be sorted so that observations on each cluster appear
#' as contiguous rows in data. If data is not sorted this way, the
#' function will not identify the clusters correctly. If \code{sort=TRUE} (default),
#' the dataframe from the \code{data} argument is sorted by the id column to avoid
#' this issue.
#' @param corstr a character string specifying the correlation structure. The following are permitted: "independence", "exchangeable", "ar1", and "unstructured"
#' @param small.samp.method a character string specifying the
#' small sample method. The following are permitted: \code{"pan"} for the
#' Pan (2001) method, \code{"md"} for the Mancl and Derouen (2001) method, and \code{"wl"} for the Wang and Long (2011) method.
#' If \code{small.samp.method} is null, small sample variance estimates are not computed.
#' The resulting object will be identical to the object created if
#' \code{geeglm} from the \pkg{geepack}-package was used.
#' @param ... additional arguments passed on to \code{geepack::geeglm}.
#'
#' @import geepack geesmv
#' @return This function returns a \code{geeglm} object with one additional items:
#' \item{small.samp.var}{small sample variance estimators using the specified method}
#'
#'
#'@author Elizabeth Wynn, \pkg{geesmv}-authors for underlying code used in small sample size variance estimators.
#'
#' @seealso \code{\link[geepack]{geeglm}}, \code{\link[geesmv]{GEE.var.pan}}, \code{\link[geesmv]{GEE.var.md}} and \code{\link[geesmv]{GEE.var.wl}}
#' @references
#' Mancl LA, DeRouen TA (2001). "A covariance estimator for GEE with improved small-sample properties." Biometrics, 57(1), 126-134. ISSN 0006341X. doi:10. 1111/j.0006-341X.2001.00126.x.
#'
#' Pan W (2001). "On the robust variance estimator in generalised estimating equations." \emph{Biometrika}, 88(3), 901-906. ISSN 00063444. doi:10.1093/biomet/88.3.901.
#'
#' Wang M, Long Q (2011). "Modified robust variance estimator for generalized estimating equations with improved small-sample performance." \emph{Statistics in Medicine}, 30(11), 1278-1291. ISSN 02776715. doi:10.1002/sim.4150
#'
#'
#'
#' @examples
#' data("simdata")
#' sample_meta_data <- simdata$metadata
#'
#' #Subset down to one observation (i.e. gene)
#' counts=simdata$counts[1,]
#'
#' #Combine counts, metadata into dataframe
#' df=cbind(counts, sample_meta_data)
#'
#' #Sort data by id (Function also does this if sort=T)
#' df=df[order(df$ids),]
#'
#' #Fit the Model-use Pan method for small sample variance
#' fit.gee.pan<-geeglm_small_samp(formula =counts ~ group * time,
#'                               family=poisson, data=df, id=ids,
#'                               corstr="exchangeable",
#'                               small.samp.method="pan", sort=T)
#' @export
#'

#NOTE: Added parameter mu throughout all functions because function doesn't account for offset
#Use geepack instead of gee model
#Doesn't support binomial model

geeglm_small_samp<-function (formula,
                             family = poisson,
                             data,
                             id,
                             corstr = "exchangeable",
                             small.samp.method=NULL,
                             sort=T,
                             ...)
{
  call <- match.call(expand.dots = TRUE)
  if (is.null(data$id)) {
    index <- which(names(data) == id)
    data$id <- data[, index]
  }
  init <- model.frame(formula, data)
  init$num <- 1:length(init[, 1])
  if (any(is.na(init))) {
    index <- na.omit(init)$num
    data <- data[index, ]
    m <- model.frame(formula, data)
    mt <- attr(m, "terms")
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  else {
    m <- model.frame(formula, data)
    mt <- attr(m, "terms")
    data$response <- model.response(m, "numeric")
    mat <- as.data.frame(model.matrix(formula, m))
  }
  #################Make sure data is sorted################
  if(sort) {
    data<-data[order(data[[paste(call$id)]]),]
    call$data=quote(data)
  }

  ###Use geepack model (geesmv fits "gee" package model)###
  #Calling geeglm
  call_list<-as.list(call)
  if(is.null(call_list$family)) call_list$family=quote(poisson)
  if(is.null(call_list$corstr)) call_list$corstr="exchangeable"
  call_list[[1]]<-call_list$small.samp.method<-call_list$sort<-NULL
  gee.fit<-do.call("geeglm", call_list)

  #If small samp method is provided
  if(!is.null(small.samp.method)){
  beta_est <- gee.fit$coefficients
  alpha <- summary(gee.fit)$corr[[1]]

  #########Drawing largely from code from geesmv package############
  len <- length(beta_est)
  len_vec <- len^2
  data$id <- as.numeric(gee.fit$id)
  cluster <- cluster.size(data$id)
  ncluster <- max(cluster$n)
  size <- cluster$m
  mat$subj <- rep(unique(data$id), cluster$n)
  if (is.character(corstr)) {
    var <- switch(corstr, independence = cormax.ind(ncluster),
                  exchangeable = cormax.exch(ncluster, alpha), `AR-M` = cormax.ar1(ncluster,
                                                                                   alpha), unstructured = summary(gee.fit)$working.correlation
    )
  }
  else {
    print(corstr)
    stop("'working correlation structure' not recognized")
  }
  if (is.character(family)) {
    family <- switch(family, gaussian = "gaussian",
                     poisson = "poisson")
  }
  else {
    if (is.function(family)) {
      family <- family()[[1]]
    }
    else {
      print(family)
      stop("'family' not recognized")
    }
  }
  if(small.samp.method=="wl"){
    m <- model.frame(formula, data)
    mat <- as.data.frame(model.matrix(formula, m))
    mat$subj <- rep(unique(data$id), cluster$n)
  }
  cov.beta <- unstr <- matrix(0, nrow = len, ncol = len)

  if(small.samp.method=="pan") step <- matrix(0, nrow = cluster$n[1], ncol = cluster$n[1])
  else if(small.samp.method=="md") step11 <- matrix(0, nrow = len, ncol = len)
  else if(small.samp.method=="wl")step01 <- matrix(0, nrow = len, ncol = len)
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                  mat$subj == unique(data$id)[i]))
    #EAW added fitted values so don't have to solve using link
    mu=as.matrix(gee.fit$fitted.values[data$id == unique(data$id)[i]])
    if(small.samp.method=="pan"){
      if (family == "gaussian") {
        resid <- (y - mu) %*% t(y - mu)
        step <- step + resid
      }
      else if (family == "poisson") {
        resid <- (y - mu) %*% t(y -
                                  mu)
        B <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
        diag(B) <- 1/sqrt(mu)
        step <- step + B %*% resid %*% B
      }
    }else if(small.samp.method=="md"){
      var_i = var[1:cluster$n[i], 1:cluster$n[i]]
      if (family == "gaussian") {
        xx <- t(covariate) %*% solve(var_i) %*% covariate
        step11 <- step11 + xx
      }
      else if (family == "poisson") {
        D <- mat.prod(covariate, mu)
        Vi <- diag(sqrt(c(mu)),
                   cluster$n[i]) %*% var_i %*% diag(sqrt(c(mu)), cluster$n[i])
        xx <- t(D) %*% solve(Vi) %*% D
        step11 <- step11 + xx
      }
    }else if(small.samp.method=="wl"){
      var_i = var[1:cluster$n[i], 1:cluster$n[i]]
      if (family == "gaussian") {
        xx <- t(covariate) %*% solve(var_i) %*% covariate
        step01 <- step01 + xx
      }
      else if (family == "poisson") {
        D <- mat.prod(covariate, mu)
        Vi <- diag(sqrt(c(mu)),
                   cluster$n[i]) %*% var_i %*% diag(sqrt(c(mu)), cluster$n[i])
        xx <- t(D) %*% solve(Vi) %*% D
        step01 <- step01 + xx
      }
    }

  }
  if(small.samp.method=="wl"){
    step <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
    for (i in 1:size) {
      y <- as.matrix(data$response[data$id == unique(data$id)[i]])
      covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                    mat$subj == unique(data$id)[i]))
      #EAW added fitted values so don't have to solve using link
      mu=as.matrix(gee.fit$fitted.values[data$id == unique(data$id)[i]])
      var_i = var[1:cluster$n[i], 1:cluster$n[i]]
      if (family == "gaussian") {
        resid <- solve(cormax.ind(cluster$n[i]) - covariate %*%
                         solve(step01) %*% t(covariate) %*% solve(var_i)) %*%
          (y - mu)
        step <- step + resid %*% t(resid)
      }
      else if (family == "poisson") {
        B <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
        diag(B) <- 1/sqrt(mu)
        D <- mat.prod(covariate, mu)
        Vi <- diag(sqrt(c(mu)),
                   cluster$n[i]) %*% var_i %*% diag(sqrt(c(mu)), cluster$n[i])
        resid <- B %*% solve(cormax.ind(cluster$n[i]) - D %*%
                               solve(step01) %*% t(D) %*% solve(Vi)) %*% (y -mu)
        step <- step + resid %*% t(resid)
      }

    }
  }
  if(small.samp.method%in% c("pan", "wl")){
    unstr <- step/size
    step11 <- matrix(0, nrow = len, ncol = len)
  }
  step12 <- matrix(0, nrow = len, ncol = len)
  step13 <- matrix(0, nrow = len_vec, ncol = 1)
  step14 <- matrix(0, nrow = len_vec, ncol = len_vec)
  p <- matrix(0, nrow = len_vec, ncol = size)
  for (i in 1:size) {
    y <- as.matrix(data$response[data$id == unique(data$id)[i]])
    covariate <- as.matrix(subset(mat[, -length(mat[1, ])],
                                  mat$subj == unique(data$id)[i]))
    #EAW added fitted values so don't have to solve using link
    mu=as.matrix(gee.fit$fitted.values[data$id == unique(data$id)[i]])
    var_i = var[1:cluster$n[i], 1:cluster$n[i]]
    if(small.samp.method=="pan"){
      if (family == "gaussian") {
        xy <- t(covariate) %*% solve(var_i) %*% unstr %*%
          solve(var) %*% covariate
        xx <- t(covariate) %*% solve(var_i) %*% covariate
        step11 <- step11 + xx
        step12 <- step12 + xy
        step13 <- step13 + matrixcalc::vec(xy)
        p[, i] <- matrixcalc::vec(xy)
      }
      else if (family == "poisson") {
        A <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
        diag(A) <- mu
        D <- mat.prod(covariate, mu)
        Vi <- diag(sqrt(c(mu)),
                   cluster$n[i]) %*% var_i %*% diag(sqrt(c(mu)), cluster$n[i])
        xy <- t(D) %*% solve(Vi) %*% sqrt(A) %*% unstr %*%
          sqrt(A) %*% solve(Vi) %*% D
        xx <- t(D) %*% solve(Vi) %*% D
        step12 <- step12 + xy
        step11 <- step11 + xx
        step13 <- step13 + matrixcalc::vec(xy)
        p[, i] <- matrixcalc::vec(xy)
      }
    }else if(small.samp.method=="md"){
      if (family == "gaussian") {
        xy <- t(covariate) %*% solve(var_i) %*% solve(cormax.ind(cluster$n[i]) -
                                                        covariate %*% solve(step11) %*% t(covariate) %*%
                                                        solve(var_i)) %*% (y - mu)
        step12 <- step12 + xy %*% t(xy)
        step13 <- step13 + matrixcalc::vec(xy %*% t(xy))
        p[, i] <- matrixcalc::vec(xy %*% t(xy))
      }
      else if (family == "poisson") {
        D <- mat.prod(covariate, mu)
        Vi <- diag(sqrt(c(mu)),
                   cluster$n[i]) %*% var_i %*% diag(sqrt(c(mu)), cluster$n[i])
        xy <- t(D) %*% solve(Vi) %*% solve(cormax.ind(cluster$n[i]) -
                                             D %*% solve(step11) %*% t(D) %*% solve(Vi)) %*%
          (y - mu)
        step12 <- step12 + xy %*% t(xy)
        step13 <- step13 + matrixcalc::vec(xy %*% t(xy))
        p[, i] <- matrixcalc::vec(xy %*% t(xy))
      }
    }else if(small.samp.method=="wl"){
      if (family == "gaussian") {
        xy <- t(covariate) %*% solve(var_i) %*% unstr %*%
          solve(var) %*% covariate
        xx <- t(covariate) %*% solve(var_i) %*% covariate
        step11 <- step11 + xx
        step12 <- step12 + xy
        step13 <- step13 + matrixcalc::vec(xy)
        p[, i] <- matrixcalc::vec(xy)
      }
      else if (family == "poisson") {
        B <- matrix(0, nrow = cluster$n[i], ncol = cluster$n[i])
        diag(B) <- mu
        D <- mat.prod(covariate, mu)
        Vi <- diag(sqrt(c(mu)),
                   cluster$n[i]) %*% var_i %*% diag(sqrt(c(mu)), cluster$n[i])
        xy <- t(D) %*% solve(Vi) %*% sqrt(B) %*% unstr %*%
          sqrt(B) %*% solve(Vi) %*% D
        xx <- t(D) %*% solve(Vi) %*% D
        step11 <- step11 + xx
        step12 <- step12 + xy
        step13 <- step13 + matrixcalc::vec(xy)
        p[, i] <- matrixcalc::vec(xy)
      }
    }
  }
  for (i in 1:size) {
    dif <- (p[, i] - step13/size) %*% t(p[, i] - step13/size)
    step14 <- step14 + dif
  }
  cov.beta <- solve(step11) %*% (step12) %*% solve(step11)
  cov.var <- size/(size - 1) * kronecker(solve(step11), solve(step11)) %*%
    step14 %*% kronecker(solve(step11), solve(step11))
  gee.fit$small.samp.var<-diag(cov.beta)
  names(gee.fit$small.samp.var)<-names(gee.fit$coefficients)
  }
  gee.fit$call<-call
  return(gee.fit)
}

