###############################################

#' Class for Pseudo-Likelihood Negative Binomial Mixed Models
#'
#' @importClassesFrom lme4 lmerMod
#' @importClassesFrom lme4 merMod
#' @importClassesFrom lmerTest lmerModLmerTest
#' @import methods
#'
#'
#' @slot dispersion Final negative binomial model dispersion estimate
#' @slot niter Number of iterations until model convergence (or iteration limit was reached)
#'
#' @return An object of class \code{glmm_nb_mod} with slots as in
#' \code{\link[lmerTest]{lmerModLmerTest}} and a few
#' additional slots as described below.
#'
#' @author Elizabeth Wynn
#'
#'
#' @export
#'
#'
#'

glmm_nb_mod<-setClass("glmm_nb_mod", contains = c("lmerModLmerTest"),
                      slots=c(dispersion="numeric", iter="numeric", converged="logical"))




#############################################################
#' Pseudo-Likelihood Negative Binomial Mixed Model
#'
#' This function fits a negative binomial mixed model using a
#' psuedo-likelihood approach
#'
#'
#' @return An object of class \code{\link{glmm_nb_mod}}
#'
#' @seealso \code{\link{glmm_nb_mod}}, \code{\link[lmerTest]{lmer}}
#'
#' @param formula A two-sided linear formula describing both the fixed-effects and random-effects parts of the model using the syntax of the lme4 package.
#' @param data An optional data frame containing the model variables.
#' @param niter Maximum number of iterations.
#' @param epsilon	Positive convergence tolerance.
#' @param verbose Logical. Should the number of iterations and computational time be printed?
#' @param REML Logical. Should the models be fit with REML or regular ML?
#'
#' @details #' This function is similar to the function \code{\link[NBZIMM]{glmm.nb}}
#' from the \pkg{NBZIMM} package, though it utilizes the \pkg{lmerTest} package
#' rather than \pkg{nlme} package and returns and object compatible with the \pkg{pkrtest} package
#' which allows for the calculation of Kenward-Roger degrees of freedom.
#' (see \code{\link[pbkrtest]{krmodcomp}},
#' \code{\link[lmerTest]{contest.lmerModLmerTest}})
#'
#' @author Elizabeth Wynn and Camille Moore, underlying code drawn from code by \pkg{NBZIMM}-authors.
#'
#' @references
#'
#' Zhang X, Mallick H, Tang Z, Zhang L, Cui X, Benson AK, Yi N (2017). "Negative binomial mixed models for analyzing microbiome count data." \emph{BMC Bioinformatics}, 18(1), 4. ISSN 14712105. doi:10.1186/s12859-016-1441-7
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
#' #Fit the Model
#' fit.glmm.nb <- glmm_nb_lmer(formula =counts ~ group * time + (1|ids),
#'                             data=df, REML = T)
#' @export
#'


glmm_nb_lmer<-function (formula, data, niter = 40, epsilon = 1e-08, verbose = FALSE, REML=TRUE){
  niter_theta=5
  niter_theta.ml = 5
  start.time <- Sys.time()

  #Function inside of packag
  family <- NegBin()

  m <- mcall <- Call <- match.call()


  nm <- names(m)[-1L]
  keep <- is.element(nm, c("weights", "data", "subset", "na.action")) #?
  for (i in nm[!keep]) m[[i]] <- NULL

  # Save all variable names, random and fixed
  allvars<-all.vars(formula)

  # Save formula of just fixed effects
  fixed<-lme4::nobars(formula)
  Terms <- if (missing(data)) terms(fixed) else terms(fixed, data = data)

  # Which term is the offset?
  off <- attr(Terms, "offset")
  off_name<-all.vars(fixed)[off]

  # If there is an offset, then it's saving on offset variable name
  if (length(off <- attr(Terms, "offset")))
    allvars <- c(allvars, as.character(attr(Terms, "variables"))[off + 1])


  #Save formula for output later
  Call$formula<-formula
  m$formula <- as.formula(paste(allvars[1],"~", paste(allvars[-c(1)], collapse = "+")))

  #Make sure m$formula environment the fixed environment (don't know why that's necessary)
  environment(m$formula) <- environment(fixed)

  #Set up call for model.frame function
  m$drop.unused.levels <- TRUE
  m[[1L]] <- quote(stats::model.frame)

  #Run model.frame -Returns df with counts, time, group, log_offset, id for each participant
  mf <- eval.parent(m)

  # Offset for the model
  off <- model.offset(mf)

  #sets offset to 0 if null
  if (is.null(off)) off <- 0

  # extracts weights if there are any
  wts <- model.weights(mf)

  # If weights is null, make all weights equal to 1 and save to mf
  if (is.null(wts)) wts <- rep(1, nrow(mf))
  mf$wts <- wts

  #Fit initial model
  fit0<-suppressWarnings(MASS::glm.nb(formula = fixed, data=mf))

  # Saving prior weights and linear predictors from model
  w <- fit0$prior.weights
  eta <- fit0$linear.predictors

  # Starting values for pseudo data
  zz <- eta + fit0$residuals

  # Starting values for weights
  wz <- fit0$weights

  #Don't know what this is for (again)
  nm <- names(mcall)[-1L]
  keep <- is.element(nm, c("fixed", "random", "data", "subset", "na.action", "control"))
  for (i in nm[!keep]) mcall[[i]] <- NULL

  # Replace the outcome is replaced with the transformation for LMM models on pseudo data
  formula[[2L]] <- quote(zz)
  mcall[["formula"]]<-formula #

  # Saves lmer function call (use lmerTest so we can get p-values)
  mcall[[1L]] <- quote(lmerTest::lmer)

  # Set REML to True or false based on function argument
  mcall$REML <- REML

  # Input updated weights
  mcall$weights<-wz

  # Add pseudo data to model data frame
  mf$zz<-zz

  #saves mf df (variables, weights, zz) as data
  mcall$data <- mf

  # Get outcome data
  y <- fit0$y

  # Vector to hold theta for each iteration
  my_theta <- vector()

  # Initialize number of iterations on calc_disp function
  j <- 0

  # Use theta from glm.nb in as first theta estimate
  fam <- NegBin(theta = fit0$theta)
  th <- th_old <- fam$theta
  vc<-NULL

  # for loop goes set number of iterations (specified in the call)
  for (i in seq_len(niter)) {
    #Evaluate function created from above info
    fit <- suppressMessages(eval(mcall))

    Z <- t(as.matrix(attributes(fit)$pp$Zt))
    X <- model.matrix(fit)

    # preveious model linear predictors
    etaold <- eta

    # new model linear predictors
    eta <- fitted(fit)

    # new model linear predictors (fixed effects only)
    eta_fixed <- X%*%attributes(fit)$beta + off


    # CMM: These lines calculate the expected value for each observation based on the model and the variance
    # uses formulas from the negative binomial distribution

    mu <- fam$linkinv(eta) #same as exp(eta)
    mu.eta.val <- fam$mu.eta(eta) #Always equal to mu: mu.eta is just doing pmax(exp(eta), .Machine$double.eps)

    #set mu.eta.val to small number if equal to 0
    mu.eta.val <- ifelse(mu.eta.val == 0, 1e-04, mu.eta.val)

    # Update the dispersion parameter
    th_old <- th
    # Calculate parameters to find dispersion
    my_m <- nrow(model.matrix(fit))-ncol(model.matrix(fit))

    # Get variance components for random effects
    vc_old<-vc
    if (j==0) vc_old_val<-0
    else vc_old_val<-vc_old[,4]
    vc <- lme4::VarCorr(fit)
    vc <- as.data.frame(vc,order="lower.tri")
    vc <- vc[vc$grp != "Residual",]
    vc_val<-vc[,4]

    # Create variance covariance matrix for all random effects (block diagonal matrix)
    G.mat <- matrix(NA,length(vc$var2[vc$var2=="<NA>"]),length(vc$var2[vc$var2=="<NA>"]))
    G.mat[lower.tri(G.mat,diag=TRUE)] <- vc$vcov
    G.mat <- makeSymm(G.mat)
    mat_list <- rep(G.mat, lme4::ngrps(fit))
    my_G <- Matrix::bdiag(  ## make a block-diagonal matrix
      lapply(
        split(mat_list, rep(1:lme4::ngrps(fit), each=length(G.mat))),
        matrix,dim(G.mat)))



    # Calculate dispersion
    error<-0

    tryCatch(
      th <- calc_dispersion(my_m, eta_fixed, my_G, mf$zz, mu, Z),
      error=function(err){error<<-1})
    if(error==0) {
      th<-th
      # Add 1 to number of iterations on new equation
      j <- j+1
    }
    if(error==1){ th<- suppressWarnings(
      MASS::theta.ml(y=y, mu=mu, n=sum(wts), weights=wts, limit=niter_theta.ml,
               trace=FALSE))
    # If need to use theta.ml, then j=0.
    j=0}

    fam <- NegBin(theta = th)

    #Output all thetas
    my_theta[i]<-th

    #stop loop if old and new fitted values not that different and old and new theta's not that different
    #And has run in the new equation 5 or more times
    #if (sum((eta - etaold)^2) < epsilon * sum(eta^2) & abs((th_old-th)^2)<epsilon*th^2 & ((vc_old_val-vc_val)^2)<epsilon*vc_val^2 & j>=niter_theta) break
    if (sum((eta - etaold)^2) < epsilon * sum(eta^2) & abs((th_old-th)^2)<epsilon*th^2) break


    # Update pseudo-response variable
    mf$zz <- eta + (y - mu)/mu.eta.val

    # Calculate variance
    varmu <- fam$variance(mu)

    # If varmu is 0 set to small number
    varmu <- ifelse(varmu == 0, 1e-04, varmu)


    # Update weights for LMM on pseudo responses
    wz <- w * mu.eta.val^2/varmu
    wz <- ifelse(wz == 0, 1e-04, wz)
    mcall$weights<-wz

    # update data (since zz changed)
    mcall$data <- mf
  }
  converged=T
  if (!(sum((eta - etaold)^2) < epsilon * sum(eta^2) & abs((th_old-th)^2)<epsilon*th^2)){
    warning("Model did not converge")
    converged=F
  }

  fit<-as(fit, "glmm_nb_mod")
  #fit$call outputs model formula
  attributes(fit)$call <- Call #CMM


  #save number of iterations
  attributes(fit)$iter <- i #CMM
  #logLik as NA
  ###fit$logLik <- as.numeric(NA) #I'd like to keep the LL or we could include it in a new slot called pseudo-likelihood

  #save theta value
  attributes(fit)$dispersion <- fam$theta #CMM

  #Save convergence information
  attributes(fit)$converged<-converged

  stop.time <- Sys.time()
  minutes <- round(difftime(stop.time, start.time, units = "min"), 3)

  if (verbose) {
    cat("Computational iterations:", attributes(fit)$iter, "\n") #CMM
    cat("Computational time:", minutes, "minutes \n")
  }

  #returns fitted model
  fit

}

###########################################################################################
# Helper Functions taken from NBZIMM package
###########################################################################################

# Negative Binomial Family Function
# theta = division dispersion
NegBin <- function (theta = 5, link = "log")
{
  linktemp <- substitute(link)
  if (!is.character(linktemp))
    linktemp <- deparse(linktemp)
  if (linktemp %in% c("log", "identity", "sqrt"))
    stats <- make.link(linktemp)
  else if (is.character(link)) {
    stats <- make.link(link)
    linktemp <- link
  }
  else {
    if (inherits(link, "link-glm")) {
      stats <- link
      if (!is.null(stats$name))
        linktemp <- stats$name
    }
    else stop(gettextf("\"%s\" link not available for negative binomial family; available links are \"identity\", \"log\" and \"sqrt\"",
                       linktemp))
  }
  .Theta <- theta
  env <- new.env(parent = .GlobalEnv)
  assign(".Theta", theta, envir = env)
  variance <- function(mu) mu + mu^2/.Theta
  validmu <- function(mu) all(mu > 0)
  dev.resids <- function(y, mu, wt) 2 * wt * (y * log(pmax(1, y)/mu) - (y + .Theta) * log((y + .Theta)/(mu + .Theta)))
  aic <- function(y, n, mu, wt, dev) {
    term <- (y + .Theta) * log(mu + .Theta) - y * log(mu) +
      lgamma(y + 1) - .Theta * log(.Theta) + lgamma(.Theta) -
      lgamma(.Theta + y)
    2 * sum(term * wt)
  }
  initialize <- expression({
    if (any(y < 0)) stop("negative values not allowed for the negative binomial family")
    n <- rep(1, nobs)
    mustart <- y + (y == 0)/6
  })
  simfun <- function(object, nsim) {
    ftd <- fitted(object)
    rnegbin(nsim * length(ftd), ftd, .Theta)
  }
  environment(variance) <- environment(validmu) <- environment(dev.resids) <- environment(aic) <- environment(simfun) <- env
  famname <- paste("NegBin(", format(round(theta, 4)), ")", sep = "")
  structure(list(family = famname, link = linktemp, linkfun = stats$linkfun,
                 linkinv = stats$linkinv, variance = variance, dev.resids = dev.resids,
                 aic = aic, mu.eta = stats$mu.eta, initialize = initialize,
                 validmu = validmu, valideta = stats$valideta, simulate = simfun, theta = .Theta),
            class = "family")
}

# Dispersion Calculation Function
# m = number observations - number parameters
# eta_fixed = linear predictor for each observation, including fixed effects and offsets
# G = covariance matrix of the random effects
# p = pseudo data
# Z random effects design matrix
# uses uniroot to solve for theta.  k = multiplicative dispersion. theta = 1/k

calc_dispersion <- function(m, eta_fixed, G, p, mu, Z){
  eq<-function(k){
    r <- p - eta_fixed
    V <- diag((mu+k*mu^2)/(mu^2)) + (Z %*%G %*%t(Z))
    return(
      as.numeric(m- t(r)%*%solve(V)%*%r)
    )
  }

  rootobj<-uniroot(eq, interval = c(0, 100))
  1/rootobj$`root`}

# Function to create a symmetrix matrix from a lower diagonal matrix
# Used in getting G matrix for variance components
# m = a lower diagonal matrix

makeSymm <- function(m) {
  m[upper.tri(m)] <- t(m)[upper.tri(m)]
  return(m)
}



