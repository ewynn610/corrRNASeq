#' Calculate Sigma and G matrices for Pseudo-likelihood Negative Binomial Mixed Model
#'
#' Calculates the Sigma and G components described in Kenward and
#' Roger (1997) while accounting for \code{glmm_nb_mod} object
#' weights. Function is used in calculating Kenward-Roger degrees of
#' freedom in functions from the \pkg{pbkrtest}-package such as \code{KRmodcomp}.
#' @importFrom pbkrtest get_SigmaG
#' @param object A \code{glmm_nb_mod} model object from the \code{glmm_nb_lmer} function
#' @param details If larger than 0 some timing details are printed.
#' @return
#' \item{Sigma}{The covariance matrix of Y}
#' \item{G}{The G matrices that sum up to Sigma}
#' \item{n.ggamma}{The number of G matrices (Called M in Kenward and Roger 1997)}
#'
#' @author Elizabeth Wynn and Camille Moore, underlying code drawn from code by \pkg{pbkrtest}-authors.
#'
#' @seealso \code{\link{glmm_nb_mod}}, \code{\link[pbkrtest]{get_SigmaG}}, \code{\link[pbkrtest]{KRmodcomp}}
#'
#' @references
#' Kenward MG, Roger JH (1997). "Small Sample Inference for Fixed Effects from Restricted Maximum Likelihood." \emph{Biometrics}, 53(3), 983. ISSN 0006341X. doi: 10.2307/2533558.
#'
#'@export
#'@importFrom methods as signature
#'


#Code from get_SigmaG.lmerMod function except where indicated
get_SigmaG.glmm_nb_mod  <- function(object, details=0) {
  DB     <- details > 0 ## For debugging only

  GGamma <- VarCorr(object)
  SS     <- .shgetME( object )

  ## Put covariance parameters for the random effects into a vector:
  ## Fixme: It is a bit ugly to throw everything into one long vector here; a list would be more elegant
  ggamma <- NULL
  for ( ii in 1:( SS$n.RT )) {
    Lii    <- GGamma[[ii]]
    ggamma <- c(ggamma, Lii[ lower.tri( Lii, diag=TRUE ) ] )
  }
  ggamma   <- c( ggamma, sigma( object )^2 ) ## Extend ggamma by the residuals variance
  n.ggamma <- length(ggamma)

  ## Find G_r:
  G  <- NULL
  Zt <- getME( object, "Zt" )
  for (ss in 1:SS$n.RT) {
    ZZ    <- .shget_Zt_group( ss, Zt, SS$Gp )
    n.lev <- SS$n.lev.by.RT2[ ss ] ## ; cat(sprintf("n.lev=%i\n", n.lev))
    Ig    <- sparseMatrix(1:n.lev, 1:n.lev, x=1)
    for (rr in 1:SS$n.parm.by.RT[ ss ]) {
      ## This is takes care of the case where there is random regression and several matrices have to be constructed.
      ## FIXME: I am not sure this is correct if there is a random quadratic term. The '2' below looks suspicious.
      ii.jj <- .index2UpperTriEntry( rr, SS$n.comp.by.RT[ ss ] ) ##; cat("ii.jj:"); print(ii.jj)
      ii.jj <- unique(ii.jj)
      if (length(ii.jj)==1){
        EE <- sparseMatrix(ii.jj, ii.jj, x=1, dims=rep(SS$n.comp.by.RT[ ss ], 2))
      } else {
        EE <- sparseMatrix(ii.jj, ii.jj[2:1], dims=rep(SS$n.comp.by.RT[ ss ], 2))
      }
      EE <- Ig %x% EE  ## Kronecker product
      G  <- c( G, list( t(as.matrix(ZZ)) %*% EE %*% ZZ ) ) #changed to as.matrix because t() won't work on sparse matrix
    }
  }

  ## Extend by the indentity for the residual
  n.obs <- nrow(getME(object,'X'))

  ##########Adjusted from pbkrtest package to include weights##########
  #G    <- c( G, list(sparseMatrix(1:n.obs, 1:n.obs, x=1 )) ) #### CMM Old Code
  G    <- c( G, list(Matrix(diag(1/attributes(object)$frame$"(weights)"), sparse=T)) ) ### CMM new line to use weights

  Sigma <- ggamma[1] * G[[1]]
  for (ii in 2:n.ggamma) {
    Sigma <- Sigma + ggamma[ii] * G[[ii]]
  }

  SigmaG <- list(Sigma=Sigma, G=G, n.ggamma=n.ggamma)
  SigmaG
}

##############################################
#Helper functions from the pbkrtest package
##############################################

.shgetME<-function (object)
{
  Gp <- getME(object, "Gp")
  n.RT <- length(Gp) - 1
  n.lev.by.RT <- sapply(getME(object, "flist"), function(x) length(levels(x)))
  n.comp.by.RT <- .get.RT.dim.by.RT(object)
  n.parm.by.RT <- (n.comp.by.RT + 1) * n.comp.by.RT/2
  n.RE.by.RT <- diff(Gp)
  n.lev.by.RT2 <- n.RE.by.RT/n.comp.by.RT
  list(Gp = Gp, n.RT = n.RT, n.lev.by.RT = n.lev.by.RT, n.comp.by.RT = n.comp.by.RT,
       n.parm.by.RT = n.parm.by.RT, n.RE.by.RT = n.RE.by.RT,
       n.lev.by.RT2 = n.lev.by.RT2, n_rtrms = getME(object,
                                                    "n_rtrms"))
}




.get.RT.dim.by.RT<-function (object)
{
  qq <- if (inherits(object, "mer")) {
    sapply(object@ST, function(X) nrow(X))
  }
  else {
    sapply(object@cnms, length)
  }
  qq
}

.shget_Zt_group<-function (ii.group, Zt, Gp, ...)
{
  zIndex.sub <- (Gp[ii.group] + 1):Gp[ii.group + 1]
  ZZ <- Zt[zIndex.sub, ]
  return(ZZ)
}
.index2UpperTriEntry<-function (k, N)
{
  aa <- cumsum(N:1)
  aaLow <- c(0, aa[-length(aa)])
  i <- which(aaLow < k & k <= aa)
  j <- k - N * i + N - i * (3 - i)/2 + i
  return(c(i, j))
}
