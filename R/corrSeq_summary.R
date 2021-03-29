#' Function to summarize individual regression coefficients
#'
#' Conducts t-tests on individual regression coefficients from models fit from the \code{\link{corrSeq_fit}} function.
#'
#' @param corrSeq_results Results object from running \code{\link{corrSeq_fit}}.
#' @param coefficient Character string or numeric indicator of which coefficient to summarize
#' @param p_adj_method Method for adjusting for multiple comparisons (default is Benjamini-Hochberg). See \code{\link[stats]{p.adjust.methods}}.
#' @param ddf Method for computing degrees of freedom and t-statistics.
#' The options "Satterthwaite" and "Kenward-Roger" can only be used for
#' models fit using nbmm_pl or lmm. Options "containment" and
#' "residual" can be used for models fit using any method. Alternatively, a single numeric value representing the ddf for all tests can also be given.
#' @param sort_results Should the results table be sorted by adjusted p-value?
#'
#' @return This function returns a list object with the following components:
#'
#' \item{coefficient}{Name of the coefficient being summarized.}
#' \item{summary_table}{A summary table including the gene name, estimate, standard error, degrees of freedom, test statistic, and raw and adjusted p-value.}
#' \item{ddf}{Method for computing the degrees of freedom.}
#' \item{p_adj_method}{Method for adjusting the raw p-values.}
#' \item{singular_fits}{Gene names for genes that resulted in singular model fits. The summary information for these genes will be NA. Not applicable for models fit using \code{"gee"}.}
#' \item{method}{Method used to fit the models.}
#'
#'
#'@author Elizabeth Wynn
#'
#' @seealso \code{\link{corrSeq_fit}} \code{\link{geeglm_small_samp}}, \code{\link{glmm_nb_lmer}}, \code{\link[lmerTest]{lmer}}, \code{\link[glmmADMB]{glmmadmb}}
#' @examples
#' data("simdata")
#' sample_meta_data <- simdata$metadata
#'
#' ## Subset down to 10 observation (i.e. gene)
#' counts=simdata$counts[1:10,]
#
#'
#' ## Fit NBMM-PL models
#' ## Use log(library size) as an offset
#' nbmm_pl_fit <- corrSeq_fit(formula = ~ group * time+(1|ids)+offset(log(lib_size)),
#'                            expr_mat = counts,
#'                            sample_data = sample_meta_data,
#'                            method="nbmm_pl")
#'
#'
#' ## Summarize the group coefficient with Satterthwaite degrees of freedom
#' model_sum <-corrSeq_summary(corrSeq_results = nbmm_pl_fit,
#'                              coefficient = "group",
#'                              p_adj_method = 'BH',
#'                              ddf = 'Satterthwaite',
#'                              sort_results = T)
#'
#'
#' @export
#'

corrSeq_summary <- function(corrSeq_results = NULL, # Results object from running lmerSeq.fit
                            coefficient = NULL, # Character string or numeric indicator of which coefficient to summarize
                            p_adj_method = "BH", # Method for adjusting for multiple comparisons (default is Benjamini-Hochberg)
                            ddf = "Satterthwaite", # Method for computing degrees of freedom and t-statistics. Options are "Satterthwaite" and "Kenward-Roger"
                            sort_results = T # Should the results table be sorted by adjusted p-value?
){
  if(identical(names(corrSeq_results[[1]]),c("fit", "gene"))){
    method="lmm"
    ddf_methods=c("Satterthwaite", "Kenward-Roger", "containment", "residual")
    if(!(ddf %in% ddf_methods)&!is.numeric(ddf)){
      stop("Invalid ddf method")
    }
    #If using a different method, first just analyze using satterthwait
    #This gives us t-statistic, etc.
    ddf2=ddf
    if(!(ddf %in% c("Satterthwaite", "Kenward-Roger"))) ddf2="Satterthwaite"
    ret2=lmerSeq.summary(corrSeq_results, coefficient = coefficient,
                         p_adj_method = p_adj_method,
                         ddf=ddf2, sort_results = sort_results)
    #First non-null model
    if(ddf!=ddf2){
      if(!is.numeric(ddf)){
        idx_non_null_1<-which(sapply(corrSeq_results, function(x) !is.null(x$fit)))[1]
        ddf=calc_ddf(corrSeq_results[[idx_non_null_1]], ddf = ddf, method=method)
      }
      ret2$summary_table<-ret2$summary_table%>%dplyr::mutate(df=ifelse(is.na(df), NA, ddf), p_val_raw=2*pt(-abs(t.value),df=df),
                                                      p_val_adj=p.adjust(p_val_raw,
                                                                         method=p_adj_method))
      ret2$ddf=match.call()$ddf
    }

  }else{
    #Check if models are Null
    if( sum(sapply(corrSeq_results, is.null))){
      stop("Model fits for all genes are null")
    }
    #Get first non-null model
    idx_non_null_1=which(!sapply(corrSeq_results, is.null))[1]
    #For gee models
    if("geeglm"%in%class(corrSeq_results[[idx_non_null_1]])){
      method="gee"
      #Get coef names
      coef_names <- names(corrSeq_results[[idx_non_null_1]]$coefficients)
      #DDF methods for gee
      ddf_methods=c("containment", "residual")
      idx_not_singular=1:length(corrSeq_results)
    }else if(class(corrSeq_results[[idx_non_null_1]])=="glmm_nb_mod"){
      method="nbmm_pl"
      coef_names<-names(lme4::fixef(corrSeq_results[[idx_non_null_1]]))
      ddf_methods=c("containment", "residual", "Satterthwaite", "Kenward-Roger")
      idx_singular<-which(sapply(corrSeq_results, lme4::isSingular))
      idx_not_singular <- which(!sapply(corrSeq_results, lme4::isSingular))
    }else if(class(corrSeq_results[[idx_non_null_1]])=="glmmadmb"){
      method="nbmm_ml"
      coef_names=names(coef(corrSeq_results[[idx_non_null_1]]))
      idx_singular<-which(sapply(corrSeq_results, function(x) x$S$ids[1]<1e-05))
      idx_not_singular<-which(!sapply(corrSeq_results, function(x) x$S$ids[1]<1e-05))
      ddf_methods=c("containment", "residual")
    }


    gene_names <- names(corrSeq_results)
    ############################################################################################################
    #Error Messages for insufficient or inconsistent information
    ############################################################################################################
    if(is.numeric(coefficient)){
      if((coefficient > length(coef_names)) | coefficient < 1){
        stop("Coefficient number is invalid")
      }
      coef_out <- coef_names[coefficient]
    }
    else{
      if(!(coefficient %in% coef_names)){
        stop("Coefficient name is invalid")
      }
      coef_out <- coefficient
    }

    if(!(p_adj_method %in% p.adjust.methods)){
      stop("Invalid p_adj_method")
    }

    if(!(ddf %in% ddf_methods)&!is.numeric(ddf)){
      stop("Invalid ddf method")
    }

    if(ddf %in% c("Satterthwaite", "Kenward-Roger")){
      ret <- do.call(rbind, lapply(corrSeq_results[idx_not_singular], function(x){
        # x = corrSeq_results$fitted_models[[1]]
        res_sub <- summary(x, ddf = ddf)$coefficients[coefficient, ]
        return(res_sub)
      }))%>%data.frame()%>%dplyr::rename(p_val_raw ="Pr...t..")%>%dplyr::mutate(Gene=gene_names[idx_not_singular],
                                                                  p_val_adj = p.adjust(p_val_raw, method = p_adj_method))%>%
        dplyr::select(Gene, Estimate, Std.Error="Std..Error", "df","t.value",
               p_val_raw, p_val_adj)


    }else{
      #If df will be same for all models, calculate now
      constant_ddf_methods=c("residual", "containment")
      if(ddf %in% constant_ddf_methods){
        ddf=calc_ddf(model = corrSeq_results[[1]],ddf = ddf, method=method)
      }
      ret=lapply(corrSeq_results[idx_not_singular], function(x){

        #Estimate and std. error for gee
        if(method=="gee"){
          Estimate=x$coefficients[coefficient]
          #if small sample method was used
          if(!is.null(x$small.samp.va)){
            Std.Error=sqrt(x$small.samp.var[coefficient])
          }else Std.Error=summary(x)$coefficients[coefficient,"Std.err"]
        }else if(method=="nbmm_pl"){
          Estimate=summary(x)$coefficients[coefficient,"Estimate"]
          Std.Error=summary(x)$coefficients[coefficient,"Std. Error"]
        }else if(method=="nbmm_ml"){
          Estimate=summary(x)$coefficient[coefficient, "Estimate"]
          Std.Error=summary(x)$coefficient[coefficient, "Std. Error"]
        }
        if(!is.numeric(ddf)){
          ddf=calc_ddf(model=x, ddf=ddf, method=method)
        }
        t.value=Estimate/Std.Error
        p_val_raw=2*pt(-abs(t.value),
                       df=ddf)

        data.frame(Estimate=Estimate, Std.Error=Std.Error, df=ddf,
                   t.value=t.value, p_val_raw=p_val_raw)
      })%>%dplyr::bind_rows()%>%
        dplyr::mutate(Gene=gene_names[idx_not_singular], p_val_adj=p.adjust(p_val_raw, method = p_adj_method))%>%
        dplyr::select(Gene, Estimate, Std.Error, df, t.value, p_val_raw, p_val_adj)

    }



    if(sort_results){
      ret <- ret %>%
        dplyr::arrange(p_val_adj)
    }
    rownames(ret) <- NULL
    ret2 <- list(coefficient = coef_out,
                 summary_table = ret,
                 ddf = match.call()$ddf,
                 p_adj_method = p_adj_method)
    if(method !="gee"){
      browser()
      genes_singular_fits <- as.character(gene_names[idx_singular])
      ret2$singular_fits = genes_singular_fits
      ret2$summary_table$Gene<-as.character(ret2$summary_table$Gene)
      ret2$summary_table <- gtools::smartbind(ret2$summary_table, data.frame(Gene = genes_singular_fits))

    }
  }
  ret2$model_method=method
  rownames(ret2$summary_table)<-NULL
  return(ret2)
}


########################################
#Calculate containment and residual ddf
########################################
calc_ddf<-function(model, ddf, method){
  if(method=="lmm"|method=="nbmm_pl"|method=="nbmm_ml"){
    #if nbmm_ml, fit lmm in order to get zmat and xmat
    if(method=="nbmm_ml") model=lme4::lmer(formula = model$formula, data=model$data)
    if(method=="lmm") model=model$fit
    ids<-names(ranef(model))
    if(length(unique(table(model@frame[,ids])))!=1) stop("Must have the same number of measurements for each subject")
    x_mat=matrix(getME(model, "X"), ncol =ncol(getME(model, "X")))
    if(ddf=="residual"){
      x_rank=qr(x_mat)$rank
      ddf_val=nrow(x_mat)-x_rank
    }else if(ddf=="containment"){
      z_mat<-matrix(getME(model, "Z"), ncol =ncol(getME(model, "Z")) )
      xz_mat<-cbind(z_mat, x_mat)
      xz_rank<-qr(xz_mat)$rank
      ddf_val=nrow(model@frame)-xz_rank
    }
  }else if(method=="gee"){
    if(length(unique(table(model$id)))!=1) stop("Must have the same number of measurements for each subject")
    x_mat=model.matrix(model$formula, model$data)
    if(ddf=="residual"){
      x_rank=qr(x_mat)$rank
      ddf_val=nrow(model$data)-x_rank
    }else if(ddf=="containment"){
      z_mat=model.matrix(~-1+as.factor(model$id))
      xz_mat<-cbind(z_mat, model$model)
      xz_rank<-qr(xz_mat)$rank
      ddf_val=nrow(model$data)-xz_rank
    }
  }
  return(ddf_val)
}






