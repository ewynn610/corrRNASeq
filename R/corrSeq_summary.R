#' Function to summarize individual regression coefficients
#'
#' Conducts hypothesis testing on models fit using the \code{\link{corrSeq_fit}} function.
#'
#' @param corrSeq_results Results object from \code{\link{corrSeq_fit}}.
#' @param corrSeq_results_reduced Results object from \code{\link{corrSeq_fit}} containing reduced models for LRT tests.
#' Only applicable for nbmm_agq models or when using multi-row contrasts (matrix) for nbmm_lp models.
#' @param coefficient Character string or numeric indicator of which coefficient to summarize. Ignored if contrast is specified.
#' @param contrast Numeric vector or matrix specifying a contrast of the linear model coefficients to be tested.
#' Number of columns must equal the number of coefficients in the model. If specified, then takes precedence over coefficient.
#' @param p_adj_method Method for adjusting for multiple comparisons (default is Benjamini-Hochberg). See \code{\link[stats]{p.adjust.methods}}.
#' @param df Method for computing degrees of freedom for t- and F-tests.
#' The options "Satterthwaite" and "Kenward-Roger" can only be used for
#' models fit using nbmm_pl or lmm. Options "containment" and
#' "residual" can be used for models fit using any method except nbmm_agq (which does not use degrees of freedom, so use df=NA).
#' Alternatively, a single numeric value representing the df for all tests can also be given.
#' If testing a multi-row contrast for nbmm_lp, nbmm_agq, or gee, use df=NA since these tests do not use degrees of freedom.
#' If using a multi-row contrast for nbmm_pl or lmm, only df="Satterthwaite" and df="Kenward-Roger" are available.
#' @param sort_results Should the results table be sorted by adjusted p-value?
#' @param include_singular Should singular genes be included in results table?
#'
#' @details For single DF tests (single line contrasts or testing a single coefficient) all methods use a t-test except nbmm_agq.
#' For multiple DF tests (multi-row contrasts), nbmm_pl and lmm use an F-test, nbmm_lp uses a likelihood ratio test, and gee uses a Wald test.
#' For both single DF and multiple DF tests, nbmm_agq will perform a Wald test if no reduced model is provided. Otherwise a LRT is used.
#'
#' No DF is required for Wald of LRT tests, so df=NA should be used. Satterthwaite and Kenward-Rogers are only available for lmm and nbmm_pl. For multi-row contrasts
#' for nbmm_pl and lmm, only Satterthwaite and Kenward-Rogers can be used.
#'
#' @return This function returns a list object with the following components:
#'
#' \item{coefficient}{Name of the coefficient being summarized (if given).}
#' \item{contrast_mat}{Contrast matrix (if given)}
#' \item{summary_table}{A summary table including the gene name, estimate, standard error, degrees of freedom, test statistic, and raw and adjusted p-value.}
#' \item{df}{Method for computing the degrees of freedom.}
#' \item{p_adj_method}{Method for adjusting the raw p-values.}
#' \item{singular_fits}{Gene names for genes that resulted in singular model fits. The summary information for these genes will be NA unless include_singular is set to TRUE. Not applicable for models fit using \code{"gee"}.}
#' \item{method}{Method used to fit the models.}
#'
#'
#'@author Elizabeth Wynn
#'
#' @seealso \code{\link{corrSeq_fit}} \code{\link{geeglm_small_samp}}, \code{\link{glmm_nb_lmer}}, \code{\link[lmerTest]{lmer}}, \code{\link[glmmADMB]{glmmadmb}}, \code{\link[GLMMadaptive]{mixed_model}}
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
#' ## Note, group and time are factors
#' nbmm_pl_fit <- corrSeq_fit(formula = ~ group * time+(1|ids)+offset(log_offset),
#'                            expr_mat = counts,
#'                            sample_data = sample_meta_data,
#'                            method="nbmm_pl")
#'
#'
#'
#' ## Test for differential expression between groups at any timepoints
#' contrast_mat<-rbind(c(0, 1, 0, 0, 0, 0, 0, 0), #Difference in groups at time1
#'                     c(0, 1, 0, 0, 0, 1, 0, 0), #Difference in groups at time2
#'                     c(0, 1, 0, 0, 0, 0, 1, 0), #Difference in groups at time3
#'                     c(0, 1, 0, 0, 0, 0, 0, 1)  #Difference in groups at time4
#'                     )
#'
#'
#' group_sum <-corrSeq_summary(corrSeq_results = nbmm_pl_fit,
#'                              contrast = contrast_mat,
#'                              p_adj_method = 'BH',
#'                              df = 'Satterthwaite',
#'                              sort_results = T)
#'
#' ## Test for differential expression between groups at time 1
#' group_diff_time1 <-corrSeq_summary(corrSeq_results = nbmm_pl_fit,
#'                              coefficient = "group1",
#'                              p_adj_method = 'BH',
#'                              df = 'Satterthwaite',
#'                              sort_results = T)
#'
#'
#' @export
#'

corrSeq_summary <- function(corrSeq_results = NULL, # Results object from running lmerSeq.fit
                            corrSeq_results_reduced=NULL, #
                            coefficient = NULL, # Character string or numeric indicator of which coefficient to summarize
                            contrast=NULL, #Matrix with matrix to be tested
                            p_adj_method = "BH", # Method for adjusting for multiple comparisons (default is Benjamini-Hochberg)
                            df = "residual", # Method for computing degrees of freedom and t-statistics. Options are "Satterthwaite" and "Kenward-Roger"
                            sort_results = T, # Should the results table be sorted by adjusted p-value?
                            include_singular=F
){
  if(!is.null(contrast)){
    corrSeq_results_reduced=NULL
  }
  if(!is.null(contrast)|!is.null(corrSeq_results_reduced)){
    coefficient=NULL
  }

  #Save df name in new variable (df is used to save values later)
  df_name=df

  #Get genenames (names from list)
  gene_names <- names(corrSeq_results)

  # Is it contrast?
  contrast_tf=!is.null(contrast)

  # Are reduced models provided?
  reduced_tf=!is.null(corrSeq_results_reduced)

  ## Is it a multiple coefficient test?
  joint_flag=ifelse(contrast_tf, nrow(contrast)>1, F)
  joint_flag=ifelse(is.na(joint_flag), F, joint_flag)

  ## First non-null index
  idx_non_null_1=which(!sapply(corrSeq_results, is.null))[1]

  ## First non-null index for reduced mods
  if(reduced_tf){
    idx_non_null_1_reduced=which(!sapply(corrSeq_results_reduced, is.null))[1]
  }

  ##LMM
  if(identical(names(corrSeq_results[[idx_non_null_1]]),c("fit", "gene"))){
    method="lmm"
    if(reduced_tf){
      stop("Reduced models should only be provided for nbmm_agq or nbmm_lp")
    }
    if(joint_flag){
      df_methods=c("Satterthwaite", "Kenward-Roger")
    }else df_methods=c("Satterthwaite", "Kenward-Roger", "containment", "residual")
    if(!(df %in% df_methods)&!is.numeric(df)){
      stop("Invalid df method for given the model fitting method")
    }
    #If using a different method, first just analyze using satterthwait
    #This gives us t-statistic, etc.
    df2=df
    if(!(df %in% c("Satterthwaite", "Kenward-Roger"))) df2="Satterthwaite"

    if(contrast_tf){
      ret2=lmerSeq.contrast(corrSeq_results, contrast_mat = rbind(contrast),
                           p_adj_method = p_adj_method,
                           ddf=df2, sort_results = sort_results,
                           include_singular = include_singular)
      ## Remove upper, lower if one dimensional contrast
      # Just so results match
      if(!joint_flag) ret2$summary_table=ret2$summary_table%>%dplyr::select(-upper, -lower)%>%dplyr::rename(Std.Error="Std..Error")
    }else{
      ret2=lmerSeq.summary(corrSeq_results, coefficient = coefficient,
                           p_adj_method = p_adj_method,
                           ddf=df2, sort_results = sort_results)
      ret2$summary_table=ret2$summary_table%>%dplyr::rename(Std.Error="Std..Error")
    }

    #First non-null model
    if(df!=df2){
      if(!is.numeric(df)){
        idx_non_null_1<-which(sapply(corrSeq_results, function(x) !is.null(x$fit)))[1]
        df_new=calc_df(corrSeq_results[[idx_non_null_1]], df = df, method=method)
      }
      ret2$summary_table<-ret2$summary_table%>%dplyr::mutate(df=ifelse(is.na(df), NA, df_new), p_val_raw=2*pt(-abs(t.value),df=df_new),
                                                             p_val_adj=p.adjust(p_val_raw,
                                                                                method=p_adj_method))
    }
    #Capitolize column name so it matches with others
    colnames(ret2$summary_table)[colnames(ret2$summary_table)=="gene"]<-"Gene"
    #get rid of ddf (artifact from lmerseq, gets added later)
    ret2$ddf<-NULL
  }else{
    #Check if models are Null
    if( sum(sapply(corrSeq_results, function(x) !(is.null(x))))==0){
      stop("Model fits for all genes are null")
    }
    #Get first non-null model
    idx_non_null_1=which(!sapply(corrSeq_results, is.null))[1]
    #For gee models
    if("geeglm"%in%class(corrSeq_results[[idx_non_null_1]])){
      method="gee"
      #Get coef names
      coef_names <- names(corrSeq_results[[idx_non_null_1]]$coefficients)
      #df methods for gee
      if(joint_flag){
        df_methods=NA
      }else df_methods=c("containment", "residual")
      idx_not_converged<-which(sapply(corrSeq_results, is.null))
      idx_converged_not_singular=which(!(1:length(corrSeq_results)%in% idx_not_converged))
    }else if(class(corrSeq_results[[idx_non_null_1]])=="glmm_nb_mod"){
      method="nbmm_pl"
      coef_names<-names(lme4::fixef(corrSeq_results[[idx_non_null_1]]))
      if(joint_flag){
        df_methods=c("Satterthwaite", "Kenward-Roger")
      }else df_methods=c("containment", "residual", "Satterthwaite", "Kenward-Roger")
      idx_singular<-which(sapply(corrSeq_results, function(x) if(is.null(x)) F else lme4::isSingular(x)))
      idx_not_converged<-which(sapply(corrSeq_results, function(x) if(!is.null(x)) x@converged==F else T))
      if(include_singular){
        idx_converged_not_singular <- which(!(1:length(corrSeq_results)%in% c(idx_not_converged)))

      }else{
        idx_converged_not_singular <- which(!(1:length(corrSeq_results)%in% c(idx_not_converged, idx_singular)))

      }
    }else if(class(corrSeq_results[[idx_non_null_1]])=="glmmadmb"){
      method="nbmm_lp"
      coef_names=names(coef(corrSeq_results[[idx_non_null_1]]))
      idx_not_converged<-which(sapply(corrSeq_results, is.null))
      idx_singular<-which(sapply(corrSeq_results, function(x){
        any(sapply(x$S, function(x) any(diag(x)<1e-05)))
      })
      )
      if(include_singular){
        idx_converged_not_singular <- which(!(1:length(corrSeq_results)%in% c(idx_not_converged)))

      }else{
        idx_converged_not_singular <- which(!(1:length(corrSeq_results)%in% c(idx_not_converged, idx_singular)))

      }
      if(reduced_tf){
        if(class(corrSeq_results_reduced[[idx_non_null_1]])!="glmmadmb"){
          stop("Method for full and reduced models do not match")
        }
        df_methods=NA
      }else df_methods=c("containment", "residual")
    }else if(class(corrSeq_results[[idx_non_null_1]])=="MixMod"){
        method="nbmm_agq"
        if(reduced_tf){
        if(class(corrSeq_results_reduced[[idx_non_null_1_reduced]])!="MixMod"){
          stop("Method for full and reduced models do not match")
        }}
        coef_names=names(corrSeq_results[[idx_non_null_1]]$coefficients)
        idx_not_converged<-which(sapply(corrSeq_results, function(x) ifelse(is.null(x$converged), T, !x$converged)))
        idx_singular<-which(sapply(corrSeq_results, function(x){
          any(diag(x$D)<1e-05)
        }))
        if(include_singular){
          idx_converged_not_singular <- which(!(1:length(corrSeq_results)%in% c(idx_not_converged)))

        }else{
          idx_converged_not_singular <- which(!(1:length(corrSeq_results)%in% c(idx_not_converged, idx_singular)))

        }
        df_methods=NA
        }


    ############################################################################################################
    #Error Messages for insufficient or inconsistent information
    ############################################################################################################
    if(!contrast_tf& !reduced_tf){
      if(is.numeric(coefficient)){
        if((coefficient > length(coef_names)) | coefficient < 1){
          stop("Coefficient number is invalid")
        }
        coef_out <- coef_names[coefficient]
      }else{
        if(!(coefficient %in% coef_names)){
          stop("Coefficient name is invalid")
        }
        coef_out <- coefficient
      }
    }

    if(!(method%in% c("nbmm_agq", "nbmm_lp"))& reduced_tf){
      stop("Reduced models should only be provided for nbmm_agq or nbmm_lp")
    }

    if(!(p_adj_method %in% p.adjust.methods)){
      stop("Invalid p_adj_method")
    }

    if(method=="nbmm_lp"&contrast_tf){
      stop("Contrast testing not available for method nbmm_lp")
    }

    if(!(df %in% df_methods)&!is.numeric(df)){
      stop("Invalid df method")
    }

    if(df %in% c("Satterthwaite", "Kenward-Roger")){
      ret <- do.call(rbind, lapply(corrSeq_results[idx_converged_not_singular], function(x){
        # x = corrSeq_results$fitted_models[[1]]
        if(!contrast_tf){
          res_sub <- summary(x, ddf = df)$coefficients[coefficient, ]
          names(res_sub)[names(res_sub)=="Pr(>|t|)"]<-"p_val_raw"
          }else{
            res_sub<-tryCatch({
              res_sub=lmerTest::contest(x, contrast, ddf=df, joint=joint_flag)
              names(res_sub)[names(res_sub)%in%c("Pr(>F)","Pr(>|t|)")]<-"p_val_raw"
              res_sub
            },
            error=function(e){
              c("Sum Sq"=NA, "Mean Sq"=NA, "NumDF"=NA, "DenDF"=NA,
                         "F value"=NA, p_val_raw=NA)
            })
          }
        return(res_sub)
      }))%>%data.frame()%>%
        dplyr::mutate(Gene=gene_names[idx_converged_not_singular],
                      p_val_adj = p.adjust(p_val_raw, method = p_adj_method))
      if(joint_flag){
        ret=ret%>%dplyr::select(Gene, Sum.Sq, Mean.Sq,NumDF, DenDF,F.value,
                                p_val_raw, p_val_adj)
      }else{
        ret=ret%>%dplyr::select(Gene, Estimate, Std.Error="Std..Error", "df","t.value",
                                p_val_raw, p_val_adj)
      }

    }else{
      #If df will be same for all models, calculate now
      constant_df_methods=c("residual", "containment")
      if(df %in% constant_df_methods){
        df=calc_df(model = corrSeq_results[[idx_non_null_1]],df = df, method=method)
      }
      ret=lapply(idx_converged_not_singular, function(x){
        #### If using degrees of freedom (GEE, NBMM-LP, NBMM-PL singe effects tests)
        if(method=="gee"&!joint_flag){
          ## For single line contrasts
          if(contrast_tf){
            test=doBy::esticon(corrSeq_results[[x]], contrast, joint.test = F)
            Estimate=test$estimate
            Std.Error=test$std.error
          }else{ ## If they gave a coefficient
            Estimate=corrSeq_results[[x]]$coefficients[coefficient]
            #if small sample method was used
            if(!is.null(corrSeq_results[[x]]$small.samp.va)){
              Std.Error=sqrt(corrSeq_results[[x]]$small.samp.var[coefficient])
            }else Std.Error=summary(corrSeq_results[[x]])$coefficients[coefficient,"Std.err"]
          }
        }else if(method=="nbmm_pl"){
          if(contrast_tf){
            cont=lmerTest::contest(corrSeq_results[[x]], L=contrast, joint=F)$Estimate
            Estimate=cont$Estimate
            Std.Error=cont$`Std. Error`
          }else{
            Estimate=summary(corrSeq_results[[x]])$coefficients[coefficient,"Estimate"]
            Std.Error=summary(corrSeq_results[[x]])$coefficients[coefficient,"Std. Error"]
          }
        }else if(method=="nbmm_lp"& !reduced_tf){
          Estimate=summary(corrSeq_results[[x]])$coefficient[coefficient, "Estimate"]
          Std.Error=summary(corrSeq_results[[x]])$coefficient[coefficient, "Std. Error"]
        }
        if(!is.numeric(df)&!is.na(df)&method!="nbmm_agq"){
          res_df=calc_df(model=corrSeq_results[[x]], df=df, method=method)
        }

        ## NA degrees of freedom
        if(method!="nbmm_agq"&!reduced_tf&!(method=="gee"&joint_flag)){
          t.value=Estimate/Std.Error
          p_val_raw=2*pt(-abs(t.value),
                         df=df)
          res_df=data.frame(Estimate=Estimate, Std.Error=Std.Error, df=df,
                        t.value=t.value, p_val_raw=p_val_raw)
        }else if(method=="nbmm_agq"){
          if(contrast_tf){
            res_df=tryCatch({
              res_df=GLMMadaptive::anova(corrSeq_results[[x]], L=rbind(contrast))$aovTab.L
              colnames(res_df)[colnames(res_df)=="Pr(>|Chi|)"]<-"p_val_raw"
              res_df
            },error=function(e){
              data.frame(Chisq=NA, df=NA, p_val_raw=NA)
            })

          }else if(reduced_tf){
            res_df=tryCatch({
              my_aov=GLMMadaptive::anova(corrSeq_results[[x]], corrSeq_results_reduced[[x]])
              res_df=data.frame(LRT=my_aov$LRT, df=my_aov$df, p_val_raw=my_aov$p.value)
              res_df
            },error=function(e){
              data.frame(LRT=NA, df=NA, p_val_raw=NA)
            })
          }else{
            res_df=summary(corrSeq_results[[x]])$coef_table[coefficient,]
            names(res_df)[names(res_df)=="p-value"]<-"p_val_raw"
          }
        }else if(method=="nbmm_lp"&reduced_tf){
          res_df=tryCatch({
            my_aov=anova(corrSeq_results_reduced[[x]],corrSeq_results[[x]])
            res_df=data.frame(df=my_aov$Df[2], p_val_raw=my_aov$`Pr(>Chi)`[2])
            res_df
          },error=function(e){
            data.frame(df=NA, p_val_raw=NA)
          })
        }else if (method=="gee"& joint_flag){
            res_df=tryCatch({
              test=doBy::esticon(corrSeq_results[[x]], contrast, joint.test = T)
              data.frame(Chisq=test$X2.stat,df=test$DF, p_val_raw=test$`Pr(>|X^2|)`)
            },
            error=function(e){
              data.frame(Chisq=NA,df=NA, p_val_raw=NA)
            })

        }
        res_df
      })%>%dplyr::bind_rows()%>%
        dplyr::mutate(Gene=gene_names[idx_converged_not_singular], p_val_adj=p.adjust(p_val_raw, method = p_adj_method))
        if(method!="nbmm_agq"&!reduced_tf&(method=="gee"&!joint_flag)){
          ret=ret%>%dplyr::select(Gene, Estimate, Std.Error, df, t.value,
                                  p_val_raw, p_val_adj)
        }else if(method=="nbmm_agq"|(method=="gee"&joint_flag)){
          if(contrast_tf){
            ret=ret%>%dplyr::select(Gene, Chisq, df, p_val_raw, p_val_adj)
          }else if(reduced_tf){
            ret=ret%>%dplyr::select(Gene, LRT, df, p_val_raw, p_val_adj )
          }else ret=ret%>%dplyr::select(Gene, Estimate, Std.Err, "z-value", p_val_raw, p_val_adj)
        }else if(method=="nbmm_lp"& reduced_tf){
          ret=ret%>%dplyr::select(Gene,df, p_val_raw, p_val_adj )
        }

    }






    rownames(ret) <- NULL
    if(contrast_tf){
      ret2<-list(contrast_mat=contrast,summary_table = data.frame(ret),
                 p_adj_method = p_adj_method)
    }else if(reduced_tf){
      ret2<-list(summary_table = data.frame(ret),
                 p_adj_method = p_adj_method)
    } else{
      ret2 <- list(coefficient = coef_out,
                   summary_table = data.frame(ret),
                   p_adj_method = p_adj_method)
    }
    if(method !="gee"){
      genes_singular_fits <- as.character(gene_names[idx_singular])
      if(include_singular){
        genes_null=c(as.character(gene_names[idx_not_converged]))%>%unique()
      }else{
        genes_null=c(genes_singular_fits, as.character(gene_names[idx_not_converged]))%>%unique()
      }
      ret2$singular_fits = genes_singular_fits
    }else genes_null=as.character(gene_names[idx_not_converged])
    ret2$summary_table$Gene<-as.character(ret2$summary_table$Gene)
    ret2$summary_table <- gtools::smartbind(ret2$summary_table, data.frame(Gene = genes_null))
  }
  ret2$summary_table<-ret2$summary_table[match(gene_names,
                                               ret2$summary_table$Gene),]
  if(sort_results){
    ret2$summary_table <- ret2$summary_table %>%
      dplyr::arrange(p_val_adj)
  }
  ret2$model_method=method
  if(method!="nbmm_agq") ret2$df=df_name
  rownames(ret2$summary_table)<-NULL
  return(ret2)
}


########################################
#Calculate containment and residual df
########################################
calc_df<-function(model, df, method){
  if(method=="lmm"|method=="nbmm_pl"|method=="nbmm_lp"){
    #if nbmm_lp, fit lmm in order to get zmat and xmat
    if(method=="nbmm_lp") model=lme4::lmer(formula = model$formula, data=model$data)
    if(method=="lmm") model=model$fit
    ids<-names(ranef(model))
    #Mistake-not arequirement to have same number of measurements for each subject
    #if(length(unique(table(model@frame[,ids])))!=1) stop("Must have the same number of measurements for each subject")
    x_mat=matrix(getME(model, "X"), ncol =ncol(getME(model, "X")))
    if(df=="residual"){
      x_rank=qr(x_mat)$rank
      df_val=nrow(x_mat)-x_rank
    }else if(df=="containment"){
      z_mat<-matrix(getME(model, "Z"), ncol =ncol(getME(model, "Z")) )
      xz_mat<-cbind(z_mat, x_mat)
      xz_rank<-qr(xz_mat)$rank
      df_val=nrow(model@frame)-xz_rank
    }
  }else if(method=="gee"){
    #if(length(unique(table(model$id)))!=1) stop("Must have the same number of measurements for each subject")
    x_mat=model.matrix(model$formula, model$data)
    if(df=="residual"){
      x_rank=qr(x_mat)$rank
      df_val=nrow(model$data)-x_rank
    }else if(df=="containment"){
      z_mat=model.matrix(~-1+as.factor(model$id))
      xz_mat<-cbind(z_mat, model$model)
      xz_rank<-qr(xz_mat)$rank
      df_val=nrow(model$data)-xz_rank
    }
  }
  return(df_val)
}






