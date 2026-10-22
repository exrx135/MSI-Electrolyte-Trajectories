# Run with Rscript --vanilla; requires mice, sandwich, jsonlite, ggplot2, ragg.
# New, explicitly authorized covariate-adjusted analysis. Frozen cohorts,
# phenotype assignments and MAKE30 outcomes are validated and never re-estimated.
if (.Platform$OS.type=='windows') Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8')
args0 <- commandArgs(trailingOnly=FALSE)
file_arg <- sub('^--file=', '', args0[grepl('^--file=', args0)])
P <- normalizePath(file.path(dirname(file_arg), '..'), winslash='/', mustWork=TRUE)
if (!length(file_arg)) stop('Run this file using Rscript --vanilla.')
options(stringsAsFactors=FALSE, warn=1)
library(mice)
library(sandwich)
library(jsonlite)
spec <- fromJSON(file.path(P,'model_specification.json'))
M2 <- spec$model2; M3 <- spec$model3
db_run <- commandArgs(trailingOnly=TRUE)
stopifnot(length(db_run)<=1,all(db_run %in% c('MIMIC','eICU','SICdb')))
csv <- function(x, filename) {
  combined <- c('pooled_RR_results.csv','model_fit_audit.csv','per_imputation_coefficients.csv','imputed_distributions.csv','crude_direct_ratio_check.csv')
  if(length(db_run) && basename(filename)%in%combined) filename <- file.path(dirname(filename),paste0(db_run,'_',basename(filename)))
  write.csv(x,file.path(P,filename),row.names=FALSE,na='')
}
all_results <- list(); fits_audit <- list(); imputation_summary <- list()
coef_draws <- list(); warning_rows <- list(); crude_checks <- list()
for (db in if(length(db_run)) db_run else c('MIMIC','eICU','SICdb')) {
  cat('\nSTART ',db,' ',format(Sys.time()),'\n',sep=''); flush.console()
  raw <- read.csv(file.path(P,'data',paste0(db,'_analysis_input.csv')))
  d <- raw[,c('make30',M3)]
  # Encode observed auxiliary information, not a filled-in analysis weight.
  # The centered weight term is zero when weight is missing; the accompanying
  # indicator distinguishes this pattern from an observed median weight.
  weight_center <- median(raw$weight_kg,na.rm=TRUE)
  d$weight_observed_centered <- ifelse(is.na(raw$weight_kg),0,raw$weight_kg-weight_center)
  d$weight_missing <- factor(as.integer(is.na(raw$weight_kg)),levels=c(0,1))
  d$phenotype <- factor(paste0('C',raw$assigned_class),levels=c('C1','C2','C3'))
  for (v in c('male','hypertension','diabetes','ckd_or_renal_disease')) d[[v]] <- factor(d[[v]],levels=c(0,1))
  stopifnot(!anyNA(d$make30), !anyNA(d$phenotype), max(colMeans(is.na(d)))<.30)
  meth <- mice::make.method(d); meth[] <- ''
  incomplete <- names(d)[colSums(is.na(d))>0]
  for(v in incomplete) meth[v] <- if(is.factor(d[[v]])) 'logreg' else 'pmm'
  pred <- mice::make.predictorMatrix(d); diag(pred)<-0
  pred[setdiff(names(d),incomplete),]<-0
  seed <- spec$seed + match(db,c('MIMIC','eICU','SICdb'))*1000
  midfile <- file.path(P,'diagnostics',paste0(db,'_mice_final.rds'))
  # Continue longer chains only if the saved inputs and imputation specification
  # match exactly. This permits convergence-directed extension without refitting
  # trajectories or repeatedly starting the same chain from scratch.
  if(file.exists(midfile)) {
    imp <- readRDS(midfile)
    stopifnot(identical(imp$data,d),identical(imp$method,meth),
              identical(imp$predictorMatrix,pred),imp$m==spec$m,imp$seed==seed)
    if(imp$iteration<spec$iterations) imp<-mice::mice.mids(imp,maxit=spec$iterations-imp$iteration,printFlag=FALSE)
  } else {
    imp <- mice(d,m=spec$m,maxit=spec$iterations,method=meth,predictorMatrix=pred,
                seed=seed,printFlag=FALSE,donors=5,ridge=1e-5)
  }
  saveRDS(imp,midfile)
  csv(data.frame(variable=names(meth),method=meth),paste0('diagnostics/',db,'_imputation_methods.csv'))
  if(!is.null(imp$loggedEvents)) csv(imp$loggedEvents,paste0('diagnostics/',db,'_mice_logged_events.csv'))
  # Full trace panels and dependency-free split-Rhat are produced by
  # check_imputations.R; optional rstan diagnostics are not required.
  cat('IMPUTATION COMPLETE ',db,' ',format(Sys.time()),'\n');flush.console()
  completed <- lapply(seq_len(spec$m),function(i)mice::complete(imp,i))
  for (j in seq_along(completed)) {
    stopifnot(!anyNA(completed[[j]]))
    for(v in names(d)) {
      observed <- !is.na(d[[v]])
      stopifnot(all(as.character(completed[[j]][[v]][observed])==as.character(d[[v]][observed])))
    }
  }
  for(v in incomplete) {
    obs <- as.numeric(as.character(d[[v]][!is.na(d[[v]])]))
    vals <- unlist(lapply(completed,function(x)as.numeric(as.character(x[[v]][is.na(d[[v]])]))))
    imputation_summary[[length(imputation_summary)+1L]] <- data.frame(dataset=db,variable=v,missing_n=sum(is.na(d[[v]])),observed_min=min(obs),observed_median=median(obs),observed_max=max(obs),imputed_min=min(vals),imputed_median=median(vals),imputed_max=max(vals))
    if(meth[v]=='pmm') stopifnot(all(vals %in% obs))
  }
  fit_one <- function(dat, vars, model, j) {
    # Scaling only improves numerical conditioning; effects and predictions
    # are identical to an unpenalized model using the original units.
    continuous <- vars[!vapply(dat[,vars,drop=FALSE],is.factor,logical(1))]
    for(v in continuous) {
      center <- mean(d[[v]],na.rm=TRUE); sc <- sd(d[[v]],na.rm=TRUE)
      dat[[v]] <- (dat[[v]]-center)/sc
    }
    form <- reformulate(c('phenotype',vars),response='make30')
    warn <- character()
    fit <- withCallingHandlers(glm(form,data=dat,family=poisson(link='log'),control=glm.control(maxit=100,epsilon=1e-9)),warning=function(w){warn <<- c(warn,conditionMessage(w));invokeRestart('muffleWarning')})
    if(!fit$converged || anyNA(coef(fit)) || any(!is.finite(coef(fit)))) stop(db,' ',model,' nonconvergent or rank deficient')
    vc <- sandwich::vcovHC(fit,type='HC1')
    stopifnot(all(is.finite(vc)),all(diag(vc)>0),nobs(fit)==nrow(raw))
    mm <- model.matrix(fit)
    # Continuous/factor column collinearity diagnostic, independent of robust SE.
    VIF <- if(ncol(mm)>2) diag(solve(cor(mm[,-1,drop=FALSE]))) else 1
    lev <- hatvalues(fit)
    fits_audit[[length(fits_audit)+1L]] <<- data.frame(dataset=db,model=model,imputation=j,n=nobs(fit),parameters=length(coef(fit)),converged=fit$converged,max_fitted_risk=max(fitted(fit)),fitted_risk_gt1=sum(fitted(fit)>1),max_hat=max(lev),max_vif=max(VIF),warning=paste(warn,collapse='; '))
    data.frame(dataset=db,model=model,imputation=j,comparison=c('C2 vs C1','C3 vs C1'),Q=coef(fit)[c('phenotypeC2','phenotypeC3')],U=diag(vc)[c('phenotypeC2','phenotypeC3')],n=nobs(fit),k=length(coef(fit)))
  }
  for(model in c('Model 1','Model 2','Model 3')) {
    vars <- switch(model,'Model 1'=character(),'Model 2'=M2,'Model 3'=M3)
    draws <- if(model=='Model 1') fit_one(completed[[1]],vars,model,0) else do.call(rbind,lapply(seq_len(spec$m),function(j)fit_one(completed[[j]],vars,model,j)))
    coef_draws[[length(coef_draws)+1L]] <- draws
    for (comparison in c('C2 vs C1','C3 vs C1')) {
      z <- draws[draws$comparison==comparison,]
      if(model=='Model 1') {
        q <- z$Q; total <- z$U; df <- Inf; b <- 0; fmi <- 0; ubar <- z$U
        cl <- strsplit(comparison,' ')[[1]][1]
        direct <- mean(raw$make30[raw$assigned_class==as.numeric(sub('C','',cl))])/mean(raw$make30[raw$assigned_class==1])
        stopifnot(abs(exp(q)-direct)<1e-7)
        crude_checks[[length(crude_checks)+1L]] <- data.frame(dataset=db,comparison=comparison,model_rr=exp(q),direct_risk_ratio=direct,difference=exp(q)-direct)
      } else {
        pooled <- mice::pool.scalar(z$Q,z$U,n=z$n[1],k=z$k[1])
        q<-pooled$qbar;total<-pooled$t;df<-pooled$df;b<-pooled$b;fmi<-pooled$fmi;ubar<-pooled$ubar
      }
      se <- sqrt(total); crit <- qt(.975,df)
      all_results[[length(all_results)+1L]] <- data.frame(dataset=db,model=model,comparison=comparison,n=nrow(raw),events=sum(raw$make30),RR=exp(q),CI95_low=exp(q-crit*se),CI95_high=exp(q+crit*se),p=2*pt(-abs(q/se),df),log_RR=q,SE_log_RR=se,df=df,within_variance=ubar,between_variance=b,fmi=fmi,m=ifelse(model=='Model 1',1,spec$m),MCSE_log_RR=ifelse(model=='Model 1',0,sqrt(b/spec$m)))
    }
  }
  csv(do.call(rbind,all_results),'pooled_RR_results.csv')
  csv(do.call(rbind,fits_audit),'diagnostics/model_fit_audit.csv')
  csv(do.call(rbind,coef_draws),'diagnostics/per_imputation_coefficients.csv')
  csv(do.call(rbind,imputation_summary),'diagnostics/imputed_distributions.csv')
  csv(do.call(rbind,crude_checks),'diagnostics/crude_direct_ratio_check.csv')
  cat('FINISHED ',db,' ',format(Sys.time()),'\n');flush.console()
}
capture.output(sessionInfo(),file=file.path(P,'diagnostics',paste0(if(length(db_run)) paste0(db_run,'_') else '', 'R_sessionInfo.txt')))
cat('\nALL ANALYSES COMPLETE\n')
