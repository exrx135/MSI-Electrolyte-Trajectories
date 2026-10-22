options(stringsAsFactors=FALSE, warn=1)
OUT <- "Z:/Research/MSI-1/analysis_updates/Incremental_prediction_20261021"
PRIVATE <- file.path(OUT,"private_patient_level")
dir.create(file.path(OUT,"diagnostics"),showWarnings=FALSE)
SEED <- 20261021L
B <- 1000L
K <- 10L
REPEATS <- 10L
THRESHOLDS <- seq(.05,.50,.01)
DBS <- c("MIMIC","eICU","SICdb")
BASE_MODELS <- c("SOFA","SOFA_plus_class","SOFA_plus_score","SOFA_plus_score_plus_class")
all_apparent <- list(); all_validated <- list(); all_repeat <- list()
all_calibration <- list(); all_dca <- list(); all_coefficients <- list()
all_optimism <- list(); all_auc_diagnostics <- list(); fit_failures <- list()

# DeLong placements for a set of paired continuous scores. Midranks explicitly
# preserve ties, which matter when SOFA is the only risk-model covariate.
delong <- function(y,P) {
  if(is.null(dim(P))) P <- matrix(P,ncol=1)
  n1 <- sum(y==1); n0 <- sum(y==0)
  v10 <- v01 <- vector("list",ncol(P))
  for(j in seq_len(ncol(P))) {
    ra <- rank(P[,j],ties.method="average")
    v10[[j]] <- (ra[y==1]-rank(P[y==1,j],ties.method="average"))/n0
    v01[[j]] <- 1-(ra[y==0]-rank(P[y==0,j],ties.method="average"))/n1
  }
  v10 <- do.call(cbind,v10); v01 <- do.call(cbind,v01)
  list(auc=colMeans(v10), covariance=cov(v10)/n1+cov(v01)/n0)
}
auc_rank <- function(y,p) {
  n1<-sum(y==1);n0<-sum(y==0)
  (sum(rank(p,ties.method="average")[y==1])-n1*(n1+1)/2)/(n1*n0)
}
increment <- function(y,p0,p1) {
  change <- p1-p0
  s <- ifelse(change>1e-12,1,ifelse(change< -1e-12,-1,0))
  event <- mean(s[y==1]); nonevent <- -mean(s[y==0])
  c(delta_auc=auc_rank(y,p1)-auc_rank(y,p0),
    continuous_nri=event+nonevent,nri_events=event,nri_nonevents=nonevent,
    idi=(mean(p1[y==1])-mean(p1[y==0]))-(mean(p0[y==1])-mean(p0[y==0])))
}
fit_risk <- function(X,y,weights=NULL) {
  if(is.null(weights)) weights<-rep(1,length(y))
  fit <- suppressWarnings(glm.fit(X,y,weights=weights,family=binomial(),control=glm.control(maxit=50)))
  if(!fit$converged || any(!is.finite(fit$coefficients)))
    stop("Logistic fit failed or non-finite coefficients")
  fit
}
predict_risk <- function(fit,X) pmin(pmax(plogis(drop(X%*%fit$coefficients)),1e-8),1-1e-8)
make_designs <- function(d,full=FALSE) {
  basic <- cbind(intercept=1,SOFA_per5=d$sofa/5)
  clas <- cbind(C2=as.integer(d$class==2),C3=as.integer(d$class==3))
  result <- list(SOFA=basic,SOFA_plus_class=cbind(basic,clas))
  if(!full) {
    combined <- cbind(basic,additional_score_per50=d$additional_score/50)
    result$SOFA_plus_score <- combined
    result$SOFA_plus_score_plus_class <- cbind(combined,clas)
  }
  result
}
stratified_folds <- function(y,k) {
  fold <- integer(length(y))
  for(g in 0:1) {
    idx<-which(y==g);fold[idx]<-sample(rep(seq_len(k),length.out=length(idx)))
  }
  fold
}
calibration_metrics <- function(y,p) {
  lp<-qlogis(p)
  joint<-glm(y~lp,family=binomial())
  citl<-glm(y~1+offset(lp),family=binomial())
  cs<-summary(joint)$coefficients
  ci<-summary(citl)$coefficients
  c(brier=mean((y-p)^2),calibration_in_large=coef(citl)[1],
    calibration_in_large_low=coef(citl)[1]-1.96*ci[1,2],
    calibration_in_large_high=coef(citl)[1]+1.96*ci[1,2],
    calibration_intercept=coef(joint)[1],
    calibration_slope=coef(joint)[2],
    calibration_slope_low=coef(joint)[2]-1.96*cs[2,2],
    calibration_slope_high=coef(joint)[2]+1.96*cs[2,2])
}
wilson <- function(events,n) {
  z<-1.96;prop<-events/n;den<-1+z^2/n
  centre<-(prop+z^2/(2*n))/den
  half<-z*sqrt(prop*(1-prop)/n+z^2/(4*n^2))/den
  c(low=centre-half,high=centre+half)
}
calibration_bins <- function(y,p) {
  # Quantile breakpoints keep tied predictions together rather than splitting
  # patients with exactly the same SOFA-based prediction arbitrarily.
  br<-unique(quantile(p,probs=seq(0,1,.1),names=FALSE))
  bins<-cut(p,breaks=c(-Inf,br[-c(1,length(br))],Inf),include.lowest=TRUE,labels=FALSE)
  do.call(rbind,lapply(sort(unique(bins)),function(g){
    z<-which(bins==g);ci<-wilson(sum(y[z]),length(z))
    data.frame(bin=g,n=length(z),mean_predicted=mean(p[z]),observed=mean(y[z]),
               observed_low=ci[1],observed_high=ci[2])
  }))
}
net_benefit <- function(y,p,t) {
  positive<-p>=t
  mean(positive & y==1)-mean(positive & y==0)*t/(1-t)
}

# Verify the placement implementation against pairwise definitions, including ties.
test_y<-c(1,1,1,0,0,0,0)
test_P<-cbind(c(.1,.4,.4,.1,.2,.4,.8),c(.2,.3,.6,.1,.3,.5,.9))
test_dl<-delong(test_y,test_P)
for(j in 1:2) {
  brute<-mean(outer(test_P[test_y==1,j],test_P[test_y==0,j],
                    function(a,b)as.numeric(a>b)+.5*as.numeric(a==b)))
  stopifnot(abs(brute-test_dl$auc[j])<1e-12)
}
writeLines(c(capture.output(sessionInfo()),
             paste("Seed",SEED,"Bootstrap",B,"Folds",K,"Repeats",REPEATS)),
           file.path(OUT,"diagnostics","session_info.txt"))

for(di in seq_along(DBS)) {
  db<-DBS[di]
  d<-read.csv(file.path(PRIVATE,paste0(db,"_complete_score_input.csv")),check.names=FALSE)
  y<-as.integer(d$y); n<-nrow(d)
  designs<-make_designs(d); names_m<-names(designs)
  full_fits<-lapply(designs,fit_risk,y=y)
  apparent<-sapply(names_m,function(m)predict_risk(full_fits[[m]],designs[[m]]))
  colnames(apparent)<-names_m
  dl<-delong(y,apparent)
  stopifnot(max(abs(dl$auc-apply(apparent,2,function(p)auc_rank(y,p))))<1e-12)
  for(m in names_m) {
    all_coefficients[[length(all_coefficients)+1]]<-data.frame(dataset=db,model=m,
         term=names(full_fits[[m]]$coefficients),coefficient=full_fits[[m]]$coefficients)
  }
  # Paired patient bootstrap: each sample refits all four risk models.
  # Percentile confidence intervals describe the apparent incremental metrics.
  set.seed(SEED+1000L*di)
  boot<-array(NA_real_,dim=c(B,2,5),
              dimnames=list(NULL,c("SOFA","SOFA_plus_score"),names(increment(y,apparent[,1],apparent[,2]))))
  optimism<-matrix(NA_real_,nrow=B,ncol=length(names_m),dimnames=list(NULL,names_m))
  for(b in seq_len(B)) {
    idx<-sample.int(n,n,replace=TRUE)
    fits<-lapply(designs,function(X)fit_risk(X[idx,,drop=FALSE],y[idx]))
    train_p<-sapply(names_m,function(m)predict_risk(fits[[m]],designs[[m]][idx,,drop=FALSE]))
    test_p<-sapply(names_m,function(m)predict_risk(fits[[m]],designs[[m]]))
    for(j in 1:2) {
      a<-if(j==1)1 else 3
      boot[b,j,]<-increment(y[idx],train_p[,a],train_p[,a+1])
    }
    optimism[b,]<-vapply(seq_along(names_m),function(j)
      auc_rank(y[idx],train_p[,j])-auc_rank(y,test_p[,j]),numeric(1))
    if(b%%100==0) {
      cat(sprintf("%s paired bootstrap %d/%d\n",db,b,B));flush.console()
    }
  }
  saveRDS(list(incremental_metrics=boot,auc_optimism=optimism),
          file.path(PRIVATE,paste0(db,"_bootstrap_replicates.rds")))
  for(j in 1:2) {
    a<-if(j==1)1 else 3; b<-a+1
    inc<-increment(y,apparent[,a],apparent[,b])
    se<-sqrt(dl$covariance[a,a]+dl$covariance[b,b]-2*dl$covariance[a,b])
    z<-inc["delta_auc"]/se
    ci<-apply(boot[,j,,drop=FALSE],3,quantile,probs=c(.025,.975),na.rm=TRUE)
    row<-data.frame(dataset=db,comparison=if(j==1)"SOFA" else "SOFA + additional score",
          additional_score=if(db=="SICdb")"SAPS 3" else "APS",n=n,events=sum(y),
          auc_base=dl$auc[a],auc_base_low=max(0,dl$auc[a]-1.96*sqrt(dl$covariance[a,a])),
          auc_base_high=min(1,dl$auc[a]+1.96*sqrt(dl$covariance[a,a])),
          auc_augmented=dl$auc[b],auc_augmented_low=max(0,dl$auc[b]-1.96*sqrt(dl$covariance[b,b])),
          auc_augmented_high=min(1,dl$auc[b]+1.96*sqrt(dl$covariance[b,b])),
          delta_auc=inc["delta_auc"],delta_auc_low=inc["delta_auc"]-1.96*se,
          delta_auc_high=inc["delta_auc"]+1.96*se,p_delong=2*pnorm(-abs(z)),
          continuous_nri=inc["continuous_nri"],nri_low=ci[1,"continuous_nri"],nri_high=ci[2,"continuous_nri"],
          nri_events=inc["nri_events"],nri_nonevents=inc["nri_nonevents"],
          idi=inc["idi"],idi_low=ci[1,"idi"],idi_high=ci[2,"idi"],
          auc_base_bootstrap_corrected=dl$auc[a]-mean(optimism[,a]),
          auc_augmented_bootstrap_corrected=dl$auc[b]-mean(optimism[,b]),
          delta_auc_bootstrap_corrected=inc["delta_auc"]-mean(optimism[,b]-optimism[,a]))
    rownames(row)<-NULL;all_apparent[[length(all_apparent)+1]]<-row
  }
  # Repeated stratified 10-fold CV. The held-out outcome is never used in model
  # fitting or preprocessing. All four models use identical folds.
  cv<-array(NA_real_,dim=c(n,length(names_m),REPEATS),dimnames=list(NULL,names_m,NULL))
  folds<-matrix(NA_integer_,n,REPEATS)
  for(rep in seq_len(REPEATS)) {
    set.seed(SEED+100000L*di+rep)
    fold<-stratified_folds(y,K);folds[,rep]<-fold
    for(k in seq_len(K)) {
      tr<-fold!=k;te<-!tr
      for(j in seq_along(names_m)) {
        X<-designs[[j]]
        fit<-fit_risk(X[tr,,drop=FALSE],y[tr])
        cv[te,j,rep]<-predict_risk(fit,X[te,,drop=FALSE])
      }
    }
    for(j in seq_along(names_m))
      all_repeat[[length(all_repeat)+1]]<-data.frame(dataset=db,repeat_id=rep,model=names_m[j],
                          auc=auc_rank(y,cv[,j,rep]),brier=mean((y-cv[,j,rep])^2))
  }
  stopifnot(all(is.finite(cv)))
  oof<-apply(cv,c(1,2),mean)
  saveRDS(list(id=d$id,outcome=y,folds=folds,predictions_by_repeat=cv,mean_predictions=oof),
          file.path(PRIVATE,paste0(db,"_OOF_predictions_and_folds.rds")))
  write.csv(data.frame(id=d$id,y=y,class=d$class,oof),
            file.path(PRIVATE,paste0(db,"_OOF_predictions.csv")),row.names=FALSE)
  for(j in seq_along(names_m)) {
    m<-names_m[j];p<-oof[,j];cm<-calibration_metrics(y,p);names(cm)<-sub("\\..*$","",names(cm))
    all_validated[[length(all_validated)+1]]<-data.frame(dataset=db,model=m,n=n,events=sum(y),
                   auc=auc_rank(y,p),as.list(cm),row.names=NULL,check.names=FALSE)
    cb<-calibration_bins(y,p);cb$dataset<-db;cb$model<-m
    all_calibration[[length(all_calibration)+1]]<-cb
  }
  for(j in 1:2) {
    a<-if(j==1)1 else 3;b<-a+1
    for(t in THRESHOLDS) {
      nb0<-net_benefit(y,oof[,a],t);nb1<-net_benefit(y,oof[,b],t)
      all_dca[[length(all_dca)+1]]<-data.frame(dataset=db,
        comparison=if(j==1)"SOFA" else "SOFA + additional score",threshold=t,
        nb_base=nb0,nb_augmented=nb1,nb_all=mean(y)-(1-mean(y))*t/(1-t),nb_none=0,
        delta_nb=nb1-nb0,added_net_true_positives_per100=100*(nb1-nb0))
    }
  }
  cat(db," complete; apparent and OOF AUC:\n")
  print(data.frame(model=names_m,apparent=dl$auc,OOF=apply(oof,2,function(p)auc_rank(y,p))))
  write.csv(do.call(rbind,all_apparent),file.path(OUT,"Apparent_incremental_performance.csv"),row.names=FALSE)
  write.csv(do.call(rbind,all_validated),file.path(OUT,"OOF_performance_and_calibration.csv"),row.names=FALSE)
  flush.console()
}

# Coverage check only: SOFA comparisons in all 5,709 evaluable eICU patients.
d<-read.csv(file.path(PRIVATE,"eICU_full_SOFA_input.csv"),check.names=FALSE)
y<-as.integer(d$y);X<-make_designs(d,full=TRUE)
fit<-lapply(X,fit_risk,y=y)
P<-sapply(seq_along(X),function(j)predict_risk(fit[[j]],X[[j]]))
cv<-array(NA_real_,c(nrow(d),2,REPEATS))
for(rep in seq_len(REPEATS)) {
  set.seed(SEED+900000L+rep);fold<-stratified_folds(y,K)
  for(k in seq_len(K))for(j in 1:2) {
    tr<-fold!=k;te<-!tr
    cv[te,j,rep]<-predict_risk(fit_risk(X[[j]][tr,,drop=FALSE],y[tr]),X[[j]][te,,drop=FALSE])
  }
}
oof<-apply(cv,c(1,2),mean);dl<-delong(y,P)
se<-sqrt(dl$covariance[1,1]+dl$covariance[2,2]-2*dl$covariance[1,2])
coverage_check<-data.frame(dataset="eICU full MAKE30 cohort",n=nrow(d),events=sum(y),
                           auc_base=dl$auc[1],auc_augmented=dl$auc[2],delta_auc=diff(dl$auc),
                           p_delong=2*pnorm(-abs(diff(dl$auc)/se)),
                           oof_auc_base=auc_rank(y,oof[,1]),oof_auc_augmented=auc_rank(y,oof[,2]))
write.csv(coverage_check,file.path(OUT,"eICU_full_SOFA_coverage_check.csv"),row.names=FALSE)
write.csv(do.call(rbind,all_repeat),file.path(OUT,"CV_repeat_metrics.csv"),row.names=FALSE)
write.csv(do.call(rbind,all_calibration),file.path(OUT,"OOF_calibration_bins.csv"),row.names=FALSE)
write.csv(do.call(rbind,all_dca),file.path(OUT,"OOF_decision_curve_values.csv"),row.names=FALSE)
write.csv(do.call(rbind,all_coefficients),file.path(OUT,"Apparent_logistic_model_coefficients.csv"),row.names=FALSE)
final_apparent<-do.call(rbind,all_apparent)
final_apparent$q_BH_six_AUC_comparisons<-p.adjust(final_apparent$p_delong,method="BH")
write.csv(final_apparent,file.path(OUT,"Apparent_incremental_performance.csv"),row.names=FALSE)
writeLines("All logistic fits converged. No missing score was imputed. All paired models used identical patients and CV folds.",
           file.path(OUT,"diagnostics","completion.txt"))
cat("ANALYSIS COMPLETE\n")
