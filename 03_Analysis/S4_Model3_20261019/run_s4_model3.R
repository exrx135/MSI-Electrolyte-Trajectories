if (.Platform$OS.type=='windows') Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8')
suppressPackageStartupMessages({library(mice);library(sandwich);library(jsonlite)})
options(warn=1)
root<-'Z:/Research/MSI-1'; ar<-file.path(root,'历史归档/results_整理前_20261019/results')
out<-file.path(root,'analysis_updates/S4_Model3_20261019');dir.create(out,recursive=TRUE,showWarnings=FALSE)
db<-commandArgs(TRUE)[1];stopifnot(db %in% c('MIMIC','eICU','SICdb'))
cp<-file.path(out,paste0(db,'_checkpoints'));dir.create(cp,showWarnings=FALSE)
read<-function(p)read.csv(p,check.names=FALSE,stringsAsFactors=FALSE,na.strings=c('','NA','NULL'))
spec<-fromJSON(file.path(root,'results/03_Final_Code/model_specification.json'));M3<-spec$model3
frozen_state<-readRDS(file.path(root,'RDS/V318_frozen_scoring_state_G3.rds'))
Z<-c('msi_z','sodium_z','potassium_z','chloride_z','bicarbonate_z')
rawcols<-c(paste0('msi_m',1:12),unlist(lapply(c('sodium','potassium','chloride','bicarbonate'),function(v)paste0(v,'_t',1:4))))
varindex<-rep(1:5,c(12,4,4,4,4))
scenarios<-rbind(MAR=rep(0,5),MNAR_down025=rep(-.25,5),MNAR_down050=rep(-.5,5),MNAR_up025=rep(.25,5),MNAR_up050=rep(.5,5),MNAR_pattern025=c(.25,0,.25,.25,-.25),MNAR_pattern050=c(.5,0,.5,.5,-.5),MNAR_reverse025=c(-.25,0,-.25,-.25,.25),MNAR_reverse050=c(-.5,0,-.5,-.5,.5));colnames(scenarios)<-Z
makepanel<-function(d){n<-nrow(d);p<-matrix(NA_real_,n*12,5,dimnames=list(NULL,Z));p[,1]<-as.vector(t(as.matrix(d[rawcols[1:12]])));ii<-rep(seq_len(n),each=4);tt<-rep(c(2,5,8,11),n);rr<-(ii-1)*12+tt
 for(j in 2:5)p[rr,j]<-as.vector(t(as.matrix(d[rawcols[varindex==j]])))
 sweep(sweep(p,2,frozen_state$raw_center[Z],'-'),2,frozen_state$raw_scale[Z],'/')}
complete_fast<-function(p,k){n<-nrow(p);tt<-rep(seq_len(12),n/12);mu<-frozen_state$Xmat%*%frozen_state$completion_beta[[k]];S<-frozen_state$completion_Sigma[[k]];codes<-as.vector(is.na(p)%*%(2^(0:4)));ans<-p
 for(code in unique(codes)){rr<-which(codes==code);mis<-which(as.logical(intToBits(as.integer(code))[1:5]));if(!length(mis))next;ob<-setdiff(1:5,mis);v<-mu[tt[rr],mis,drop=FALSE];if(length(ob))v<-v+(p[rr,ob,drop=FALSE]-mu[tt[rr],ob,drop=FALSE])%*%solve(S[ob,ob,drop=FALSE],S[ob,mis,drop=FALSE]);ans[rr,mis]<-v};ans}
score<-function(p){n<-nrow(p)/12;tt<-rep(1:12,n);scores<-matrix(NA_real_,n,3);for(k in 1:3){v<-complete_fast(p,k);S<-frozen_state$score_Sigma[[k]];dif<-v-(frozen_state$Xmat%*%frozen_state$score_beta[[k]])[tt,,drop=FALSE];den<--2.5*log(2*pi)-.5*log(det(S))-.5*rowSums((dif%*%solve(S))*dif);scores[,k]<-colSums(matrix(den,12,n))+log(frozen_state$prior[k])}
 bound<-getFromNamespace('pldBound','gbmt');scores[]<-vapply(as.vector(scores),bound,numeric(1));pr<-exp(scores-apply(scores,1,max));pr<-pr/rowSums(pr);list(class=max.col(pr,ties.method='first'),prob=pr)}
ari<-function(a,b){x<-table(a,b);ch<-function(z)z*(z-1)/2;A<-sum(ch(rowSums(x)));B<-sum(ch(colSums(x)));E<-A*B/ch(sum(x));(sum(ch(x))-E)/((A+B)/2-E)}

rawfile<-c(MIMIC='MIMIC_FINAL_MAKE30_ELECTROLYTEQC.csv',eICU='eICU_FINAL_RAW_ELECTROLYTEQC.csv',SICdb='SICdb_FINAL_MAKE30_ELECTROLYTEQC.csv')[[db]]
id<-c(MIMIC='stay_id',eICU='patientunitstayid',SICdb='case_id')[[db]]
d<-read(file.path(root,rawfile));a<-read(file.path(ar,'关联分析_三模型_MICE_RR/data',paste0(db,'_analysis_input.csv')))
ai<-match(a$cohort_unit_id,d[[id]]);stopifnot(!anyNA(ai),!anyDuplicated(a$cohort_unit_id),all(a$make30==d$make30[ai]))
statehash<-unname(tools::md5sum(file.path(root,'RDS/V318_frozen_scoring_state_G3.rds')));stopifnot(statehash=='89f684d00680fafc30577d677fee6a30')
ref<-score(makepanel(d))$class;stopifnot(all(ref[ai]==a$assigned_class))
traj<-readRDS(file.path(ar,'V331_MAR_MNAR敏感性_G3',paste0('V331_',db,'_MICE20_G3.rds')))
stopifnot(traj$m==20,identical(traj$data[rawcols],d[rawcols]))
covimp<-readRDS(file.path(ar,'关联分析_三模型_MICE_RR/diagnostics',paste0(db,'_mice_final.rds')))
stopifnot(covimp$m==30,nrow(covimp$data)==nrow(a),all(covimp$data$make30==a$make30),all(as.character(covimp$data$phenotype)==paste0('C',a$assigned_class)))
mask<-is.na(as.matrix(d[rawcols]));lims<-t(sapply(1:5,function(j)range(as.matrix(d[rawcols[varindex==j]]),na.rm=TRUE)))
oldstat<-read(file.path(ar,'V331_MAR_MNAR敏感性_G3/V331_分类逐次结果_G3.csv'));oldstat<-oldstat[oldstat$dataset==db,]
fitrr<-function(dat){
 continuous<-M3[!vapply(dat[,M3,drop=FALSE],is.factor,logical(1))]
 for(v in continuous)dat[[v]]<-(dat[[v]]-mean(covimp$data[[v]],na.rm=TRUE))/sd(covimp$data[[v]],na.rm=TRUE)
 f<-glm(reformulate(c('phenotype',M3),'make30'),dat,family=poisson('log'),control=glm.control(maxit=100,epsilon=1e-9),na.action=na.fail)
 stopifnot(f$converged,!anyNA(coef(f)),nobs(f)==nrow(a))
 V<-vcovHC(f,type='HC1');ix<-match(c('phenotypeC2','phenotypeC3'),names(coef(f)))
 stopifnot(all(is.finite(V)),all(diag(V)>0))
 data.frame(comparison=c('C2 vs C1','C3 vs C1'),Q=as.numeric(coef(f)[ix]),U=diag(V)[ix],n=nobs(f),k=length(coef(f)),max_fitted=max(fitted(f)))
}
all<-list();stats<-list()
for(j in 1:30){z<-fitrr(complete(covimp,j));z$dataset<-db;z$scenario<-'Original';z$draw<-j;all[[length(all)+1]]<-z}
cat(db,'original 30 Model 3 fits complete',format(Sys.time()),'\n');flush.console()
for(j in 1:20){baseval<-as.matrix(complete(traj,j)[rawcols]);stopifnot(all(baseval[!mask]==as.matrix(d[rawcols])[!mask]))
 for(si in seq_len(nrow(scenarios))){sn<-rownames(scenarios)[si];file<-file.path(cp,paste0(sn,'_',sprintf('%02d',j),'.rds'))
  if(file.exists(file)){res<-readRDS(file);stopifnot(identical(res$covariates,M3));all[[length(all)+1]]<-res$fit;stats[[length(stats)+1]]<-res$stats;next}
  val<-baseval
  for(v in 1:28){ix<-mask[,v];vv<-varindex[v];nv<-val[ix,v]+scenarios[sn,vv]*frozen_state$raw_scale[Z[vv]];val[ix,v]<-pmin(lims[vv,2],pmax(lims[vv,1],nv))}
  dd<-d;dd[rawcols]<-val;cl<-score(makepanel(dd))$class
  st<-data.frame(dataset=db,scenario=sn,draw=j,agreement=mean(cl==ref),ARI=ari(ref,cl))
  os<-oldstat[oldstat$scenario==sn & oldstat$draw==j,];stopifnot(nrow(os)==1,abs(os$agreement-st$agreement)<1e-12,abs(os$ARI-st$ARI)<1e-10)
  cd<-covimp$data;cd$phenotype<-factor(paste0('C',cl[ai]),levels=c('C1','C2','C3'))
  init<-complete(covimp,j);init$phenotype<-cd$phenotype
  seed<-202610190+match(db,c('MIMIC','eICU','SICdb'))*10000+si*100+j
  # One independent conditional covariate draw per saved trajectory draw.
  # No recycling or Cartesian product of 20 and 30 imputations.
  mi<-mice(cd,m=1,maxit=30,method=covimp$method,predictorMatrix=covimp$predictorMatrix,data.init=init,seed=seed,donors=5,ridge=1e-5,printFlag=FALSE)
  cc<-complete(mi,1);stopifnot(!anyNA(cc))
  for(v in names(cd)){ix<-!is.na(cd[[v]]);stopifnot(all(as.character(cc[[v]][ix])==as.character(cd[[v]][ix])))}
  fit<-fitrr(cc);fit$dataset<-db;fit$scenario<-sn;fit$draw<-j
  res<-list(fit=fit,stats=st,covariates=M3,seed=seed,imputed_values=mi$imp,chainMean=mi$chainMean,chainVar=mi$chainVar,loggedEvents=mi$loggedEvents,temporary_class=cl,analysis_ids=a$cohort_unit_id)
  saveRDS(res,file);all[[length(all)+1]]<-fit;stats[[length(stats)+1]]<-st
 }
 cat(db,'draw',j,'of 20 complete',format(Sys.time()),'\n');flush.console()
 write.csv(do.call(rbind,all),file.path(out,paste0(db,'_coefficients.csv')),row.names=FALSE)
}
rr<-do.call(rbind,all)
pooled<-do.call(rbind,lapply(split(rr,interaction(rr$scenario,rr$comparison,drop=TRUE)),function(x){p<-pool.scalar(x$Q,x$U,n=x$n[1],k=x$k[1]);se<-sqrt(p$t);cr<-qt(.975,p$df);data.frame(dataset=db,scenario=x$scenario[1],comparison=x$comparison[1],RR=exp(p$qbar),lower=exp(p$qbar-cr*se),upper=exp(p$qbar+cr*se),SE_logRR=se,df=p$df,m=nrow(x),MCSE_logRR=sqrt(p$b/nrow(x)),n=x$n[1],events=sum(a$make30))}))
main<-read(file.path(ar,'关联分析_三模型_MICE_RR/pooled_RR_results.csv'));main<-main[main$dataset==db & main$model=='Model 3',]
z<-merge(pooled[pooled$scenario=='Original',],main,by='comparison',suffixes=c('_new','_main'))
stopifnot(max(abs(z$RR_new-z$RR_main))<1e-8,max(abs(z$lower-z$CI95_low))<1e-8,max(abs(z$upper-z$CI95_high))<1e-8)
write.csv(pooled,file.path(out,paste0(db,'_pooled.csv')),row.names=FALSE);write.csv(do.call(rbind,stats),file.path(out,paste0(db,'_classification_check.csv')),row.names=FALSE)
cat(db,'COMPLETE: original reproduces published Model 3; all classification draws reproduce old sensitivity.',format(Sys.time()),'\n')

