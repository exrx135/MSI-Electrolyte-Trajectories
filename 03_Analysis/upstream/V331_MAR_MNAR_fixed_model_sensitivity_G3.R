# V331: complete MAR/MNAR sensitivity; fixed G3 classifier, temporary assignments, archived MAKE30 regression specification.
main <- function() {
root <- 'Z:/Research/MSI-1'
suppressPackageStartupMessages(library(ggplot2))
out <- file.path(root,'results','V331_MAR_MNAR敏感性_G3');dir.create(out,recursive=TRUE,showWarnings=FALSE)
fn <- function(n,e='csv')file.path(out,paste0('V331_',n,'_G3.',e))
read <- function(p)read.csv(p,check.names=FALSE,stringsAsFactors=FALSE,na.strings=c('','NA','NULL'))
csv <- function(x,n)write.csv(x,fn(n),row.names=FALSE,na='',fileEncoding='UTF-8')
raws<-c(MIMIC='MIMIC_FINAL_MAKE30_ELECTROLYTEQC.csv',eICU='eICU_FINAL_RAW_ELECTROLYTEQC.csv',SICdb='SICdb_FINAL_MAKE30_ELECTROLYTEQC.csv')
cls<-c(MIMIC='V318_MIMIC_formal_classification_G3.csv',eICU='V318_eICU_fixed_scoring_G3.csv',SICdb='V318_SICdb_fixed_scoring_G3.csv')
ids<-c(MIMIC='stay_id',eICU='patientunitstayid',SICdb='case_id')
paths<-c(file.path(root,raws),file.path(root,'STEP5_冻结分型与预后外部验证','冻结分类',cls))
md5<-unname(tools::md5sum(paths));base<-read(file.path(root,'V325_文件总清单_GALL.csv'))
stopifnot(all(md5==base$MD5[match(substring(paths,nchar(root)+2),base$文件)]))
cohorts<-list();counts<-list()
for(db in names(raws)) {
 d<-read(file.path(root,raws[[db]]));cl<-read(file.path(root,'STEP5_冻结分型与预后外部验证','冻结分类',cls[[db]]))
 ix<-match(d[[ids[[db]]]],cl$unit_id);stopifnot(!anyNA(ix),!anyDuplicated(d[[ids[[db]]]]),!anyDuplicated(cl$unit_id),nrow(d)==nrow(cl),all(cl$assigned_class%in%1:3))
 d$phenotype<-paste0('C',cl$assigned_class[ix]);cohorts[[db]]<-d
 counts[[db]]<-data.frame(dataset=db,phenotype=paste0('C',1:3),n=as.integer(table(factor(d$phenotype,levels=paste0('C',1:3)))))
}

suppressPackageStartupMessages(library(mice))
statepath<-file.path(root,'RDS','V318_frozen_scoring_state_G3.rds');stopifnot(unname(tools::md5sum(statepath))=='89f684d00680fafc30577d677fee6a30');state<-readRDS(statepath)
Z<-c('msi_z','sodium_z','potassium_z','chloride_z','bicarbonate_z');G<-3L
rawcols<-c(paste0('msi_m',1:12),unlist(lapply(c('sodium','potassium','chloride','bicarbonate'),function(v)paste0(v,'_t',1:4))))
varindex<-rep(1:5,c(12,4,4,4,4));CORE<-c('age','sex','sepsis_sofa_score','baseline_egfr_mdrd','ckd_or_renal_disease','diabetes','mechanical_ventilation_0_24h','vasopressor_0_24h')
scenarios<-rbind(MAR=rep(0,5),MNAR_down025=rep(-.25,5),MNAR_down050=rep(-.5,5),MNAR_up025=rep(.25,5),MNAR_up050=rep(.5,5),MNAR_pattern025=c(.25,0,.25,.25,-.25),MNAR_pattern050=c(.5,0,.5,.5,-.5),MNAR_reverse025=c(-.25,0,-.25,-.25,.25),MNAR_reverse050=c(-.5,0,-.5,-.5,.5));colnames(scenarios)<-Z
csv(data.frame(scenario=rownames(scenarios),scenarios,check.names=FALSE),'情景标准差偏移');csv(data.frame(scenario=rownames(scenarios),sweep(scenarios,2,state$raw_scale[Z],'*'),check.names=FALSE),'情景原单位偏移')
logmsg<-function(...)cat(format(Sys.time()),paste(...,collapse=' '),'\n')
makepanel<-function(d){n<-nrow(d);p<-matrix(NA_real_,n*12,5,dimnames=list(NULL,Z));p[,1]<-as.vector(t(as.matrix(d[rawcols[1:12]])));ii<-rep(seq_len(n),each=4);tt<-rep(c(2,5,8,11),n);rr<-(ii-1)*12+tt
 for(j in 2:5)p[rr,j]<-as.vector(t(as.matrix(d[rawcols[varindex==j]])))
 sweep(sweep(p,2,state$raw_center[Z],'-'),2,state$raw_scale[Z],'/')}
complete_fast<-function(p,k){n<-nrow(p);tt<-rep(seq_len(12),n/12);mu<-state$Xmat%*%state$completion_beta[[k]];S<-state$completion_Sigma[[k]];codes<-as.vector(is.na(p)%*%(2^(0:4)));ans<-p
 for(code in unique(codes)){rr<-which(codes==code);mis<-which(as.logical(intToBits(as.integer(code))[1:5]));if(!length(mis))next;ob<-setdiff(1:5,mis);v<-mu[tt[rr],mis,drop=FALSE];if(length(ob))v<-v+(p[rr,ob,drop=FALSE]-mu[tt[rr],ob,drop=FALSE])%*%solve(S[ob,ob,drop=FALSE],S[ob,mis,drop=FALSE]);ans[rr,mis]<-v};ans}
score<-function(p){n<-nrow(p)/12;tt<-rep(1:12,n);scores<-matrix(NA_real_,n,3);for(k in 1:3){v<-complete_fast(p,k);S<-state$score_Sigma[[k]];dif<-v-(state$Xmat%*%state$score_beta[[k]])[tt,,drop=FALSE];den<--2.5*log(2*pi)-.5*log(det(S))-.5*rowSums((dif%*%solve(S))*dif);scores[,k]<-colSums(matrix(den,12,n))+log(state$prior[k])}
 bound<-getFromNamespace('pldBound','gbmt');scores[]<-vapply(as.vector(scores),bound,numeric(1));pr<-exp(scores-apply(scores,1,max));pr<-pr/rowSums(pr);list(class=max.col(pr,ties.method='first'),prob=pr)}
ari<-function(a,b){x<-table(a,b);ch<-function(z)z*(z-1)/2;A<-sum(ch(rowSums(x)));B<-sum(ch(colSums(x)));E<-A*B/ch(sum(x));(sum(ch(x))-E)/((A+B)/2-E)}
fitrr<-function(a,cl){a$phenotype<-factor(paste0('C',cl),levels=paste0('C',1:3));a$sex<-factor(a$sex,levels=c('Female','Male'));f<-glm(reformulate(c('phenotype',CORE),'make30'),data=a,family=poisson('log'),x=TRUE,y=TRUE,na.action=na.fail);stopifnot(f$converged,!anyNA(coef(f)));V<-vcov(f);U<-f$x*as.numeric(f$y-fitted(f));V<-V%*%crossprod(U)%*%V*nrow(f$x)/(nrow(f$x)-ncol(f$x));ix<-match(c('phenotypeC2','phenotypeC3'),names(coef(f)));list(Q=coef(f)[ix],U=diag(V)[ix],n=nrow(a),p=ncol(f$x))}
afiles<-c(MIMIC='MIMIC_预后关联/V319_MIMIC_analysis_cohort_G3.csv',eICU='eICU_外部验证/V321_analysis_cohort_G3.csv',SICdb='SICdb_外部验证/V322_analysis_cohort_G3.csv');sexdraw<-read(file.path(root,'STEP5_冻结分型与预后外部验证/eICU_外部验证/V321_sex_MICE_20x_G3.csv'))
allstats<-allrr<-alltrans<-alldesc<-allclassmeans<-list();allpool<-list();replay<-list();diagall<-list()
for(db in names(cohorts)){
 logmsg('START',db);d<-cohorts[[db]];ref<-as.integer(sub('C','',d$phenotype));p0<-makepanel(d);b0<-score(p0);stopifnot(all(b0$class==ref))
 # independent comparison of optimized conditional mean with frozen package implementation
 ss<-seq_len(min(240,nrow(p0)));for(k in 1:3){slow<-getFromNamespace('imputFun','gbmt')(p0[ss,],rep(state$time_model-min(state$time_model)+1,length.out=length(ss)),state$completion_beta[[k]],state$completion_Sigma[[k]],2L);stopifnot(max(abs(slow-complete_fast(p0[ss,],k)))<1e-10)}
 replay[[db]]<-data.frame(dataset=db,n=nrow(d),agreement=mean(ref==b0$class));csv(do.call(rbind,replay),'原分类重演核验')
 input<-d[c(rawcols,CORE)];input$sex<-factor(input$sex);stopifnot(!anyNA(input[CORE]))
 method<-setNames(rep('',ncol(input)),names(input));method[rawcols]<-'pmm';pred<-matrix(0,ncol(input),ncol(input),dimnames=list(names(input),names(input)));pred[rawcols,]<-1;diag(pred)<-0
 csv(data.frame(variable=names(input),method=method,predictors=apply(pred,1,function(x)paste(names(x)[x==1],collapse=';'))),paste0(db,'_插补模型'))
 dest<-fn(paste0(db,'_MICE20'),'rds')
 if(file.exists(dest)){imp<-readRDS(dest);stopifnot(identical(imp$data,input))}else{logmsg('MICE 20 x 20',db);imp<-mice(input,m=20,maxit=20,method=method,predictorMatrix=pred,seed=20261012+match(db,names(cohorts)),donors=5,printFlag=TRUE);saveRDS(imp,dest)}
 rhat<-function(arr){v<-arr[(floor(nrow(arr)/2)+1):nrow(arr),,drop=FALSE];half<-floor(nrow(v)/2);v<-cbind(v[1:half,,drop=FALSE],v[(nrow(v)-half+1):nrow(v),,drop=FALSE]);W<-mean(apply(v,2,var));B<-nrow(v)*var(colMeans(v));if(W==0)return(if(B==0)1 else Inf);sqrt(((nrow(v)-1)/nrow(v)*W+B/nrow(v))/W)}
 calcdiag<-function(imp)do.call(rbind,lapply(rawcols,function(v)data.frame(variable=v,missing_n=sum(is.na(input[[v]])),rhat=if(sum(is.na(input[[v]])))rhat(imp$chainMean[v,,,drop=TRUE]) else NA_real_,sd_rhat=if(sum(is.na(input[[v]])))rhat(sqrt(imp$chainVar[v,,,drop=TRUE])) else NA_real_)))
 dg<-calcdiag(imp)
 minimum_iterations<-c(MIMIC=100,eICU=100,SICdb=60)[[db]];while((imp$iteration<minimum_iterations || any(c(dg$rhat,dg$sd_rhat)>1.1,na.rm=TRUE)) && imp$iteration<160){extra<-if(imp$iteration<minimum_iterations)minimum_iterations-imp$iteration else min(30,160-imp$iteration);logmsg('Extending chains',db,'by',extra);imp<-mice.mids(imp,maxit=extra,printFlag=TRUE);saveRDS(imp,dest);dg<-calcdiag(imp)};if(any(c(dg$rhat,dg$sd_rhat)>1.1,na.rm=TRUE))warning(db,': chain diagnostic remains above 1.1; results provisional')
 dg$dataset<-db;dg$iterations<-imp$iteration;diagall[[db]]<-dg;csv(do.call(rbind,diagall),'链诊断')
 if(!is.null(imp$loggedEvents))csv(imp$loggedEvents,paste0(db,'_MICE事件'))
 trace<-do.call(rbind,lapply(rawcols,function(v)do.call(rbind,lapply(1:20,function(j)data.frame(variable=v,draw=j,iteration=seq_len(imp$iteration),mean=imp$chainMean[v,,j],sd=sqrt(imp$chainVar[v,,j]))))));csv(trace,paste0(db,'_链轨迹'))
 a<-read(file.path(root,'STEP5_冻结分型与预后外部验证',afiles[[db]]));ai<-match(a[[ids[[db]]]],d[[ids[[db]]]]);stopifnot(!anyNA(ai),!anyNA(a$make30));mask<-is.na(as.matrix(d[rawcols]));lims<-t(sapply(1:5,function(j)range(as.matrix(d[rawcols[varindex==j]]),na.rm=TRUE)))
 if(db=='MIMIC' && all(a$sex%in%c('F','M')))a$sex<-ifelse(a$sex=='F','Female','Male')
 if(db=='SICdb' && all(a$sex%in%c('female','male')))a$sex<-ifelse(a$sex=='female','Female','Male')
 for(j in 1:20){filled<-complete(imp,j);stopifnot(!anyNA(filled[rawcols]));baseval<-as.matrix(filled[rawcols]);stopifnot(all(baseval[!mask]==as.matrix(d[rawcols])[!mask]));aj<-a
  if(db=='eICU'){sx<-sexdraw[sexdraw$imputation==j,];ix<-match(sx$unit_id,aj[[ids[[db]]]]);stopifnot(!anyNA(ix));aj$sex[ix]<-sx$imputed_binary_sex}
  for(sn in c('Original',rownames(scenarios))){cl<-ref;clipped<-0L
   if(sn!='Original'){val<-baseval;for(v in 1:28){ix<-mask[,v];jj<-varindex[v];shift<-scenarios[sn,jj]*state$raw_scale[Z[jj]];nv<-val[ix,v]+shift;clipped<-clipped+sum(nv<lims[jj,1]|nv>lims[jj,2]);val[ix,v]<-pmin(lims[jj,2],pmax(lims[jj,1],nv))};dd<-d;dd[rawcols]<-val;cl<-score(makepanel(dd))$class
    if(j==1 || sn=='MAR')for(v in 1:28){ix<-mask[,v];alldesc[[length(alldesc)+1]]<-data.frame(dataset=db,scenario=sn,draw=j,variable=rawcols[v],observed_n=sum(!ix),missing_n=sum(ix),observed_mean=mean(d[[rawcols[v]]],na.rm=TRUE),imputed_mean=if(any(ix))mean(val[ix,v]) else NA_real_,imputed_sd=if(sum(ix)>1)sd(val[ix,v]) else NA_real_)}
   }
   stt<-data.frame(dataset=db,scenario=sn,draw=j,agreement=mean(cl==ref),ARI=ari(ref,cl),clipped_cells=clipped)
   for(k in 1:3){stt[[paste0('n_C',k)]]<-sum(cl==k);stt[[paste0('Jaccard_C',k)]]<-sum(cl==k & ref==k)/sum(cl==k|ref==k)}
   allstats[[length(allstats)+1]]<-stt
   values<-if(sn=='Original')as.matrix(d[rawcols]) else val
   for(k in 1:3){vv<-values[cl==k,,drop=FALSE];allclassmeans[[length(allclassmeans)+1]]<-data.frame(dataset=db,scenario=sn,draw=j,phenotype=paste0('C',k),variable=rawcols,n_class=nrow(vv),observed_n=colSums(!is.na(vv)),mean=colMeans(vv,na.rm=TRUE))}
   tr<-as.data.frame(table(original=factor(ref,levels=1:3),sensitivity=factor(cl,levels=1:3)));tr$dataset<-db;tr$scenario<-sn;tr$draw<-j;alltrans[[length(alltrans)+1]]<-tr
   rr<-fitrr(aj,cl[ai]);allrr[[length(allrr)+1]]<-data.frame(dataset=db,scenario=sn,draw=j,comparison=c('C2 vs C1','C3 vs C1'),logRR=as.numeric(rr$Q),variance=rr$U,n=rr$n,p=rr$p,events=sum(aj$make30))
  };logmsg(db,'draw',j,'complete');csv(do.call(rbind,allstats),'分类逐次结果');csv(do.call(rbind,allrr),'RR逐次结果')
 }
 csv(do.call(rbind,allclassmeans),'敏感性类别逐窗均值');csv(do.call(rbind,alltrans),'分类转移计数');csv(do.call(rbind,alldesc),'插补分布描述')
}
rr<-do.call(rbind,allrr);pooled<-do.call(rbind,lapply(split(rr,interaction(rr$dataset,rr$scenario,rr$comparison,drop=TRUE)),function(x){p<-pool.scalar(x$logRR,x$variance,n=x$n[1],k=x$p[1],rule='rubin1987');se<-sqrt(p$t);crit<-qt(.975,p$df);data.frame(dataset=x$dataset[1],scenario=x$scenario[1],comparison=x$comparison[1],RR=exp(p$qbar),lower=exp(p$qbar-crit*se),upper=exp(p$qbar+crit*se),SE_logRR=se,df=p$df,MCSE_logRR=sqrt(p$b/20),MCSE_over_SE=sqrt(p$b/20)/se,between_variance=p$b,n=x$n[1],events=x$events[1],m=20)}));csv(pooled,'RR合并结果')
stats<-do.call(rbind,allstats);summ<-do.call(rbind,lapply(split(stats,interaction(stats$dataset,stats$scenario,drop=TRUE)),function(x){z<-data.frame(dataset=x$dataset[1],scenario=x$scenario[1]);for(v in setdiff(names(x),c('dataset','scenario','draw'))){z[[paste0(v,'_mean')]]<-mean(x[[v]]);z[[paste0(v,'_min')]]<-min(x[[v]]);z[[paste0(v,'_max')]]<-max(x[[v]])};z}));csv(summ,'分类汇总')
stopifnot(all(md5==unname(tools::md5sum(paths))),unname(tools::md5sum(statepath))=='89f684d00680fafc30577d677fee6a30');csv(data.frame(file=paths,md5_before=md5,md5_after=unname(tools::md5sum(paths))),'输入身份核验');writeLines(capture.output(sessionInfo()),fn('R运行环境','txt'));logmsg('COMPLETE')
}
main()


