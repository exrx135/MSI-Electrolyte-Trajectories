if(.Platform$OS.type=='windows') Sys.setlocale('LC_CTYPE','Chinese (Simplified)_China.utf8')
library(mice);library(sandwich);library(mitml);library(jsonlite)
options(stringsAsFactors=FALSE,warn=1)
fa<-sub('^--file=','',commandArgs(FALSE)[grepl('^--file=',commandArgs(FALSE))]);P<-normalizePath(file.path(dirname(fa),'..'),winslash='/',mustWork=TRUE)
S<-file.path(dirname(P),'关联分析_三模型_MICE_RR');db<-commandArgs(TRUE)[1];stopifnot(db%in%c('MIMIC','eICU','SICdb'))
spec<-fromJSON(file.path(S,'model_specification.json'));vars<-spec$model3
csv<-function(x,nm)write.csv(x,file.path(P,paste0(db,'_',nm,'.csv')),row.names=FALSE,na='')
old<-readRDS(file.path(S,'diagnostics',paste0(db,'_mice_final.rds')));d<-old$data
raw<-read.csv(file.path(S,'data',paste0(db,'_analysis_input.csv')))
stopifnot(nrow(d)==nrow(raw),all(d$make30==raw$make30),all(as.character(d$phenotype)==paste0('C',raw$assigned_class)))
meth<-old$method;meth[meth!='']<-'cart_safe';pred<-old$predictorMatrix
# Correct numeric singleton sampling in installed mice 3.19.0 CART.
# sample(c(27),1) samples 1:27; sample.int selects a donor index instead.
mice.impute.cart_safe <- mice::mice.impute.cart
btext <- paste(deparse(body(mice.impute.cart_safe)),collapse='\n')
btext <- sub('sample(donor[[s]], 1)', 'donor[[s]][sample.int(length(donor[[s]]), 1)]',btext,fixed=TRUE)
# deparse can split the original call across lines; use whitespace-aware match.
btext <- gsub('sample\\(donor\\[\\[s\\]\\],\\s*1\\)', 'donor[[s]][sample.int(length(donor[[s]]), 1)]',btext,perl=TRUE)
stopifnot(grepl('sample.int',btext,fixed=TRUE))
body(mice.impute.cart_safe) <- parse(text=btext)[[1]]
path<-file.path(P,'diagnostics',paste0(db,'_cart_safe_mice.rds'))
cat('START',db,format(Sys.time()),'\n');flush.console()
if(file.exists(path)) {imp<-readRDS(path);stopifnot(identical(imp$data,d),identical(imp$method,meth))} else {
 imp<-mice(d,m=30,maxit=30,method=meth,predictorMatrix=pred,seed=20261019+match(db,c('MIMIC','eICU','SICdb'))*1000,printFlag=FALSE,minbucket=5,cp=1e-4,xval=0)
 saveRDS(imp,path)
}
cat('IMPUTATION DONE',db,format(Sys.time()),'\n');flush.console()
if(!is.null(imp$loggedEvents)) csv(imp$loggedEvents,'imputation_events')
split_rhat<-function(x){x<-tail(x,20);n<-floor(nrow(x)/2);y<-cbind(x[1:n,,drop=FALSE],x[(n+1):(2*n),,drop=FALSE]);W<-mean(apply(y,2,var));B<-n*var(colMeans(y));if(!is.finite(W))return(NA_real_);if(W==0)return(ifelse(B==0,1,Inf));sqrt(((n-1)/n*W+B/n)/W)}
conv<-list()
for(v in names(meth)[meth!='']) for(st in c('mean','sd')) {x<-if(st=='mean')imp$chainMean[v,,] else sqrt(imp$chainVar[v,,]);conv[[length(conv)+1]]<-data.frame(variable=v,statistic=st,split_Rhat=split_rhat(x))}
csv(do.call(rbind,conv),'convergence')
completed<-lapply(1:30,function(i)complete(imp,i))
for(dat in completed)for(v in names(d)){obs<-!is.na(d[[v]]);stopifnot(!anyNA(dat[[v]]),all(as.character(dat[[v]][obs])==as.character(d[[v]][obs])));if(meth[v]=='cart_safe')stopifnot(all(as.character(dat[[v]][!obs])%in%as.character(d[[v]][obs])))}
getgroup<-function(dat,sg){switch(sg,Age=as.numeric(dat$age>=65),Sex=as.numeric(as.character(dat$male)),eGFR=as.numeric(dat$baseline_egfr_mdrd<60),Diabetes=as.numeric(as.character(dat$diabetes)))}
ans<-list();inter<-list();audit<-list();cnt<-list();draws<-list();icount<-list();joint<-list()
for(sg in c('Age','Sex','eGFR','Diabetes')) {
 cat('FIT',db,sg,'\n');flush.console()
 graw<-getgroup(d,sg)
 for(g in c(0,1,NA))for(cl in c('C1','C2','C3')){ix<-if(is.na(g))is.na(graw) else !is.na(graw)&graw==g;ix<-ix&d$phenotype==cl;cnt[[length(cnt)+1]]<-data.frame(subgroup=sg,stratum=ifelse(is.na(g),'Unknown',g),class=cl,n=sum(ix),events=sum(d$make30[ix]))}
 if(length(unique(graw[!is.na(graw)]))<2){
  inter[[length(inter)+1]]<-data.frame(subgroup=sg,F=NA,df1=2,df2=NA,p_interaction=NA)
  for(g in 0:1)for(cl in c('C2','C3')) ans[[length(ans)+1]]<-data.frame(subgroup=sg,stratum=g,comparison=paste(cl,'vs C1'),RR=NA,low=NA,high=NA,p=NA,df=NA,FMI=NA,MCSE_SE=NA)
  cat('NOT ESTIMABLE',db,sg,'only one observed stratum\n');next
 }
 Q<-list();U<-list();qn<-list();un<-list()
 for(j in 1:30){dat<-completed[[j]];dat$sg<-getgroup(dat,sg)
  for(g in 0:1)for(cl in c('C1','C2','C3')){ix<-dat$sg==g&dat$phenotype==cl;icount[[length(icount)+1]]<-data.frame(subgroup=sg,imputation=j,stratum=g,class=cl,n=sum(ix),events=sum(dat$make30[ix]))}
  adj<-setdiff(vars,switch(sg,Sex='male',Diabetes='diabetes',character()))
  continuous<-adj[!vapply(dat[,adj,drop=FALSE],is.factor,logical(1))]
  for(v in continuous)dat[[v]]<-(dat[[v]]-mean(d[[v]],na.rm=TRUE))/sd(d[[v]],na.rm=TRUE)
  dat$int2<-as.numeric(dat$phenotype=='C2')*dat$sg;dat$int3<-as.numeric(dat$phenotype=='C3')*dat$sg
  form<-reformulate(c('phenotype','sg','int2','int3',adj),response='make30')
  fit<-glm(form,data=dat,family=poisson('log'),control=glm.control(maxit=100,epsilon=1e-9))
  b<-coef(fit);V<-vcovHC(fit,type='HC1');stopifnot(fit$converged,all(is.finite(b)),all(is.finite(V)),nobs(fit)==nrow(d))
  audit[[length(audit)+1]]<-data.frame(subgroup=sg,imputation=j,n=nobs(fit),parameters=length(b),converged=fit$converged,max_fitted=max(fitted(fit)),fitted_gt1=sum(fitted(fit)>1),max_hat=max(hatvalues(fit)))
  Q[[j]]<-b[c('int2','int3')];U[[j]]<-V[c('int2','int3'),c('int2','int3')]
  for(g in 0:1)for(cl in c('C2','C3')){L<-setNames(rep(0,length(b)),names(b));L[paste0('phenotype',cl)]<-1;L[paste0('int',substr(cl,2,2))]<-g
   key<-paste(g,cl);qn[[key]]<-c(qn[[key]],sum(L*b));un[[key]]<-c(un[[key]],as.numeric(t(L)%*%V%*%L))
   draws[[length(draws)+1]]<-data.frame(subgroup=sg,imputation=j,stratum=g,class=cl,Q=sum(L*b),U=as.numeric(t(L)%*%V%*%L))
  }
 }
 test<-mitml::testConstraints(qhat=Q,uhat=U,constraints=c('int2','int3'),method='D1',df.com=nrow(d)-length(b))
 inter[[length(inter)+1]]<-data.frame(subgroup=sg,F=test$test[1,1],df1=test$test[1,2],df2=test$test[1,3],p_interaction=test$test[1,4])
 joint[[sg]]<-list(Q=Q,U=U,test=test)
 for(g in 0:1)for(cl in c('C2','C3')){key<-paste(g,cl);p<-pool.scalar(qn[[key]],un[[key]],n=nrow(d),k=length(b));se<-sqrt(p$t);crit<-qt(.975,p$df)
  ans[[length(ans)+1]]<-data.frame(subgroup=sg,stratum=g,comparison=paste(cl,'vs C1'),RR=exp(p$qbar),low=exp(p$qbar-crit*se),high=exp(p$qbar+crit*se),p=2*pt(-abs(p$qbar/se),p$df),df=p$df,FMI=p$fmi,MCSE_SE=sqrt(p$b/30)/se)
 }
}
csv(do.call(rbind,ans),'RR');csv(do.call(rbind,inter),'interaction');csv(do.call(rbind,cnt),'observed_counts');csv(do.call(rbind,icount),'imputed_counts');csv(do.call(rbind,audit),'model_audit');csv(do.call(rbind,draws),'draws')
saveRDS(joint,file.path(P,'diagnostics',paste0(db,'_joint_tests.rds')))
capture.output(sessionInfo(),file=file.path(P,'diagnostics',paste0(db,'_session.txt')))
cat('COMPLETE',db,format(Sys.time()),'\n')

