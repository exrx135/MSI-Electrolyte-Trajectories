# V331 optional independent preparation of a cohort MICE checkpoint.
# Same input, seed, methods, and diagnostic extension as the main sensitivity script.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1,args%in%c('eICU','SICdb'));db<-args[1]
root<-'Z:/Research/MSI-1';out<-file.path(root,'results','V331_MAR_MNAR敏感性_G3')
raws<-c(eICU='eICU_FINAL_RAW_ELECTROLYTEQC.csv',SICdb='SICdb_FINAL_MAKE30_ELECTROLYTEQC.csv')
rawcols<-c(paste0('msi_m',1:12),unlist(lapply(c('sodium','potassium','chloride','bicarbonate'),function(v)paste0(v,'_t',1:4))))
CORE<-c('age','sex','sepsis_sofa_score','baseline_egfr_mdrd','ckd_or_renal_disease','diabetes','mechanical_ventilation_0_24h','vasopressor_0_24h')
d<-read.csv(file.path(root,raws[[db]]),check.names=FALSE,stringsAsFactors=FALSE,na.strings=c('','NA','NULL'));input<-d[c(rawcols,CORE)];input$sex<-factor(input$sex);stopifnot(!anyNA(input[CORE]))
method<-setNames(rep('',ncol(input)),names(input));method[rawcols]<-'pmm';pred<-matrix(0,ncol(input),ncol(input),dimnames=list(names(input),names(input)));pred[rawcols,]<-1;diag(pred)<-0
imp<-mice::mice(input,m=20,maxit=20,method=method,predictorMatrix=pred,seed=20261012+match(db,c('MIMIC','eICU','SICdb')),donors=5,printFlag=TRUE)
rhat<-function(arr){v<-arr[max(1,nrow(arr)-9):nrow(arr),,drop=FALSE];half<-floor(nrow(v)/2);v<-cbind(v[1:half,,drop=FALSE],v[(nrow(v)-half+1):nrow(v),,drop=FALSE]);W<-mean(apply(v,2,var));B<-nrow(v)*var(colMeans(v));if(W==0)return(if(B==0)1 else Inf);sqrt(((nrow(v)-1)/nrow(v)*W+B/nrow(v))/W)}
rh<-sapply(rawcols,function(v)if(anyNA(input[[v]]))rhat(imp$chainMean[v,,,drop=TRUE]) else NA_real_)
if(any(rh>1.1,na.rm=TRUE))imp<-mice::mice.mids(imp,maxit=20,printFlag=TRUE)
dest<-file.path(out,paste0('V331_',db,'_MICE20_G3.rds'));tmp<-paste0(dest,'.tmp');saveRDS(imp,tmp);stopifnot(file.rename(tmp,dest));cat('PREPARED ',db,' iterations=',imp$iteration,'\n')
