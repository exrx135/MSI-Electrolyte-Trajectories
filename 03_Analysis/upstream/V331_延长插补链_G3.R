# Continue saved V331 chains; final diagnostics use the latter half of >=60 iterations.
args<-commandArgs(trailingOnly=TRUE);stopifnot(length(args)==1,args%in%c('MIMIC','eICU','SICdb'));db<-args[1]
out<-'Z:/Research/MSI-1/results/V331_MAR_MNAR敏感性_G3';dest<-file.path(out,paste0('V331_',db,'_MICE20_G3.rds'));imp<-readRDS(dest)
rh<-function(v){h<-floor(nrow(v)/2);v<-cbind(v[1:h,,drop=FALSE],v[(nrow(v)-h+1):nrow(v),,drop=FALSE]);W<-mean(apply(v,2,var));B<-nrow(v)*var(colMeans(v));if(W==0)return(if(B==0)1 else Inf);sqrt(((nrow(v)-1)/nrow(v)*W+B/nrow(v))/W)}
diagfun<-function(x){it<-x$iteration;do.call(rbind,lapply(names(x$data)[colSums(is.na(x$data))>0],function(v)data.frame(variable=v,mean_rhat=rh(x$chainMean[v,(floor(it/2)+1):it,]),sd_rhat=rh(sqrt(x$chainVar[v,(floor(it/2)+1):it,])),recent10_rhat=rh(x$chainMean[v,(it-9):it,]),iterations=it)))}
repeat {
 dg<-diagfun(imp)
 if(imp$iteration>=60 && max(dg$mean_rhat,dg$sd_rhat)<=1.1)break
 if(imp$iteration>=160)break
 extra<-if(imp$iteration<100)100-imp$iteration else 30
 cat('CONTINUE',db,imp$iteration,'+',extra,'\n')
 imp<-mice::mice.mids(imp,maxit=extra,printFlag=TRUE)
 tmp<-paste0(dest,'.tmp');saveRDS(imp,tmp);stopifnot(file.copy(tmp,dest,overwrite=TRUE));unlink(tmp)
}
write.csv(dg,file.path(out,paste0('V331_',db,'_延长链诊断_G3.csv')),row.names=FALSE)
cat('DONE ',db,' iterations=',imp$iteration,' mean_Rhat=',max(dg$mean_rhat),' sd_Rhat=',max(dg$sd_rhat),'\n')
if(max(dg$mean_rhat,dg$sd_rhat)>1.1)warning('Diagnostic criterion not met: report results as provisional')
