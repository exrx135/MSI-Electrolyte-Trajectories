p<-'Z:/Research/MSI-1/analysis_updates/S4_Model3_20261019'
files<-list.files(p,pattern='\\.rds$',recursive=TRUE,full.names=TRUE);stopifnot(length(files)==540)
logs<-list();diagnostics<-list()
for(f in files){x<-readRDS(f);if(!is.null(x$loggedEvents)){z<-unique(x$loggedEvents[,c('dep','meth','out')]);z$file<-basename(dirname(f));logs[[length(logs)+1]]<-z};stopifnot(all(is.finite(x$fit$Q)),all(x$fit$U>0))}
z<-unique(do.call(rbind,logs));print(z);write.csv(z,file.path(p,'Imputation_logged_events_review.csv'),row.names=FALSE)
cat('540 completed conditional imputations and fits verified.\n')
