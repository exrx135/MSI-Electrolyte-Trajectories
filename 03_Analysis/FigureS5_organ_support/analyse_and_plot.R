# Descriptive analysis of frozen G3 class, early organ support and MAKE30.
# Usage: Rscript --vanilla analyse_and_plot.R [project_root] [output_directory] [--plot-only]
# Reads existing cohort data only; does not refit classes, impute or query a database.
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)>=1) args[1] else 'Z:/Research/MSI-1'
out <- if(length(args)>=2) args[2] else file.path(root,'analysis_updates/FigureS5_organ_support_20261019/outputs')
plot_only <- length(args)>=3 && args[3]=='--plot-only'
dir.create(out,recursive=TRUE,showWarnings=FALSE)
src <- file.path(root,'\u5386\u53f2\u5f52\u6863/results_\u6574\u7406\u524d_20261019/results/table1/derived')
dbs <- c('MIMIC','eICU','SICdb'); names_db <- c(MIMIC='MIMIC-IV',eICU='eICU-CRD',SICdb='SICdb')
expected_n <- c(MIMIC=9758,eICU=5709,SICdb=3635)
expected_events <- c(MIMIC=2442,eICU=2060,SICdb=949)
expected_class <- list(MIMIC=c(4263,2833,2662),eICU=c(1756,2523,1430),SICdb=c(1064,908,1663))
support_labels <- c('Neither support','Vasopressors only','Mechanical ventilation only','Both supports')
classes <- paste0('C',1:3)
cols <- c(C1='#0D4483',C2='#DD1215',C3='#39A08F')
all_counts <- all_summary <- all_distribution <- list(); audits <- list()
wilson <- function(e,n){
  z<-qnorm(.975); p<-e/n; den<-1+z*z/n
  m<-(p+z*z/(2*n))/den; h<-z*sqrt(p*(1-p)/n+z*z/(4*n*n))/den
  c(lower=max(0,m-h),upper=min(1,m+h))
}
if(plot_only){
  counts<-read.csv(file.path(out,'Sankey_path_counts.csv'),check.names=FALSE)
  for(db in dbs){
    tab<-counts[counts$cohort==db,]
    stopifnot(sum(tab$n)==expected_n[db],sum(tab$n*tab$event)==expected_events[db])
  }
} else {
for(db in dbs){
  file <- file.path(src,paste0(db,'_Table1_analysis.csv'))
  d<-read.csv(file,check.names=FALSE,fileEncoding='UTF-8-BOM')
  vars<-c('cohort_unit_id','assigned_class','mechanical_ventilation_0_24h_analysis','vasopressor_0_24h_analysis','make30_analysis')
  stopifnot(all(vars %in% names(d)),!anyDuplicated(d$cohort_unit_id))
  initial_n<-nrow(d)
  d<-d[!is.na(d$make30_analysis) & d$make30_analysis %in% c(0,1),vars]
  names(d)<-c('id','class','mv','vaso','event')
  stopifnot(!anyNA(d),all(d$class %in% 1:3),all(d$mv %in% 0:1),all(d$vaso %in% 0:1))
  stopifnot(nrow(d)==expected_n[db],sum(d$event)==expected_events[db],all(as.numeric(table(factor(d$class,levels=1:3)))==expected_class[[db]]))
  d$support<-d$vaso+2*d$mv
  tab<-as.data.frame(table(class=factor(d$class,levels=1:3),support=factor(d$support,levels=0:3),event=factor(d$event,levels=0:1)))
  names(tab)[4]<-'n';tab$class<-paste0('C',tab$class);tab$support<-as.integer(as.character(tab$support));tab$event<-as.integer(as.character(tab$event));tab$cohort<-db
  tab<-tab[,c('cohort','class','support','event','n')]
  stopifnot(sum(tab$n)==nrow(d),sum(tab$n*tab$event)==sum(d$event))
  all_counts[[db]]<-tab
  records<-list(); idx<-0
  for(s in 0:3) for(cl in c(classes,'Overall')){
    z<-tab[tab$support==s & (tab$class==cl | cl=='Overall'),]
    n<-sum(z$n);e<-sum(z$n*z$event);ci<-wilson(e,n)
    class_n<-if(cl=='Overall') nrow(d) else sum(tab$n[tab$class==cl])
    idx<-idx+1;records[[idx]]<-data.frame(cohort=db,class=cl,support=s,support_label=support_labels[s+1],n=n,events=e,event_percent=100*e/n,lower_percent=100*ci[1],upper_percent=100*ci[2],class_n=class_n,within_class_support_percent=100*n/class_n,row.names=NULL)
  }
  all_summary[[db]]<-do.call(rbind,records)
  audits[[db]]<-data.frame(cohort=db,source=enc2utf8(file),source_sha256=digest::digest(file=file,algo='sha256'),source_n=initial_n,analysis_n=nrow(d),make30_events=sum(d$event),additional_support_exclusions=0)
}
counts<-do.call(rbind,all_counts);summary<-do.call(rbind,all_summary);audit<-do.call(rbind,audits)
write.csv(counts,file.path(out,'Sankey_path_counts.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(summary,file.path(out,'Support_and_MAKE30_summary.csv'),row.names=FALSE,fileEncoding='UTF-8')
write.csv(audit,file.path(out,'Input_and_cohort_audit.csv'),row.names=FALSE,fileEncoding='UTF-8')
contrasts<-list();j<-0
for(db in dbs) for(s in 0:3){
 z<-summary[summary$cohort==db & summary$support==s & summary$class!='Overall',];z<-z[match(classes,z$class),]
 for(k in 2:3){j<-j+1;contrasts[[j]]<-data.frame(cohort=db,support=s,support_label=support_labels[s+1],comparison=paste0(classes[k],' vs C1'),risk_difference_percentage_points=z$event_percent[k]-z$event_percent[1],crude_risk_ratio=z$event_percent[k]/z$event_percent[1])}
}
write.csv(do.call(rbind,contrasts),file.path(out,'Descriptive_within_support_contrasts.csv'),row.names=FALSE)
}

# Editable vector alluvial plot. Ribbon thickness is proportional to counts;
# each cohort is normalised separately to 100%. Colors retain original class.
W<-2700;H<-1110; svg<-c(sprintf('<svg xmlns="http://www.w3.org/2000/svg" width="%s" height="%s" viewBox="0 0 %s %s">',W,H,W,H),'<title>Trajectory classes, early organ support and MAKE30 across three cohorts</title>','<rect width="100%" height="100%" fill="white"/>')
add<-function(s)svg<<-c(svg,s)
esc<-function(s){s<-gsub('&','&amp;',s,fixed=TRUE);s<-gsub('<','&lt;',s,fixed=TRUE);gsub('>','&gt;',s,fixed=TRUE)}
text<-function(x,y,label,size=23,color='#20272D',bold=FALSE,anchor='middle')add(sprintf('<text x="%.3f" y="%.3f" font-family="Arial,Helvetica,sans-serif" font-size="%.2f" fill="%s" font-weight="%s" text-anchor="%s">%s</text>',x,y,size,color,if(bold)'bold'else'normal',anchor,esc(label)))
rect<-function(x,y,w,h,fill,stroke='none')add(sprintf('<rect x="%.3f" y="%.3f" width="%.3f" height="%.3f" fill="%s" stroke="%s" stroke-width="1.3"/>',x,y,w,h,fill,stroke))
ribbon<-function(x1,x2,y1,y2,h,color){dx<-(x2-x1)*.48;add(sprintf('<path d="M %.3f %.3f C %.3f %.3f %.3f %.3f %.3f %.3f L %.3f %.3f C %.3f %.3f %.3f %.3f %.3f %.3f Z" fill="%s" fill-opacity="0.34" stroke="white" stroke-width="0.65"/>',x1,y1,x1+dx,y1,x2-dx,y2,x2,y2,x2,y2+h,x2-dx,y2+h,x1+dx,y1+h,x1,y1+h,color))}
num<-function(n)format(n,big.mark=',',scientific=FALSE,trim=TRUE)
# Keep every column between y=235 and y=1051, with the same 720/N count scale.
# Distribute the same total blank space (96) across its available node gaps.
pos<-function(ns,N){scale<-720/N;gap<-96/(length(ns)-1);start<-235;setNames(start+c(0,head(cumsum(ns*scale+gap),-1)),names(ns))}
for(di in seq_along(dbs)){
 db<-dbs[di];N<-expected_n[db];tab<-counts[counts$cohort==db,];off<-(di-1)*900;scale<-720/N
 x<-c(off+38,off+347,off+708);width<-c(112,177,153)
 text(off+38,62,LETTERS[di],32,bold=TRUE,anchor='start');text(off+450,65,names_db[db],36,bold=TRUE)
 text(off+450,104,paste0('n = ',num(N)),24,color='#59636B')
 text(x[1]+width[1]/2,173,'Trajectory',23,bold=TRUE);text(x[1]+width[1]/2,202,'class',23,bold=TRUE)
 text(x[2]+width[2]/2,173,'Early organ',23,bold=TRUE);text(x[2]+width[2]/2,202,'support',23,bold=TRUE)
 text(x[3]+width[3]/2,173,'MAKE30',23,bold=TRUE);text(x[3]+width[3]/2,202,'outcome',23,bold=TRUE)
 nc<-setNames(vapply(classes,function(cl)sum(tab$n[tab$class==cl]),numeric(1)),classes)
 ns<-setNames(vapply(0:3,function(s)sum(tab$n[tab$support==s]),numeric(1)),0:3)
 ne<-setNames(vapply(0:1,function(e)sum(tab$n[tab$event==e]),numeric(1)),0:1)
 pc<-pos(nc,N);ps<-pos(ns,N);pe<-pos(ne,N);curc<-pc;curs<-ps
 for(cl in classes)for(s in 0:3){n<-sum(tab$n[tab$class==cl & tab$support==s]);h<-n*scale;ribbon(x[1]+width[1],x[2],curc[cl],curs[as.character(s)],h,cols[cl]);curc[cl]<-curc[cl]+h;curs[as.character(s)]<-curs[as.character(s)]+h}
 curs<-ps;cure<-pe
 # Same class/support ordering at both sides of each intermediate node.
 for(cl in classes)for(s in 0:3)for(e in 0:1){n<-tab$n[tab$class==cl & tab$support==s & tab$event==e];h<-n*scale;ribbon(x[2]+width[2],x[3],curs[as.character(s)],cure[as.character(e)],h,cols[cl]);curs[as.character(s)]<-curs[as.character(s)]+h;cure[as.character(e)]<-cure[as.character(e)]+h}
 for(cl in classes){h<-nc[cl]*scale;rect(x[1],pc[cl],width[1],h,cols[cl]);cy<-pc[cl]+h/2;text(x[1]+width[1]/2,cy-20,cl,29,'white',TRUE);text(x[1]+width[1]/2,cy+12,num(nc[cl]),22,'white');text(x[1]+width[1]/2,cy+42,sprintf('(%.1f%%)',100*nc[cl]/N),21,'white')}
 mid_labels<-c('Neither','Vaso only','MV only','Both')
 for(s in 0:3){key<-as.character(s);h<-ns[key]*scale;rect(x[2],ps[key],width[2],h,'#EEF1F3','#BDC6CD');cy<-ps[key]+h/2;fs<-if(h<45)17 else 22;text(x[2]+width[2]/2,cy-3,mid_labels[s+1],fs,bold=TRUE);text(x[2]+width[2]/2,cy+if(h<45)12 else 21,sprintf('%s (%.1f%%)',num(ns[key]),100*ns[key]/N),if(h<45)15 else 20)}
 for(e in 0:1){key<-as.character(e);h<-ne[key]*scale;rect(x[3],pe[key],width[3],h,if(e==0)'#E3EBE9'else'#607182');cy<-pe[key]+h/2;col<-if(e==0)'#253A35'else'white';text(x[3]+width[3]/2,cy-21,if(e==0)'No MAKE30'else'MAKE30',22,col,TRUE);text(x[3]+width[3]/2,cy+12,num(ne[key]),23,col);text(x[3]+width[3]/2,cy+44,sprintf('(%.1f%%)',100*ne[key]/N),21,col)}
}
svg<-c(svg,'</svg>')
writeLines(svg,file.path(out,'Figure S5. Trajectory classes, early organ support and MAKE30.svg'),useBytes=TRUE)
if(!plot_only){
  writeLines(capture.output(sessionInfo()),file.path(out,'R_session_info.txt'))
  cat('Completed. Frozen cohort sizes and outcome totals reconciled.\n')
  print(audit[,c('cohort','analysis_n','make30_events','additional_support_exclusions')],row.names=FALSE)
} else cat('Redrawn from frozen aggregate counts only; all three columns have aligned top and bottom edges.\n')
