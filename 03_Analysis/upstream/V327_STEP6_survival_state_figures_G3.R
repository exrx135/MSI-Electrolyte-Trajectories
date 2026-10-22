# V327: complete standalone death KM and bicarbonate-state alluvial figures.
# Run with Rscript --vanilla or source(..., encoding='UTF-8').
# Requires ggplot2 and survival. Frozen inputs and G3 classifications are never modified.
main <- function() {
  root <- 'Z:/Research/MSI-1'
  stopifnot(dir.exists(root))
  for (pkg in c('ggplot2','survival')) if (!requireNamespace(pkg,quietly=TRUE)) stop('Install required package: ',pkg)
  suppressPackageStartupMessages(library(ggplot2))
  out <- file.path(root,'results','V327_生存与状态转换_G3')
  dir.create(out,recursive=TRUE,showWarnings=FALSE)
  fn <- function(content,ext='csv') file.path(out,paste0('V327_',content,'_G3.',ext))
  csv <- function(x,name) write.csv(x,fn(name),row.names=FALSE,na='',fileEncoding='UTF-8')
  read <- function(p) read.csv(p,check.names=FALSE,stringsAsFactors=FALSE,na.strings=c('','NA','NULL'))
  raw_names <- c(MIMIC='MIMIC_FINAL_MAKE30_ELECTROLYTEQC.csv',eICU='eICU_FINAL_RAW_ELECTROLYTEQC.csv',SICdb='SICdb_FINAL_MAKE30_ELECTROLYTEQC.csv')
  class_names <- c(MIMIC='V318_MIMIC_formal_classification_G3.csv',eICU='V318_eICU_fixed_scoring_G3.csv',SICdb='V318_SICdb_fixed_scoring_G3.csv')
  ids <- c(MIMIC='stay_id',eICU='patientunitstayid',SICdb='case_id')
  expected_n <- c(MIMIC=9804L,eICU=5839L,SICdb=3831L)
  paths <- c(file.path(root,raw_names),file.path(root,'STEP5_冻结分型与预后外部验证','冻结分类',class_names))
  baseline <- read(file.path(root,'V325_文件总清单_GALL.csv'))
  before <- unname(tools::md5sum(paths))
  rel <- substring(paths,nchar(root)+2)
  stopifnot(!anyNA(before),all(before==baseline$MD5[match(rel,baseline$文件)]))
  cohorts <- list(); input_audit <- data.frame(file=rel,md5_before=before)
  for(db in names(raw_names)) {
    d <- read(file.path(root,raw_names[[db]])); cl <- read(file.path(root,'STEP5_冻结分型与预后外部验证','冻结分类',class_names[[db]]))
    stopifnot(nrow(d)==expected_n[[db]],nrow(cl)==nrow(d),!anyDuplicated(d[[ids[[db]]]]),!anyDuplicated(cl$unit_id))
    ix <- match(as.character(d[[ids[[db]]]]),as.character(cl$unit_id));stopifnot(!anyNA(ix),all(cl$assigned_class %in% 1:3))
    d$phenotype <- paste0('C',cl$assigned_class[ix]);d$unit_id <- as.character(d[[ids[[db]]]])
    cohorts[[db]] <- d
  }
  # Audit exact available columns. No death times imputed from censor/discharge.
  availability <- do.call(rbind,lapply(names(cohorts),function(db) {
    d<-cohorts[[db]]
    data.frame(dataset=db,classified_n=nrow(d),positive_death_n=sum(d$make30_death==1,na.rm=TRUE),
      explicit_death_time_available='death_time' %in% names(d),
      time_fields=paste(names(d)[grepl('death|censor|discharge|risk_start|landmark',names(d))],collapse=';'),
      survival_status=if(db=='MIMIC') 'Exploratory KM only; discharge censoring limitation' else 'NOT GENERATED: independent death time unavailable')
  })); csv(availability,'死亡字段可用性')
  # MIMIC uses exact stored window; endpoint-specific eligibility, not MAKE30 complete cases.
  d <- cohorts$MIMIC
  dt <- function(x) as.POSIXct(x,format='%d/%m/%Y %H:%M:%S',tz='UTC')
  start <- dt(d$landmark_time);end <- dt(d$make30_censor_time);death <- dt(d$death_time)
  stopifnot(!anyNA(start),!anyNA(end),max(abs(as.numeric(difftime(end,start,units='days'))-d$make30_followup_days_after_landmark))<1e-8)
  positive <- !is.na(d$make30_death)&d$make30_death==1
  stopifnot(!anyNA(death[positive]),all(death[positive]>start[positive]),all(death[positive]<=end[positive]))
  in_window <- !is.na(death)&death>start&death<=end
  stopifnot(all(in_window==positive),all(d$make30_risk_start==d$landmark_time))
  eligible <- d$make30_evaluable_window==1 & !is.na(d$make30_death) & end>start
  duration <- as.numeric(difftime(end,start,units='days'))
  duration[positive] <- as.numeric(difftime(death[positive],start[positive],units='days'))
  a <- data.frame(unit_id=d$unit_id[eligible],phenotype=factor(d$phenotype[eligible],levels=paste0('C',1:3)),
    time_days=duration[eligible],event=d$make30_death[eligible])
  stopifnot(all(a$time_days>0),all(a$time_days<30))
  csv(a,'MIMIC_死亡分析输入')
  csv(data.frame(classified_n=nrow(d),eligible_n=nrow(a),outside_window_n=sum(d$make30_evaluable_window!=1),unknown_death_n=sum(is.na(d$make30_death)),deaths=sum(a$event),
    unknown_MAKE30_retained_n=sum(is.na(d$make30[eligible])),time_origin='sepsis onset +24h',
    observation_end='stored censor: earlier hospital discharge / index ICU admission +30d',
    interpretation='Exploratory discharge-censored KM, NOT complete 30-day all-cause survival'), 'MIMIC_死亡分析流程')
  fit <- survival::survfit(survival::Surv(time_days,event)~phenotype,data=a,conf.type='log-log')
  sf <- summary(fit,censored=TRUE)
  curves <- data.frame(phenotype=sub('phenotype=','',as.character(sf$strata)),time=sf$time,survival=sf$surv,lower=sf$lower,upper=sf$upper,n_risk=sf$n.risk,n_event=sf$n.event,n_censor=sf$n.censor)
  curves$lower[is.na(curves$lower) & curves$survival==1] <- 1
  curves$upper[is.na(curves$upper) & curves$survival==1] <- 1
  csv(curves,'MIMIC_KM数值')
  zero <- data.frame(phenotype=paste0('C',1:3),time=0,survival=1,lower=1,upper=1,n_risk=as.integer(table(a$phenotype)),n_event=0,n_censor=0)
  curves <- rbind(zero,curves);curves<-curves[order(curves$phenotype,curves$time),]
  lr <- survival::survdiff(survival::Surv(time_days,event)~phenotype,data=a)
  pv <- pchisq(lr$chisq,length(lr$n)-1,lower.tail=FALSE)
  csv(data.frame(test='Unadjusted overall log-rank; exploratory',chisq=lr$chisq,df=2,p_value=pv),'MIMIC_logrank')
  ticks <- c(0,7,14,21,28)
  risk <- do.call(rbind,lapply(levels(a$phenotype),function(k) do.call(rbind,lapply(ticks,function(t) data.frame(phenotype=k,time=t,n_risk=sum(a$phenotype==k & a$time_days>=t))))))
  csv(risk,'MIMIC_风险人数')
  totals <- do.call(rbind,lapply(levels(a$phenotype),function(k) {x<-a[a$phenotype==k,];data.frame(phenotype=k,n=nrow(x),deaths=sum(x$event),observed_death_pct=100*mean(x$event))}));csv(totals,'MIMIC_死亡计数')
  cols <- c(C1='#3269A8',C2='#D46A45',C3='#40866B')
  theme_pub <- function(size=11) theme_minimal(base_size=size,base_family='sans')+theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),panel.grid.major.y=element_line(colour='#E5E9ED',linewidth=.3),plot.title=element_text(face='bold',colour='#152C40'),plot.subtitle=element_text(colour='#536473'),axis.title=element_text(colour='#283E4A'),legend.title=element_blank(),legend.position='top',plot.margin=margin(8,15,8,12))
  # Duplicate step boundaries so the CI ribbon is a step function, not interpolated.
  band <- do.call(rbind,lapply(split(curves,curves$phenotype),function(x) {
    if(nrow(x)<2)return(x)
    y<-x[rep(seq_len(nrow(x)-1),each=2),]; y$time<-as.vector(rbind(head(x$time,-1),tail(x$time,-1)));rbind(y,tail(x,1))
  }))
  km <- ggplot(curves,aes(time,survival,colour=phenotype,group=phenotype))+
    geom_ribbon(data=band,aes(ymin=lower,ymax=upper,fill=phenotype),colour=NA,alpha=.12,show.legend=FALSE)+
    geom_step(linewidth=1.05,direction='hv')+
    scale_colour_manual(values=cols,labels=setNames(sprintf('%s  |  n = %s; deaths = %s',totals$phenotype,format(totals$n,big.mark=','),totals$deaths),totals$phenotype))+
    scale_fill_manual(values=cols)+scale_x_continuous(breaks=ticks,limits=c(-.7,29.5),expand=c(0,0))+
    scale_y_continuous(limits=c(0,1),breaks=seq(0,1,.2),labels=function(x)paste0(round(x*100),'%'))+
    labs(title='Death after the 24-hour landmark',subtitle='MIMIC-IV  |  Frozen G3 phenotypes  |  Exploratory Kaplan-Meier estimates',x=NULL,y='Kaplan-Meier survival estimate')+
    annotate('text',x=1,y=.12,label=if(pv<.001)'Overall log-rank P < 0.001' else sprintf('Overall log-rank P = %.3f',pv),hjust=0,size=3.6,colour='#374A58')+theme_pub(12)
  rp<-ggplot(risk,aes(time,phenotype,label=format(n_risk,big.mark=','),colour=phenotype))+geom_text(size=3.8)+scale_colour_manual(values=cols)+
    scale_x_continuous(breaks=ticks,limits=c(-.7,29.5),expand=c(0,0))+scale_y_discrete(limits=rev(paste0('C',1:3)))+
    labs(x='Days since the 24-hour landmark',y=NULL,title='Number at risk')+theme_pub(11)+theme(legend.position='none',panel.grid=element_blank(),plot.title=element_text(size=10,face='bold'))
  kg <- ggplotGrob(km); rg <- ggplotGrob(rp); widths <- grid::unit.pmax(kg$widths,rg$widths); kg$widths <- widths; rg$widths <- widths
  draw_km<-function(){grid::grid.newpage();grid::pushViewport(grid::viewport(layout=grid::grid.layout(3,1,heights=grid::unit(c(.69,.22,.09),'null'))));grid::pushViewport(grid::viewport(layout.pos.row=1));grid::grid.draw(kg);grid::popViewport();grid::pushViewport(grid::viewport(layout.pos.row=2));grid::grid.draw(rg);grid::popViewport();grid::grid.text('Shading: pointwise 95% CI. Non-deaths are censored at the stored observation end.\nHospital discharge may be informative; these are not complete 30-day all-cause survival probabilities.',x=.055,y=.48,just='left',gp=grid::gpar(fontsize=9,col='#536473'),vp=grid::viewport(layout.pos.row=3));grid::popViewport()}
  save_draw<-function(name,fun,w,h){grDevices::cairo_pdf(fn(name,'pdf'),width=w,height=h);tryCatch(fun(),finally=dev.off());grDevices::png(fn(name,'png'),width=w,height=h,units='in',res=300,type='cairo');tryCatch(fun(),finally=dev.off())}
  save_draw('Figure_死亡KM_MIMIC_探索性',draw_km,10,8)
  # HCO3 state thresholds: adult serum bicarbonate reference 22-29 mmol/L.
  # Fixed before inspecting state-outcome associations; no outcome data enter this plot.
  states<-c('Low','Reference','High','Missing'); statecols<-c(Low='#CF6848',Reference='#429A89',High='#557DB5',Missing='#C9CFD5')
  state<-function(v) ifelse(is.na(v),'Missing',ifelse(v<22,'Low',ifelse(v<=29,'Reference','High')))
  hour<-c(3,9,15,21);nodes<-list();edges<-list();polygons<-list();nn<-ee<-pp<-0L
  for(db in names(cohorts))for(k in paste0('C',1:3)) {
    x<-cohorts[[db]][cohorts[[db]]$phenotype==k,]; n<-nrow(x)
    stopifnot(all(paste0('bicarbonate_t',1:4)%in%names(x)))
    st<-sapply(x[paste0('bicarbonate_t',1:4)],state)
    counts<-sapply(1:4,function(t)as.integer(table(factor(st[,t],levels=states))))
    bottoms<-apply(counts/n*100,2,function(v)c(0,head(cumsum(v),-1)))
    for(t in 1:4)for(i in 1:4){nn<-nn+1L;nodes[[nn]]<-data.frame(dataset=db,phenotype=k,class_n=n,time_hour=hour[t],time_index=t,state=states[i],n=counts[i,t],pct=counts[i,t]/n*100,ymin=bottoms[i,t],ymax=bottoms[i,t]+counts[i,t]/n*100)}
    for(t in 1:3) {
      mat<-table(factor(st[,t],levels=states),factor(st[,t+1],levels=states));stopifnot(all(rowSums(mat)==counts[,t]),all(colSums(mat)==counts[,t+1]),sum(mat)==n)
      src<-bottoms[,t];dst<-bottoms[,t+1]
      for(i in 1:4)for(j in 1:4) {
        count<-as.integer(mat[i,j]);ee<-ee+1L;edges[[ee]]<-data.frame(dataset=db,phenotype=k,from_hour=hour[t],to_hour=hour[t+1],from_state=states[i],to_state=states[j],n=count,pct=100*count/n)
        if(count==0)next
        h<-100*count/n; u<-seq(0,1,length.out=45);smooth<-u*u*(3-2*u);xs<-t+.12+(1-.24)*u
        lo<-src[i]+(dst[j]-src[i])*smooth;hi<-lo+h
        pp<-pp+1L;polygons[[pp]]<-data.frame(dataset=db,phenotype=k,state=states[i],group=pp,x=c(xs,rev(xs)),y=c(lo,rev(hi)))
        src[i]<-src[i]+h;dst[j]<-dst[j]+h
      }
    }
  }
  nd<-do.call(rbind,nodes);ed<-do.call(rbind,edges);poly<-do.call(rbind,polygons)
  csv(nd,'HCO3_状态占比');csv(ed,'HCO3_相邻状态流量')
  for(obj in c('nd','poly')) {z<-get(obj);z$dataset<-factor(z$dataset,levels=c('MIMIC','eICU','SICdb'),labels=c('MIMIC-IV','eICU-CRD','SICdb'));z$phenotype<-factor(z$phenotype,levels=paste0('C',1:3));assign(obj,z)}
  ann<-unique(nd[c('dataset','phenotype','class_n')]);ann$label<-paste0('n = ',format(ann$class_n,big.mark=','))
  p<-ggplot()+geom_polygon(data=poly,aes(x,y,group=group,fill=state),alpha=.34,colour=NA)+
    geom_rect(data=nd,aes(xmin=time_index-.12,xmax=time_index+.12,ymin=ymin,ymax=ymax,fill=state),colour='white',linewidth=.15)+
    geom_text(data=nd[nd$pct>=7,],aes(time_index,(ymin+ymax)/2,label=paste0(round(pct),'%')),size=2.8,colour='white',fontface='bold')+
    geom_text(data=ann,aes(x=1,y=107,label=label),hjust=0,size=3.3,colour='#3D5160')+
    facet_grid(dataset~phenotype,axes="all_x")+scale_fill_manual(values=statecols,breaks=states,labels=c('Low: <22','Reference: 22-29','High: >29','Missing'),drop=FALSE)+
    scale_x_continuous(breaks=1:4,labels=paste0(hour,' h'),limits=c(.7,4.3),expand=c(0,0))+
    scale_y_continuous(breaks=c(0,25,50,75,100),labels=function(x)paste0(x,'%'),limits=c(0,112),expand=c(0,0))+
    labs(title='Bicarbonate state transitions within fixed joint phenotypes',subtitle='Three independent cohorts  |  Observed 6-hour windows  |  No imputation',x='Hours after sepsis onset (window midpoint)',y='Patients within each fixed phenotype',
      caption='HCO3 in mmol/L. Ribbons connect the same patients between adjacent windows; colour denotes the source state.\nAll classified patients are retained. Missing measurements are a separate state. C1-C3 do not change over time.')+
    theme_pub(11)+theme(strip.text=element_text(face='bold',size=11,colour='#243F53'),strip.background=element_rect(fill='#EDF2F6',colour=NA),panel.spacing=grid::unit(1,'lines'),panel.grid=element_blank(),legend.position='top',plot.caption=element_text(hjust=0,size=9,colour='#536473',margin=margin(t=12)),axis.text.x=element_text(size=9))
  save_draw('Figure_HCO3状态转换_三库',function()print(p),12,10.5)
  csv(data.frame(state=states,definition=c('<22 mmol/L','22-29 mmol/L inclusive','>29 mmol/L','NA / not measured'),source='Mayo Clinic Laboratories HCO3 adult reference: https://www.mayocliniclabs.com/test-catalog/overview/876',role='Harmonized display bins; not database-specific reference intervals or acid-base diagnoses'),'HCO3_状态定义')
  # Conservation of every patient across windows and every interval is checked above.
  # Verify interval totals against the immutable classified population.
  for(db in names(cohorts)) for(k in paste0('C',1:3)) {
    n<-sum(cohorts[[db]]$phenotype==k)
    for(t in hour[-4]) stopifnot(sum(ed$n[ed$dataset==db & ed$phenotype==k & ed$from_hour==t])==n)
  }
  input_audit$md5_after<-unname(tools::md5sum(paths));input_audit$unchanged<-input_audit$md5_before==input_audit$md5_after
  stopifnot(all(input_audit$unchanged));csv(input_audit,'输入身份核验')
  writeLines(capture.output(sessionInfo()),fn('R运行环境','txt'),useBytes=TRUE)
  csv(data.frame(figure=c('MIMIC death KM','Three-cohort HCO3 state alluvial'),status=c('GENERATED_EXPLORATORY_NOT_MAIN_READY','GENERATED_DESCRIPTIVE'),
    limitation=c('Only MIMIC has explicit death time; competing live discharge cannot be separated from administrative censoring with current fields. External curves not generated.','Common display reference 22-29; does not establish acid-base disorder or treatment response.')),'完成状态')
  files<-list.files(out,full.names=TRUE);files<-files[!file.info(files)$isdir & basename(files)!=basename(fn('输出清单'))]
  csv(data.frame(filename=basename(files),bytes=file.info(files)$size,md5=unname(tools::md5sum(files))),'输出清单')
  print(totals);cat('COMPLETE: ',out,'\n',sep='')
}
main()
