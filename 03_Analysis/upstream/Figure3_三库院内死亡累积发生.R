# Standalone reproducible analysis. R >=4.5; survival, ggplot2, patchwork.
# Uses frozen G3 classes and user-adjudicated NEW survival endpoints.
# NEVER overwrites raw data, MAKE30, classification or prior model outputs.
# Live-discharge endpoint uses the retained frozen rule; see delivery README.
root <- 'Z:/Research/MSI-1'
out <- file.path(root,'results','生存分析_三库图与代码')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
library(survival);library(ggplot2);library(patchwork)
read <- function(p) read.csv(p,stringsAsFactors=FALSE,na.strings=c('','NA','NULL'))
csv <- function(x,n) write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE,na='',fileEncoding='UTF-8')
dt <- function(x) as.numeric(as.POSIXct(x,format='%d/%m/%Y %H:%M:%S',tz='UTC'))
src <- file.path(root,'results','生存分析_补充数据导出')
ids <- c(MIMIC='stay_id',eICU='patientunitstayid',SICdb='case_id')
rawfiles <- c(MIMIC='MIMIC_FINAL_MAKE30_ELECTROLYTEQC.csv',eICU='eICU_FINAL_RAW_ELECTROLYTEQC.csv',SICdb='SICdb_FINAL_MAKE30_ELECTROLYTEQC.csv')
clsfiles <- c(MIMIC='V318_MIMIC_formal_classification_G3.csv',eICU='V318_eICU_fixed_scoring_G3.csv',SICdb='V318_SICdb_fixed_scoring_G3.csv')
expected <- c(MIMIC=9804,eICU=5839,SICdb=3831)
flows<-list();summaries<-list();plots<-list();validations<-list();hashes<-list()
for(db in names(ids)){
 sf<-list.files(src,pattern=paste0('^0[123]_',db,'_.*csv$'),full.names=TRUE)
 rp<-file.path(root,rawfiles[[db]]);cp<-file.path(root,'STEP5_冻结分型与预后外部验证','冻结分类',clsfiles[[db]])
 stopifnot(length(sf)==1)
 hashes[[db]]<-data.frame(path=c(sf,rp,cp),md5=unname(tools::md5sum(c(sf,rp,cp))))
 s<-read(sf);r<-read(rp);cl<-read(cp);id<-ids[[db]]
 stopifnot(!anyDuplicated(s[[id]]),!anyDuplicated(r[[id]]),!anyDuplicated(cl$unit_id),nrow(r)==expected[[db]])
 ix<-match(r[[id]],s[[id]]);cx<-match(r[[id]],cl$unit_id);stopifnot(!anyNA(ix),!anyNA(cx))
 extras<-s[!s[[id]]%in%r[[id]],,drop=FALSE];csv(extras,paste0(db,'_非冻结队列记录'))
 d<-s[ix,];d$phenotype<-paste0('C',cl$assigned_class[cx]);stopifnot(all(d$phenotype%in%c('C1','C2','C3')),all(d$make30_death==r$make30_death))
 rule<-rep('original_fields',nrow(d))
 if(db=='MIMIC'){
  start<-dt(d$landmark_time);adm<-dt(d$administrative_end_time);dis<-dt(d$admission_discharge_time);dead<-dt(d$frozen_death_time);status<-ifelse(d$hospital_expire_flag==1,'Expired',ifelse(d$hospital_expire_flag==0,'Alive',NA));unit<-86400
  stopifnot(all(dt(d$sepsis_time)+86400==start),all(dt(d$icu_intime)+30*86400==adm),all(dt(d$make30_risk_start)==start),all(pmin(dis,adm)==dt(d$make30_censor_time)))
 }else if(db=='eICU'){
  start<-d$landmark_offset;adm<-d$administrative_end_icu_offset_min;dis<-d$hospitaldischargeoffset;dead<-d$death_icu_offset_min_from_status;status<-d$hospitaldischargestatus;unit<-1440
  stopifnot(all(d$sepsis_offset+1440==start),all(adm==43200),all(d$make30_risk_start_hosp_offset+d$hospitaladmitoffset==start),all(pmin(dis,adm)==d$make30_censor_hosp_offset+d$hospitaladmitoffset),all(start*60==r$landmark_offset_sec),all(pmin(dis,adm)*60==r$make30_censor_offset_sec))
  ep<-file.path(src,'05_eICU_44例同次住院核查.csv');ex<-read(ep)
  hashes[[paste0(db,'_extra')]]<-data.frame(path=ep,md5=unname(tools::md5sum(ep)))
  for(ii in which(is.na(status)|!status%in%c('Alive','Expired'))){
   z<-ex[ex$index_patientunitstayid==d[[id]][ii],]
   zd<-z[z$unitdischargestatus=='Expired'&z$unitdischargelocation=='Death',]
   za<-z[z$unitdischargestatus=='Alive'&z$unitdischargelocation%in%c('Home','Other Hospital'),]
   if(nrow(zd)) {status[ii]<-'Expired';dead[ii]<-min(zd$unit_discharge_index_icu_min);rule[ii]<-'ICU_Expired_and_Death_destination'}
   else if(nrow(za)) {status[ii]<-'Alive';rule[ii]<-'Alive_with_Home_or_OtherHospital_destination'}
  }
  stopifnot(sum(rule=='ICU_Expired_and_Death_destination')==5,sum(rule=='Alive_with_Home_or_OtherHospital_destination')==4)
 }else{
  sp<-file.path(src,'04_SICdb_补充住院出院状态.csv');sx<-read(sp)
  stopifnot(!anyDuplicated(sx$case_id));si<-match(d$case_id,sx$case_id);stopifnot(!anyNA(si));sx<-sx[si,]
  hashes[[paste0(db,'_extra')]]<-data.frame(path=sp,md5=unname(tools::md5sum(sp)))
  start<-d$landmark_offset_sec;adm<-d$administrative_end_offset_sec;dis<-d$hospital_discharge_offset_sec;dead<-d$offset_of_death_sec;unit<-86400
  status<-ifelse(sx$hospital_discharge_type_ref==3130,'Expired',ifelse(sx$hospital_discharge_type_ref%in%c(3129,3131,3132),'Alive',NA))
  # Individual adjudications are private data and are not embedded in public code.
  adj_path <- Sys.getenv('MSI_SICDB_ADJUDICATIONS', unset=file.path(root,'private_inputs','sicdb_adjudications.csv'))
  if(!file.exists(adj_path)) stop('Required private adjudication file is absent: ',adj_path)
  adj <- read(adj_path)
  stopifnot(all(c('case_id','death_override','discharge_override','ignore_negative_discharge_day') %in% names(adj)), !anyDuplicated(adj$case_id))
  death_override <- adj$case_id[adj$death_override==1]
  discharge_override <- adj$case_id[adj$discharge_override==1]
  negative_discharge_override <- adj$case_id[adj$ignore_negative_discharge_day==1]
  stopifnot(length(death_override)==12, length(discharge_override)==6, length(negative_discharge_override)==1)
  stopifnot(all(death_override%in%d$case_id),all(discharge_override%in%d$case_id))
  status[d$case_id%in%death_override]<-'Expired'
  rule[]<-'hospital_discharge_status; frozen_discharge_endpoint_retained'
  rule[d$case_id%in%death_override]<-'user_adjudicated_12_conflicting_records_as_death'
  rule[d$case_id%in%discharge_override]<-'user_adjudicated_6_postdischarge_deaths_as_discharge'
  rule[d$case_id%in%negative_discharge_override]<-'user_adjudicated_death; negative_discharge_day_ignored'
  # Explicit hospital death uses its timestamp, NOT coarse day-derived discharge cutoff.
  # Live discharge retains frozen endpoint per instruction not to redesign prior rules.
  stopifnot(all(start==d$make30_risk_start_offset_sec),all(adm==d$icu_offset_sec+2592000))
 }
 end<-pmin(dis,adm)
 if(db=='SICdb'){isdeath<-!is.na(status)&status=='Expired'&!is.na(dead);end[isdeath]<-pmin(dead[isdeath],adm[isdeath])}
 window<-!is.na(start)&!is.na(end)&end>start
 unknown<-is.na(status)|!status%in%c('Alive','Expired')
 inconsistent<-!unknown&status=='Expired'&(is.na(dead)|dead<=start|if(db=='SICdb') FALSE else dead>dis)
 inconsistent[is.na(inconsistent)]<-TRUE
 event<-rep(0L,nrow(d));event[!unknown&status=='Expired'&!is.na(dead)&dead>start&dead<=end]<-1L
 event[!unknown&status=='Alive'&!is.na(dis)&dis<=adm&dis>start]<-2L
 delta<-as.integer(event==1)-d$make30_death
 time<-(end-start)/unit;time[event==1]<-(dead[event==1]-start[event==1])/unit
 eligible<-window&!unknown&!inconsistent
 reason<-ifelse(!window,'nonpositive_or_missing_window',ifelse(unknown,'unknown_hospital_status',ifelse(inconsistent,'inconsistent_death_time','included')))
 if(db=='SICdb'){
  reason[unknown]<-'unknown_hospital_status'
  reason[!unknown&status=='Expired'&is.na(dead)]<-'death_time_missing'
  reason[!unknown&status=='Expired'&!is.na(dead)&dead<=start]<-'death_at_or_before_landmark'
  eligible<-reason=='included'
 }
 audit<-data.frame(unit_id=d[[id]],phenotype=d$phenotype,time_days=time,event=event,eligible=eligible,reason=reason,rule=rule,frozen_make30_death=d$make30_death,new_minus_frozen_death_flag=delta)
 csv(audit,paste0(db,'_逐患者纳排审计'))
 csv(cbind(d[!eligible,,drop=FALSE],exclusion_reason=audit$reason[!eligible]),paste0(db,'_排除记录原始字段'))
 a<-audit[eligible,];stopifnot(all(a$time_days>0),all(a$time_days<=30));csv(a,paste0(db,'_竞争风险分析数据'))
 flows[[db]]<-data.frame(dataset=db,exported=nrow(s),classified=nrow(d),excluded_window=sum(reason=='nonpositive_or_missing_window'),excluded_unknown_status=sum(reason=='unknown_hospital_status'),excluded_inconsistent=sum(reason=='inconsistent_death_time'),excluded_death_before_landmark=sum(reason=='death_at_or_before_landmark'),excluded_missing_death_time=sum(reason=='death_time_missing'),analyzed=nrow(a),deaths=sum(a$event==1),live_discharges=sum(a$event==2),administrative_censors=sum(a$event==0),MAKE30_unknown_retained=sum(is.na(r$make30[eligible])),status='COMPUTED under documented adjudication and retained discharge-time rules')
 curve<-list();risk<-list();summaryrows<-list();errors<-c()
 for(k in c('C1','C2','C3')){
  b<-a[a$phenotype==k,];b$state<-factor(b$event,levels=0:2,labels=c('censor','death','discharge'))
  fit<-survfit(Surv(time_days,state)~1,data=b,conf.type='log-log');j<-match('death',fit$states)
  cv<-data.frame(phenotype=k,time=c(0,fit$time),cif=c(0,fit$pstate[,j]),lower=c(0,fit$lower[,j]),upper=c(0,fit$upper[,j]));cv$lower[is.na(cv$lower)&cv$cif==0]<-0;cv$upper[is.na(cv$upper)&cv$cif==0]<-0
  # Independent product-integral recurrence, including competing events and ties.
  ss<-1;ff<-0;manual<-numeric(length(fit$time))
  for(i in seq_along(fit$time)){tt<-fit$time[i];nr<-sum(b$time_days>=tt);nd<-sum(b$time_days==tt&b$event==1);nl<-sum(b$time_days==tt&b$event==2);ff<-ff+ss*nd/nr;ss<-ss*(1-(nd+nl)/nr);manual[i]<-ff}
  errors<-c(errors,max(abs(manual-fit$pstate[,j])));stopifnot(tail(errors,1)<1e-10,all(diff(cv$cif)>=-1e-12))
  curve[[k]]<-cv;risk[[k]]<-data.frame(phenotype=k,time=c(0,7,14,21,28),n_risk=sapply(c(0,7,14,21,28),function(t)sum(b$time_days>=t)))
  summaryrows[[k]]<-data.frame(dataset=db,phenotype=k,n=nrow(b),deaths=sum(b$event==1),live_discharges=sum(b$event==2),administrative_censors=sum(b$event==0),observed_death_pct=100*mean(b$event==1))
 }
 cv<-do.call(rbind,curve);ri<-do.call(rbind,risk);csv(cv,paste0(db,'_AJ曲线及95CI'));csv(ri,paste0(db,'_风险人数'));summaries[[db]]<-do.call(rbind,summaryrows)
 validations[[db]]<-data.frame(dataset=db,check=c('Independent AJ recurrence max error','Adjudicated death flags different from frozen','Missing matched IDs','Export death fallback count'),value=c(max(errors),sum(delta!=0),0,if(db=='MIMIC')sum(d$frozen_death_used_discharge_fallback) else NA))
 pal<-c(C1='#31566E',C2='#C56E52',C3='#438D80')
 # Duplicate x coordinates to render stepwise pointwise CI bands without smoothing.
 band<-do.call(rbind,lapply(split(cv,cv$phenotype),function(z){i<-rep(seq_len(nrow(z)),each=2);zz<-z[i,];zz$time<-as.vector(rbind(z$time,c(z$time[-1],tail(z$time,1))));zz}))
 p<-ggplot(cv,aes(time,cif*100,colour=phenotype))+geom_ribbon(data=band,aes(ymin=lower*100,ymax=upper*100,fill=phenotype),alpha=.13,colour=NA)+geom_step(linewidth=.85)+scale_color_manual(values=pal)+scale_fill_manual(values=pal,guide='none')+scale_x_continuous(breaks=c(0,7,14,21,28),limits=c(0,30))+scale_y_continuous(limits=c(0,45),breaks=seq(0,40,10))+labs(title=c(MIMIC='A  MIMIC-IV',eICU='B  eICU-CRD',SICdb='C  SICdb')[[db]],x=NULL,y=if(db=='MIMIC') 'Cumulative incidence of in-hospital death (%)' else NULL,colour='Phenotype')+theme_classic(base_size=11)+theme(legend.position='top',plot.title=element_text(face='bold'))
 pr<-ggplot(ri,aes(time,phenotype,label=n_risk,colour=phenotype))+geom_text(size=3.2)+scale_color_manual(values=pal)+scale_x_continuous(breaks=c(0,7,14,21,28),limits=c(0,30))+scale_y_discrete(limits=c('C3','C2','C1'))+labs(x='Days after the 24-h landmark',y=NULL,title='Number at risk')+theme_classic(base_size=10)+theme(legend.position='none',axis.line.y=element_blank(),axis.ticks.y=element_blank(),plot.title=element_text(size=10))
 plots[[db]]<-p/pr+plot_layout(heights=c(4,1.2))
 ggsave(file.path(out,paste0(db,'_院内死亡累积发生.png')),plots[[db]],width=7,height=6,dpi=300)
}
csv(do.call(rbind,flows),'三库分析流程');csv(do.call(rbind,summaries),'类别事件汇总');csv(do.call(rbind,validations),'数值核验');csv(do.call(rbind,hashes),'输入文件MD5')
combined<-wrap_plots(plots,ncol=3)
ggsave(file.path(out,'Figure3_三库院内死亡累积发生.png'),combined,width=15,height=6,dpi=400,bg='white')
ggsave(file.path(out,'Figure3_三库院内死亡累积发生.svg'),combined,width=15,height=6,device=grDevices::svg,bg='white')
hash_table<-do.call(rbind,hashes);stopifnot(all(unname(tools::md5sum(hash_table$path))==hash_table$md5))
writeLines(capture.output(sessionInfo()),file.path(out,'R_sessionInfo.txt'))
print(do.call(rbind,flows));print(do.call(rbind,summaries));print(do.call(rbind,validations))

