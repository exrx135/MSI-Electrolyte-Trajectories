# V330: publication figures, fixed G3 S01 labels; descriptive quadratic fits with clustered CI; no GBMTM refit.
main <- function() {
root <- 'Z:/Research/MSI-1'
suppressPackageStartupMessages(library(ggplot2))
out <- file.path(root,'results','V330_发表图终版_G3');dir.create(out,recursive=TRUE,showWarnings=FALSE)
fn <- function(n,e='csv')file.path(out,paste0('V330_',n,'_G3.',e))
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
cn<-do.call(rbind,counts);rownames(cn)<-NULL;cn$pct<-100*cn$n/rep(sapply(cohorts,nrow),each=3);csv(cn,'类别人数')
marks<-c('MSI','Na','K','Cl','HCO3'); labs<-c(MSI='MSI (bpm/mmHg)',Na='Sodium (mmol/L)',K='Potassium (mmol/L)',Cl='Chloride (mmol/L)',HCO3='Bicarbonate (mmol/L)')
cols<-c(C1='#31566E',C2='#C56E52',C3='#438D80')
proposed<-c(C1='Relatively balanced',C2='Higher MSI / low bicarbonate',C3='Lower chloride / higher bicarbonate')
groups<-list(MSI=paste0('msi_m',1:12),Na=paste0('sodium_t',1:4),K=paste0('potassium_t',1:4),Cl=paste0('chloride_t',1:4),HCO3=paste0('bicarbonate_t',1:4))
hours<-list(MSI=seq(1,23,2),Na=c(3,9,15,21),K=c(3,9,15,21),Cl=c(3,9,15,21),HCO3=c(3,9,15,21))
# MSI bins use observed first-window MIMIC values, one measurement per available patient.
q<-as.numeric(quantile(cohorts$MIMIC$msi_m1,c(1/3,2/3),na.rm=TRUE,type=7))
bounds<-list(MSI=q,Na=c(135,145),K=c(3.6,5.2),Cl=c(98,107),HCO3=c(22,29))
limits_table<-do.call(rbind,lapply(marks,function(m)data.frame(indicator=m,lower=bounds[[m]][1],upper=bounds[[m]][2],rule=if(m=='MSI')'Lower <=q33; middle >q33 and <=q67; upper >q67' else 'Low <lower; reference inclusive; high >upper',source=if(m=='MSI')'MIMIC observed msi_m1; quantile type 7; frozen for all cohorts' else 'https://www.mayocliniclabs.com/test-catalog/overview/113632; HCO3 /876')));csv(limits_table,'状态界值')
theme_pub<-function()theme_minimal(base_size=11)+theme(panel.grid.minor=element_blank(),panel.grid.major.x=element_blank(),panel.grid.major.y=element_line(colour='#E8EDF1',linewidth=.3),strip.background=element_rect(fill='#EDF2F6',colour=NA),strip.text=element_text(face='bold',colour='#243F53'),plot.title=element_text(face='bold',size=17,colour='#173449'),plot.subtitle=element_text(size=10,colour='#526675'),legend.position='top',legend.title=element_blank(),plot.caption=element_text(hjust=0,size=9,colour='#526675',margin=margin(t=12)),panel.spacing=grid::unit(.8,'lines'),plot.margin=margin(12,14,12,14))
save<-function(p,n,w=12,h=11){ggsave(fn(n,'pdf'),p,width=w,height=h,device=grDevices::cairo_pdf);ggsave(fn(n,'png'),p,width=w,height=h,dpi=220,bg='white')}


# Descriptive within-class observed-data polynomial, not a new GBMTM or frozen model uncertainty.
curves<-list();fitstats<-list();allnodes<-list();alledges<-list()
states<-c('Low','Reference','High','Missing')
for(db in names(cohorts))for(k in paste0('C',1:3))for(m in marks){
 x<-cohorts[[db]][cohorts[[db]]$phenotype==k,];n<-nrow(x);mat<-as.matrix(x[groups[[m]]]);hr<-hours[[m]];nt<-length(hr)
 patient_key<-c(MIMIC='subject_id',eICU='uniquepid',SICdb='patient_id')[[db]];stopifnot(patient_key%in%names(x),!anyNA(x[[patient_key]]));long<-data.frame(patient=rep(x[[patient_key]],each=nt),t=rep((hr-12)/12,n),y=as.vector(t(mat)));long<-long[!is.na(long$y),]
 fit<-lm(y~t+I(t^2),data=long,x=TRUE,y=TRUE);X<-fit$x;bread<-solve(crossprod(X));scores<-rowsum(X*residuals(fit),long$patient,reorder=FALSE);nc<-nrow(scores);N<-nrow(X);P<-ncol(X)
 V<-bread%*%crossprod(scores)%*%bread*nc/(nc-1)*(N-1)/(N-P)
 # Independently verify the clustered covariance against sandwich::vcovCL.
 stopifnot(max(abs(V-sandwich::vcovCL(fit,cluster=long$patient,type='HC1',cadjust=TRUE)))<1e-8)
 tt<-seq(min(hr),max(hr),length.out=201);u<-(tt-12)/12;D<-cbind(1,u,u^2);mu<-as.numeric(D%*%coef(fit));se<-sqrt(rowSums((D%*%V)*D));crit<-qt(.975,nc-1)
 curves[[length(curves)+1]]<-data.frame(dataset=db,phenotype=k,indicator=m,hour=tt,mean=mu,lower=mu-crit*se,upper=mu+crit*se)
 fitstats[[length(fitstats)+1]]<-data.frame(dataset=db,phenotype=k,indicator=m,class_n=n,observed_patients=nc,observations=N,b0=coef(fit)[1],b1=coef(fit)[2],b2=coef(fit)[3],CI='Pointwise 95% cluster-robust t interval; patient cluster; fixed classes')
 b<-bounds[[m]];statefun<-function(v)ifelse(is.na(v),'Missing',ifelse(if(m=='MSI')v<=b[1] else v<b[1],'Low',ifelse(v<=b[2],'Reference','High')))
 st<-apply(mat,2,statefun);cnt<-sapply(1:nt,function(t)as.integer(table(factor(st[,t],levels=states))));bot<-apply(cnt/n*100,2,function(v)c(0,head(cumsum(v),-1)))
 for(t in 1:nt)for(i in 1:4)allnodes[[length(allnodes)+1]]<-data.frame(indicator=m,dataset=db,phenotype=k,class_n=n,time_hour=hr[t],state=states[i],n=cnt[i,t],pct=100*cnt[i,t]/n,ymin=bot[i,t],ymax=bot[i,t]+100*cnt[i,t]/n)
 for(t in seq_len(nt-1)){tr<-table(factor(st[,t],levels=states),factor(st[,t+1],levels=states));stopifnot(all(rowSums(tr)==cnt[,t]),all(colSums(tr)==cnt[,t+1]));for(i in 1:4)for(j in 1:4)alledges[[length(alledges)+1]]<-data.frame(indicator=m,dataset=db,phenotype=k,from_hour=hr[t],to_hour=hr[t+1],from_state=states[i],to_state=states[j],n=as.integer(tr[i,j]),pct=100*as.integer(tr[i,j])/n)}
}
cv<-do.call(rbind,curves);csv(cv,'三库平滑趋势及95CI');csv(do.call(rbind,fitstats),'描述性平滑模型信息')
cv$indicator<-factor(cv$indicator,levels=marks,labels=labs);cv$dataset<-factor(cv$dataset,levels=names(raws))
hh<-setNames(sapply(names(raws),function(db){z<-cn[cn$dataset==db,];paste0(db,'\n',paste0(z$phenotype,': ',format(z$n,big.mark=','),collapse=' | '))}),names(raws))
p<-ggplot(cv,aes(hour,mean,colour=phenotype,fill=phenotype))+geom_ribbon(aes(ymin=lower,ymax=upper),alpha=.16,colour=NA)+geom_line(linewidth=1.05)+facet_grid(indicator~dataset,scales='free_y',labeller=labeller(dataset=hh))+scale_colour_manual(values=cols)+scale_fill_manual(values=cols)+scale_x_continuous(breaks=c(0,6,12,18,24),limits=c(0,24))+labs(title='MSI-electrolyte trajectories across three cohorts',subtitle='Frozen G3 phenotype assignments | Observed-data quadratic trends with pointwise 95% confidence intervals',x='Hours after sepsis onset',y=NULL,caption='Shading: patient-cluster robust 95% confidence intervals for descriptive mean trends, conditional on fixed class assignments.\nThese are not refitted GBMTM trajectories or prediction intervals. Available observations only; no missing-value imputation.\nMSI: 12 two-hour windows; electrolytes: four six-hour windows. No point markers; each variable shares one y-axis across cohorts.')+theme_pub()+theme(strip.text.x=element_text(size=10),legend.text=element_text(size=11),panel.grid.major.y=element_line(colour='#EEF1F3',linewidth=.25))
save(p,'Main_三库五变量平滑轨迹95CI',13,12)
nd<-do.call(rbind,allnodes);ed<-do.call(rbind,alledges);csv(nd,'全变量状态占比');csv(ed,'全变量相邻状态流量')
states<-c('Low','Reference','High','Missing');sc<-c(Low='#CF6848',Reference='#429A89',High='#557DB5',Missing='#C9CFD5')
for(db in names(raws)){
 polys<-list();pp<-0
 for(m in marks)for(k in paste0('C',1:3)){
  hr<-hours[[m]];hw<-if(m=='MSI').6 else .9
  for(t in head(hr,-1)){
   z<-ed[ed$dataset==db & ed$indicator==m & ed$phenotype==k & ed$from_hour==t,];t2<-unique(z$to_hour);stopifnot(length(t2)==1)
   a<-nd[nd$dataset==db & nd$indicator==m & nd$phenotype==k & nd$time_hour==t,];b<-nd[nd$dataset==db & nd$indicator==m & nd$phenotype==k & nd$time_hour==t2,]
   a<-a[match(states,a$state),];b<-b[match(states,b$state),];src<-a$ymin;dst<-b$ymin;stopifnot(sum(z$n)==a$class_n[1])
   for(i in 1:4)for(j in 1:4){e<-z[z$from_state==states[i]&z$to_state==states[j],];stopifnot(nrow(e)==1);h<-e$pct;if(h==0)next
    u<-seq(0,1,length.out=35);s<-u*u*(3-2*u);xx<-t+hw+(t2-t-2*hw)*u;lo<-src[i]+(dst[j]-src[i])*s
    pp<-pp+1;polys[[pp]]<-data.frame(indicator=m,phenotype=k,state=states[i],group=pp,x=c(xx,rev(xx)),y=c(lo,rev(lo+h)));src[i]<-src[i]+h;dst[j]<-dst[j]+h
   }
   stopifnot(max(abs(src-a$ymax))<1e-8,max(abs(dst-b$ymax))<1e-8)
  }
 }
 po<-do.call(rbind,polys);z<-nd[nd$dataset==db,];z$hw<-ifelse(z$indicator=='MSI',.6,.9)
 rowlabs<-c(MSI='MSI (bpm/mmHg)',Na='Na (mmol/L)',K='K (mmol/L)',Cl='Cl (mmol/L)',HCO3='HCO3 (mmol/L)')
 z$indicator<-factor(z$indicator,levels=marks,labels=rowlabs);po$indicator<-factor(po$indicator,levels=marks,labels=rowlabs)
 cc<-cn[cn$dataset==db,];hh<-setNames(paste0(cc$phenotype,' | n = ',format(cc$n,big.mark=','),' (',sprintf('%.1f',cc$pct),'%)'),cc$phenotype)
 p<-ggplot()+geom_polygon(data=po,aes(x,y,group=group,fill=state),alpha=.32,colour=NA)+geom_rect(data=z,aes(xmin=time_hour-hw,xmax=time_hour+hw,ymin=ymin,ymax=ymax,fill=state),colour='white',linewidth=.13)+facet_grid(indicator~phenotype,labeller=labeller(phenotype=hh))+scale_fill_manual(values=sc,breaks=states,labels=c('Low / lower MSI','Reference / middle MSI','High / upper MSI','Missing'))+scale_x_continuous(breaks=c(3,9,15,21),limits=c(0,24),expand=c(0,0))+scale_y_continuous(breaks=c(0,50,100),labels=function(x)paste0(x,'%'),limits=c(0,100),expand=c(0,0))+labs(title=paste0(db,': five-variable state transitions'),subtitle='Columns: fixed C1-C3 phenotypes | Rows: five variables | Full original time resolution',x='Hours after sepsis onset (window midpoint)',y='Patients within each fixed phenotype',caption=sprintf('MSI: all 12 two-hour windows; cutoffs %.3f / %.3f from the MIMIC first-window tertiles. Only MSI bars are wider.\nElectrolytes: four six-hour windows. Reference (mmol/L): Na 135-145; K 3.6-5.2; Cl 98-107; HCO3 22-29.\nGrey: missing; denominator: all classified patients. Ribbons link adjacent windows and use source-state colours; classes remain fixed.',q[1],q[2]))+theme_pub()+theme(panel.grid=element_blank(),strip.text.x=element_text(size=10),strip.text.y=element_text(size=9))
 # Electrolyte labels retained; omit MSI labels to avoid overcrowding narrow bars.
 p<-p+geom_text(data=z[z$indicator!=rowlabs[['MSI']] & z$pct>=8,],aes(time_hour,(ymin+ymax)/2,label=paste0(round(pct),'%')),colour='white',size=2.5,fontface='bold')
 save(p,paste0('State_',db,'五变量状态转换'),12,13)
}

stopifnot(all(md5==unname(tools::md5sum(paths))));csv(data.frame(file=paths,md5_before=md5,md5_after=unname(tools::md5sum(paths))),'输入身份核验')
writeLines(capture.output(sessionInfo()),fn('R运行环境','txt'));cat('V330 COMPLETE\n')
}
main()


