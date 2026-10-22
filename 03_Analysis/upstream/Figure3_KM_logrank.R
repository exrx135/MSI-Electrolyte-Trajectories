# User-approved KM analysis, 2026-10-14.
# Reuses adjudicated eligible patients; live discharge is censored.
library(survival)
library(ggplot2)
library(patchwork)
out <- 'Z:/Research/MSI-1/results/生存分析_三库图与代码'
pal <- c(C1='#31566E',C2='#C56E52',C3='#438D80')
titles <- c(MIMIC='A  MIMIC-IV',eICU='B  eICU-CRD',SICdb='C  SICdb')
plots <- tests <- checks <- list()
csv <- function(x,n) write.csv(x,file.path(out,paste0(n,'.csv')),row.names=FALSE,fileEncoding='UTF-8')
for(db in names(titles)) {
 a <- read.csv(file.path(out,paste0(db,'_竞争风险分析数据.csv')))
 stopifnot(!anyDuplicated(a$unit_id),all(a$eligible),all(a$time_days>0),all(a$event %in% 0:2))
 a$death <- as.integer(a$event==1)
 a$phenotype <- factor(a$phenotype,levels=names(pal))
 csv(a,paste0(db,'_KM分析数据'))
 lr <- survdiff(Surv(time_days,death)~phenotype,data=a,rho=0)
 df <- length(lr$n)-1
 pv <- pchisq(lr$chisq,df,lower.tail=FALSE)
 tests[[db]] <- data.frame(dataset=db,n=nrow(a),deaths=sum(a$death),chisq=lr$chisq,df=df,p=pv,
                          test='Overall three-group log-rank; live discharge censored')
 label <- if(pv<.001) 'Log-rank P < 0.001' else sprintf('Log-rank P = %.3f',pv)
 curves <- risks <- list()
 for(k in names(pal)) {
  b <- a[a$phenotype==k,]
  f <- survfit(Surv(time_days,death)~1,data=b,conf.type='log-log')
  cv <- data.frame(phenotype=k,time=c(0,f$time),survival=c(1,f$surv),lower=c(1,f$lower),upper=c(1,f$upper),n_censor=c(0,f$n.censor))
  # Risk set immediately before t: includes events/censors occurring at t.
  risk_times <- c(0,5,15,20,25)
  risks[[k]] <- data.frame(phenotype=k,time=risk_times,n_risk=sapply(risk_times,function(t) sum(b$time_days>=t)))
  stopifnot(sum(f$n.censor)==sum(b$death==0),risks[[k]]$n_risk[1]==nrow(b))
  manual <- cumprod(1-f$n.event/f$n.risk)
  stopifnot(max(abs(manual-f$surv))<1e-10,all(diff(cv$survival)<=1e-10))
  checks[[paste(db,k)]] <- data.frame(dataset=db,phenotype=k,n=nrow(b),deaths=sum(b$death),
    last_observed_day=max(b$time_days),last_KM=tail(f$surv,1),manual_max_error=max(abs(manual-f$surv)))
  curves[[k]] <- cv
 }
 d <- do.call(rbind,curves);csv(d,paste0(db,'_KM曲线及95CI'))
 ri <- do.call(rbind,risks);csv(ri,paste0(db,'_KM风险人数'))
 csv(d[d$n_censor>0,c('phenotype','time','survival','n_censor')],paste0(db,'_KM删失标记'))
 leg <- data.frame(phenotype=names(pal),n=as.integer(table(a$phenotype)),y=c(.556,.514,.472))
 leg$label <- sprintf('%s (n = %s)',leg$phenotype,format(leg$n,big.mark=',',trim=TRUE))
 p <- ggplot(d,aes(time,survival,colour=phenotype))+
  geom_step(linewidth=.9,show.legend=FALSE)+
  annotate('rect',xmin=2.3,xmax=19.5,ymin=.439,ymax=.589,fill='white',colour='#777777',linewidth=.35)+
  geom_segment(data=leg,aes(x=3.2,xend=5.2,y=y,yend=y,colour=phenotype),inherit.aes=FALSE,linewidth=1,show.legend=FALSE)+
  geom_text(data=leg,aes(x=6,y=y,label=label),inherit.aes=FALSE,hjust=0,size=3.6,colour='#222222')+
  annotate('text',x=29.4,y=.982,label=label,hjust=1,vjust=1,size=3.8)+
  scale_colour_manual(values=pal)+
  scale_x_continuous(breaks=seq(0,30,5),expand=expansion(mult=c(.015,.015)))+
  scale_y_continuous(breaks=seq(.4,1,.1),labels=function(x)sprintf('%.1f',x),expand=expansion(mult=c(0,.025)))+
  coord_cartesian(xlim=c(0,30),ylim=c(.4,1))+
  labs(title=titles[[db]],x=NULL,y=if(db=='MIMIC') 'Survival probability' else NULL)+
  theme_classic(base_size=12,base_family='sans')+
  theme(plot.title=element_text(face='bold',size=14,margin=margin(b=12)),axis.title=element_text(size=11),
        axis.title.x=element_text(margin=margin(t=9)),axis.title.y=element_text(margin=margin(r=8)),
        axis.text=element_text(colour='#222222'),axis.line=element_line(linewidth=.45),plot.margin=margin(10,14,10,10))
 pr <- ggplot(ri,aes(time,phenotype,label=n_risk,colour=phenotype))+
  geom_text(size=3.2,show.legend=FALSE)+scale_colour_manual(values=pal)+
  scale_x_continuous(breaks=risk_times,expand=expansion(mult=c(.015,.015)))+
  scale_y_discrete(limits=c('C3','C2','C1'),expand=expansion(add=.55))+
  coord_cartesian(xlim=c(0,30),clip='off')+
  labs(title='Number at risk',x='Days after the 24-h landmark',y=NULL)+
  theme_classic(base_size=11,base_family='sans')+
  theme(axis.line=element_blank(),axis.ticks.y=element_blank(),axis.text=element_text(colour='#222222'),
        axis.text.y=element_text(margin=margin(r=14)),
        axis.title.x=element_text(size=11,margin=margin(t=8)),plot.title=element_text(size=11),
        plot.margin=margin(0,14,10,10))
 plots[[db]] <- p/pr+plot_layout(heights=c(4.8,1.3))
}
csv(do.call(rbind,tests),'KM_总体logrank检验')
csv(do.call(rbind,checks),'KM_人数与数值核验')
combined <- wrap_plots(plots,ncol=3)
stem <- file.path(out,'Figure3_KM_logrank_三库_风险人数表')
ggsave(paste0(stem,'.png'),combined,width=15,height=6.1,dpi=400,bg='white')
ggsave(paste0(stem,'.svg'),combined,width=15,height=6.1,device=grDevices::svg,bg='white')
ggsave(paste0(stem,'.pdf'),combined,width=15,height=6.1,device=grDevices::cairo_pdf,bg='white')
writeLines(c('当前选定版本：Figure3_KM_logrank_三库_风险人数表；用户于2026-10-14明确改选KM与log-rank。',
 '按用户要求不显示删失标记；风险表仅列0、5、15、20、25天发生事件或删失之前的在险人数。删失仍计入KM和log-rank分析。',
 'Figure legend待写：Censoring marks are omitted for visual clarity; numbers at risk at selected time points are shown below each panel.',
 '风险人数表展示仍在险的人数，不单独展示删失数量或每例删失时间。',
 '使用既有纳排及死亡判定后的患者数据。院内死亡为事件，存活出院在出院时删失，行政终点在入ICU后30天删失。',
 '观察起点维持sepsis起点后24小时landmark；横轴为landmark后天数，不外推至每人均有30天随访。',
 '各库P值为C1/C2/C3总体log-rank检验（rho=0，df=2），不是两两比较，也未调整协变量。',
 '纵轴为KM生存概率，显示范围0.4–1；不再是1－AJ累积发生率。存活出院删失需要独立删失假设；曲线不能解释为实际院内死亡累积风险的补数。',
 '保留原AJ版本供追溯，原MAKE30、分型及其他冻结分析未重跑。',
 '原纳排、SICdb出院时间字段局限仍适用，详见使用说明与Figure_legend待写要点.md。',
 'eICU 35例及SICdb 5例结局缺失的说明留待论文figure legend，不写入图片。'),
 file.path(out,'KM版本说明_当前选定.txt'),useBytes=TRUE)
print(do.call(rbind,tests));print(do.call(rbind,checks))
