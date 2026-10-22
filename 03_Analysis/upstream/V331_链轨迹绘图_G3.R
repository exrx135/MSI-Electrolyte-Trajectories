# V331: trace diagnostic pages, run after final sensitivity results.
suppressPackageStartupMessages(library(ggplot2))
root<-'Z:/Research/MSI-1';out<-file.path(root,'results','V331_MAR_MNAR敏感性_G3')
pdf<-file.path(out,'V331_MICE链诊断图_G3.pdf');cairo_pdf(pdf,width=13,height=13)
for(db in c('MIMIC','eICU','SICdb')){
 x<-read.csv(file.path(out,paste0('V331_',db,'_链轨迹_G3.csv')));x<-x[is.finite(x$mean),];x$variable<-factor(x$variable,levels=unique(x$variable));mu<-aggregate(mean~variable+iteration,x,mean)
 p<-ggplot(x,aes(iteration,mean,group=interaction(variable,draw)))+geom_line(colour='#438D80',alpha=.18,linewidth=.25)+geom_line(data=mu,aes(iteration,mean,group=variable),inherit.aes=FALSE,colour='#31566E',linewidth=.6)+facet_wrap(~variable,scales='free_y',ncol=4)+labs(title=paste0(db,': imputed-value chain means'),subtitle='20 chains | Thin lines: individual chains; dark line: across-chain mean',x='Iteration',y='Mean of imputed values',caption='Convergence screening uses split-Rhat for means and SDs over the latter half of iterations.\nTrace plots and numerical checks do not establish that the missingness assumptions are true.')+theme_minimal(base_size=9)+theme(panel.grid.minor=element_blank(),strip.text=element_text(face='bold'),plot.title=element_text(size=16,face='bold'),plot.caption=element_text(hjust=0))
 print(p)
}
dev.off()
