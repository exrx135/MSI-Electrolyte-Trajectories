OUT<-"Z:/Research/MSI-1/analysis_updates/Incremental_prediction_20261021"
PRIVATE<-file.path(OUT,"private_patient_level")
DBS<-c("MIMIC","eICU","SICdb")
LABELS<-c(MIMIC="MIMIC-IV",eICU="eICU-CRD",SICdb="SICdb")
MODELS<-c("SOFA","SOFA_plus_class","SOFA_plus_score","SOFA_plus_score_plus_class")
BLUE<-"#0D4483";RED<-"#DD1215"
auc_rank<-function(y,p) {
  n1<-sum(y==1);n0<-sum(y==0)
  (sum(rank(p,ties.method="average")[y==1])-n1*(n1+1)/2)/(n1*n0)
}

# AUC is a ranking functional. Compute it within each validation fold so that
# tiny between-fold intercept shifts cannot artificially break SOFA score ties.
fold_records<-list()
for(db in DBS) {
  v<-readRDS(file.path(PRIVATE,paste0(db,"_OOF_predictions_and_folds.rds")))
  for(r in seq_len(dim(v$predictions_by_repeat)[3]))for(k in 1:10) {
    ix<-v$folds[,r]==k
    for(j in seq_along(MODELS)) {
      fold_records[[length(fold_records)+1]]<-data.frame(dataset=db,model=MODELS[j],
        repeat_id=r,fold=k,n=sum(ix),events=sum(v$outcome[ix]),
        comparable_pairs=sum(v$outcome[ix]==1)*sum(v$outcome[ix]==0),
        auc=auc_rank(v$outcome[ix],v$predictions_by_repeat[ix,j,r]))
    }
  }
}
folds<-do.call(rbind,fold_records)
write.csv(folds,file.path(OUT,"CV_within_fold_AUC_details.csv"),row.names=FALSE)
summary_rows<-list()
for(db in DBS)for(m in MODELS) {
  z<-folds[folds$dataset==db & folds$model==m,]
  repeat_auc<-vapply(split(z,z$repeat_id),function(g)weighted.mean(g$auc,g$comparable_pairs),numeric(1))
  summary_rows[[length(summary_rows)+1]]<-data.frame(dataset=db,model=m,
     auc_cv_within_fold=weighted.mean(z$auc,z$comparable_pairs),
     repeat_auc_min=min(repeat_auc),repeat_auc_max=max(repeat_auc))
}
fold_summary<-do.call(rbind,summary_rows)
v<-read.csv(file.path(OUT,"OOF_performance_and_calibration.csv"),check.names=FALSE)
if("auc"%in%names(v))names(v)[names(v)=="auc"]<-"auc_pooled_mean_oof"
v<-v[,setdiff(names(v),c("auc_cv_within_fold","repeat_auc_min","repeat_auc_max")),drop=FALSE]
v<-merge(v,fold_summary,by=c("dataset","model"),sort=FALSE)
v<-v[order(match(v$dataset,DBS),match(v$model,MODELS)),]
write.csv(v,file.path(OUT,"OOF_performance_and_calibration.csv"),row.names=FALSE)
repeat_pairs<-list()
for(db in DBS)for(cmp in c("SOFA","SOFA_plus_score")) {
  b<-v[v$dataset==db & v$model==cmp,]
  a<-v[v$dataset==db & v$model==paste0(cmp,"_plus_class"),]
  repeat_pairs[[length(repeat_pairs)+1]]<-data.frame(dataset=db,comparison=cmp,
    auc_base_cv=b$auc_cv_within_fold,auc_augmented_cv=a$auc_cv_within_fold,
    delta_auc_cv=a$auc_cv_within_fold-b$auc_cv_within_fold,
    brier_base=b$brier,brier_augmented=a$brier,brier_reduction=b$brier-a$brier,
    calibration_slope_base=b$calibration_slope,calibration_slope_augmented=a$calibration_slope)
}
write.csv(do.call(rbind,repeat_pairs),file.path(OUT,"OOF_paired_comparisons_table_ready.csv"),row.names=FALSE)

bins<-read.csv(file.path(OUT,"OOF_calibration_bins.csv"))
dca<-read.csv(file.path(OUT,"OOF_decision_curve_values.csv"))

setup_layout<-function() {
  par(mfrow=c(2,3),mar=c(4.4,4.6,3.3,1.1),oma=c(3.2,.4,.3,.2),
      family="Arial",cex=.96,col="black",col.axis="black",col.lab="black",
      col.main="black",mgp=c(2.8,.85,0),tcl=-.25,bty="l",las=1)
}
panel_title<-function(db,row) {
  model<-if(row==1)"SOFA" else if(db=="SICdb")"SOFA + SAPS 3" else "SOFA + APS"
  paste0(LABELS[db],"\n",model)
}
outer_legend<-function(calibration=FALSE) {
  par(fig=c(0,1,0,1),new=TRUE,mar=c(0,0,0,0),oma=c(0,0,0,0),xpd=NA)
  plot.new()
  if(calibration)
    legend("bottom",inset=.012,
           legend=c("Base model","Base model + trajectory","Ideal calibration"),
           col=c(BLUE,RED,"#777777"),lty=c(1,1,2),pch=c(16,17,NA),pt.cex=.85,
           lwd=c(1.8,1.8,1.1),horiz=TRUE,bty="n",cex=1.02,x.intersp=1.5)
  else
    legend("bottom",inset=.012,
           legend=c("Base model","Base model + trajectory","Treat all","Treat none"),
           col=c(BLUE,RED,"#999999","#333333"),lty=c(1,1,2,3),
           lwd=c(1.8,1.8,1.1,1.1),horiz=TRUE,bty="n",cex=1.02,x.intersp=1.5)
}
draw_calibration<-function() {
  setup_layout()
  for(row in 1:2)for(db in DBS) {
    models<-if(row==1)MODELS[1:2] else MODELS[3:4]
    plot(NA,xlim=c(0,1),ylim=c(0,1),xaxs="i",yaxs="i",xaxt="n",yaxt="n",
         xlab="Predicted MAKE30 risk",ylab="Observed MAKE30 proportion",
         main=panel_title(db,row),cex.main=1.04)
    axis(1,at=seq(0,1,.2),labels=sprintf("%.1f",seq(0,1,.2)))
    axis(2,at=seq(0,1,.2),labels=sprintf("%.1f",seq(0,1,.2)))
    abline(0,1,col="#999999",lty=2,lwd=1.1)
    for(j in 1:2) {
      z<-bins[bins$dataset==db & bins$model==models[j],];z<-z[order(z$mean_predicted),]
      colour<-c(BLUE,RED)[j]
      segments(z$mean_predicted,z$observed_low,z$mean_predicted,z$observed_high,
               col=adjustcolor(colour,alpha.f=.5),lwd=.8)
      segments(z$mean_predicted-.006,z$observed_low,z$mean_predicted+.006,z$observed_low,
               col=adjustcolor(colour,alpha.f=.5),lwd=.8)
      segments(z$mean_predicted-.006,z$observed_high,z$mean_predicted+.006,z$observed_high,
               col=adjustcolor(colour,alpha.f=.5),lwd=.8)
      lines(z$mean_predicted,z$observed,col=colour,lwd=1.7)
      points(z$mean_predicted,z$observed,pch=c(16,17)[j],col=colour,cex=.72)
    }
  }
  outer_legend(TRUE)
}
draw_dca<-function() {
  setup_layout()
  for(row in 1:2)for(db in DBS) {
    comp<-if(row==1)"SOFA" else "SOFA + additional score"
    z<-dca[dca$dataset==db & dca$comparison==comp,]
    plot(NA,xlim=c(.05,.50),ylim=c(-.025,.37),xaxs="i",yaxs="i",xaxt="n",yaxt="n",
         xlab="Threshold probability",ylab="Net benefit",main=panel_title(db,row),cex.main=1.04)
    axis(1,at=c(.05,.10,.20,.30,.40,.50),labels=c(".05",".10",".20",".30",".40",".50"))
    axis(2,at=c(0,.1,.2,.3),labels=c("0.0","0.1","0.2","0.3"))
    lines(z$threshold,z$nb_all,col="#999999",lty=2,lwd=1.15)
    abline(h=0,col="#333333",lty=3,lwd=1.1)
    lines(z$threshold,z$nb_base,col=BLUE,lwd=1.9)
    lines(z$threshold,z$nb_augmented,col=RED,lwd=1.9)
  }
  outer_legend(FALSE)
}
for(kind in c("Calibration","Decision_curves")) {
  draw<-if(kind=="Calibration")draw_calibration else draw_dca
  stem<-file.path(OUT,paste0(kind,"_internally_validated_MAKE30"))
  ragg::agg_png(paste0(stem,".png"),width=13.5,height=8,units="in",res=200,background="white")
  draw();dev.off()
  svg(paste0(stem,".svg"),width=13.5,height=8,pointsize=12,family="Arial",bg="white")
  draw();dev.off()
  ragg::agg_tiff(paste0(stem,".tif"),width=13.5,height=8,units="in",res=600,
                 background="white",compression="lzw")
  draw();dev.off()
}

# Compact, honest decision-curve summaries: do not imply uniform benefit over
# every threshold. The grid is pre-specified as 5%-50% in steps of one point.
dca_summary<-list()
for(db in DBS)for(comp in c("SOFA","SOFA + additional score")) {
  z<-dca[dca$dataset==db & dca$comparison==comp,]
  for(t in c(.1,.2,.3,.4)) {
    k<-which.min(abs(z$threshold-t))
    dca_summary[[length(dca_summary)+1]]<-data.frame(dataset=db,comparison=comp,
      threshold=t,nb_base=z$nb_base[k],nb_augmented=z$nb_augmented[k],
      delta_nb=z$delta_nb[k],
      added_net_true_positives_per100=100*z$delta_nb[k])
  }
}
write.csv(do.call(rbind,dca_summary),file.path(OUT,"DCA_selected_thresholds_table_ready.csv"),row.names=FALSE)
cat("FINISHED TABLE-READY RESULTS AND SIX FIGURE FILES\n")
print(do.call(rbind,repeat_pairs))
