# RMST extension of the frozen Kaplan-Meier analysis.
# No cohort, event, censoring, class assignment or time origin is altered.
# 14 d is the primary restriction; 21 d is the longer-window comparison.
library(survival)
library(ragg)

out <- "Z:/Research/MSI-1/analysis_updates/RMST_20261021"
dir.create(file.path(out, "diagnostics"), recursive=TRUE, showWarnings=FALSE)
cohorts <- c("MIMIC", "eICU", "SICdb")
cohort_labels <- c(MIMIC="MIMIC-IV", eICU="eICU-CRD", SICdb="SICdb")
pal <- c(MIMIC="#0D4483", eICU="#DD1215", SICdb="#39A08F")
symbols <- c(MIMIC=16, eICU=18, SICdb=17)
taus <- c(14, 21)
z <- qnorm(.975)
csv <- function(d, name) write.csv(d, file.path(out, name), row.names=FALSE, na="")

# Independent direct integration and Greenwood covariance calculation.
# At event time t_j, the tail integral is int[t_j,tau] S_hat(u)du.
# Var(RMST) = sum_j tail_integral_j^2*d_j/[Y_j*(Y_j-d_j)].
direct_rmst <- function(fit, tau) {
  keep <- fit$time < tau
  tt <- fit$time[keep]
  ss <- fit$surv[keep]
  dd <- fit$n.event[keep]
  yy <- fit$n.risk[keep]
  if (!length(tt)) return(c(rmst=tau, variance=0))
  dt <- diff(c(0, tt, tau))
  pieces <- dt * c(1, ss)
  tail_area <- rev(cumsum(rev(pieces)))[-1L]
  event_idx <- which(dd > 0 & yy > dd)
  variance <- sum(tail_area[event_idx]^2 *
                  dd[event_idx] / (yy[event_idx]*(yy[event_idx]-dd[event_idx])))
  # If the risk set is exhausted by deaths, later KM area is zero.
  exhausted <- which(dd > 0 & yy == dd)
  if (length(exhausted)) stopifnot(all(abs(tail_area[exhausted]) < 1e-12))
  c(rmst=sum(pieces), variance=variance)
}

absolute <- contrasts <- risks <- validation <- counts <- list()
curve_list <- list()
idx <- 1L
for (cohort in cohorts) {
  d <- read.csv(file.path(out, "private_patient_level",
                         paste0(cohort, "_frozen_KM.csv")))
  stopifnot(!anyDuplicated(d$unit_id), all(d$eligible),
            all(d$time_days > 0), all(d$death == as.integer(d$event == 1)))
  for (cl in c("C1", "C2", "C3")) {
    b <- d[d$phenotype == cl, ]
    f <- survfit(Surv(time_days, death) ~ 1, data=b, conf.type="log-log")
    stopifnot(all(taus < max(b$time_days)))
    counts[[paste(cohort,cl)]] <- data.frame(
      cohort=cohort, trajectory_class=cl, n=nrow(b), deaths=sum(b$death),
      live_discharge_censored=sum(b$event==2),
      administratively_censored=sum(b$event==0),
      min_time=min(b$time_days), max_time=max(b$time_days))
    risks[[paste(cohort,cl)]] <- data.frame(
      cohort=cohort, trajectory_class=cl, time=c(0,7,14,21,25,28),
      n_risk=sapply(c(0,7,14,21,25,28), function(t) sum(b$time_days >= t)))
    curve_list[[paste(cohort,cl)]] <- data.frame(
      cohort=cohort, trajectory_class=cl,
      time=c(0,f$time), survival=c(1,f$surv))
    for (tau in taus) {
      tab <- summary(f, rmean=tau)$table
      mu <- unname(tab["rmean"])
      se <- unname(tab["se(rmean)"])
      direct <- direct_rmst(f,tau)
      err_mu <- abs(mu-direct["rmst"])
      err_var <- abs(se^2-direct["variance"])
      stopifnot(is.finite(mu), is.finite(se),
                err_mu < 1e-10, err_var < 1e-10,
                mu >= 0, mu <= tau)
      absolute[[idx]] <- data.frame(
        cohort=cohort, cohort_label=unname(cohort_labels[cohort]),
        tau_days=tau, trajectory_class=cl, n=nrow(b),
        deaths_by_tau=sum(b$death==1 & b$time_days<=tau),
        n_risk_at_tau=sum(b$time_days >= tau),
        rmst_days=mu, standard_error=se,
        lower_95=mu-z*se, upper_95=mu+z*se,
        restricted_mean_time_lost=tau-mu)
      validation[[idx]] <- data.frame(
        cohort=cohort, trajectory_class=cl, tau_days=tau,
        direct_area_error=err_mu, direct_variance_error=err_var)
      idx <- idx+1L
    }
  }
}
ab <- do.call(rbind,absolute)
ci <- 1L
for (tau in taus) for (cohort in cohorts) {
  a <- ab[ab$cohort==cohort & ab$tau_days==tau, ]
  ref <- a[a$trajectory_class=="C1", ]
  for (cl in c("C2","C3")) {
    x <- a[a$trajectory_class==cl, ]
    delta <- x$rmst_days-ref$rmst_days
    se <- sqrt(x$standard_error^2+ref$standard_error^2)
    contrasts[[ci]] <- data.frame(
      cohort=cohort, cohort_label=unname(cohort_labels[cohort]),
      tau_days=tau, comparison=paste(cl,"vs C1"),
      difference_days=delta, standard_error=se,
      lower_95=delta-z*se, upper_95=delta+z*se,
      z_value=delta/se, p_value=2*pnorm(-abs(delta/se)))
    ci <- ci+1L
  }
}
con <- do.call(rbind,contrasts)
con$BH_q_12_contrasts_audit <- p.adjust(con$p_value,method="BH")
ab <- ab[order(match(ab$cohort,cohorts),ab$tau_days,ab$trajectory_class),]
con <- con[order(con$tau_days,con$comparison,match(con$cohort,cohorts)),]
rownames(ab) <- rownames(con) <- NULL
csv(ab,"Absolute_RMST_by_cohort_class.csv")
csv(con,"RMST_contrasts_C2_C3_vs_C1.csv")
csv(do.call(rbind,risks),"Risk_set_support.csv")
csv(do.call(rbind,counts),"Frozen_KM_cohort_counts.csv")
csv(do.call(rbind,validation),"diagnostics/Independent_Greenwood_validation.csv")
csv(do.call(rbind,curve_list),"diagnostics/KM_step_curves_for_integral.csv")

# Concise table-ready layout: one row per cohort and restriction time.
rows <- list()
for (cohort in cohorts) for (tau in taus) {
  a <- ab[ab$cohort==cohort & ab$tau_days==tau,]
  c <- con[con$cohort==cohort & con$tau_days==tau,]
  fmt <- function(mu,lo,hi) sprintf("%.2f (%.2f to %.2f)",mu,lo,hi)
  val <- function(cl) {
    k <- a[a$trajectory_class==cl,]
    fmt(k$rmst_days,k$lower_95,k$upper_95)
  }
  cmp <- function(cl) {
    k <- c[c$comparison==paste(cl,"vs C1"),]
    fmt(k$difference_days,k$lower_95,k$upper_95)
  }
  pv <- function(cl) c$p_value[c$comparison==paste(cl,"vs C1")]
  rows[[paste(cohort,tau)]] <- data.frame(
    cohort=unname(cohort_labels[cohort]),tau_days=tau,
    C1_RMST_95CI=val("C1"),C2_RMST_95CI=val("C2"),C3_RMST_95CI=val("C3"),
    C2_minus_C1_95CI=cmp("C2"),C2_vs_C1_P=pv("C2"),
    C3_minus_C1_95CI=cmp("C3"),C3_vs_C1_P=pv("C3"))
}
csv(do.call(rbind,rows),"RMST_table_ready.csv")

plot_rmst <- function() {
  op <- par(no.readonly=TRUE)
  on.exit(par(op))
  par(mfrow=c(1,2), mar=c(4.2,4.6,2.6,1), oma=c(0,0,2.2,0),
      mgp=c(2.4,.65,0), family="sans", las=1,
      cex.axis=.95, cex.lab=1.05, cex.main=1.1, xpd=FALSE)
  xlim <- c(floor(min(con$lower_95)*2)/2-.10,.12)
  for (tau in taus) {
    plot(NA,xlim=xlim,ylim=c(.5,2.5),xlab="RMST difference (days)",
         ylab="",axes=FALSE,main=paste0(tau," days after the landmark"),
         xaxs="i",yaxs="i")
    rect(xlim[1],.5,xlim[2],1.5,col="#F7F7F7",border=NA)
    abline(v=0,col="#777777",lty=2,lwd=.8)
    axis(1,at=seq(ceiling(xlim[1]),0,by=.5),tck=-.02,col="#333333",col.axis="#111111")
    axis(2,at=c(2,1),labels=c("C2 vs C1","C3 vs C1"),tick=FALSE,line=.1)
    box(bty="l",col="#333333",lwd=.8)
    for (j in seq_along(cohorts)) {
      cohort <- cohorts[j]
      for (cl in c("C2","C3")) {
        v <- con[con$cohort==cohort & con$tau_days==tau &
                 con$comparison==paste(cl,"vs C1"),]
        y <- (if (cl=="C2") 2 else 1) + c(.15,0,-.15)[j]
        segments(v$lower_95,y,v$upper_95,y,col=pal[cohort],lwd=1.55)
        segments(c(v$lower_95,v$upper_95),y-.028,
                 c(v$lower_95,v$upper_95),y+.028,col=pal[cohort],lwd=1.15)
        points(v$difference_days,y,pch=symbols[cohort],col=pal[cohort],cex=.92)
      }
    }
  }
  par(fig=c(0,1,0,1),new=TRUE,mar=c(0,0,0,0),oma=c(0,0,0,0),xpd=NA)
  plot.new()
  legend("top",inset=.005,legend=unname(cohort_labels),col=unname(pal),
         pch=unname(symbols),lty=1,lwd=1.4,pt.cex=.85,
         bty="n",horiz=TRUE,cex=1.02,x.intersp=.7,seg.len=1.6)
}
stem <- file.path(out,"Restricted_mean_hospital_survival_differences")
agg_png(paste0(stem,".png"),width=2600,height=950,res=220,background="white")
plot_rmst()
dev.off()
svg(paste0(stem,".svg"),width=11.818,height=4.318,bg="white")
plot_rmst()
dev.off()
agg_tiff(paste0(stem,".tif"),width=7090,height=2590,res=600,
         background="white",compression="lzw")
plot_rmst()
dev.off()
capture.output(sessionInfo(),file=file.path(out,"diagnostics/session_info.txt"))
print(ab)
print(con)
