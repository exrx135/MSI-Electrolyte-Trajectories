# ============================================================
# V316 STEP 5, MIMIC-IV development + eICU-CRD/SICdb transport
# ============================================================
# Rscript V316_STEP5_transport_make30_G3.R "D:/STEP5_G3_INPUT"
# Or place this script beside the seven input files and run:
# source(file.choose())
# Requires R and gbmt 0.1.4. No refitting, no recalibration, no G selection.
# All input file MD5s correspond to the uploaded September 26 cohort exports
# and the V295 G3 S01 fit, except the documented eICU V316 time conversion/exclusions.
# G3 assignments and the 5 MIMIC z parameters remain frozen.
# Outputs are written under GBMTM_stage5_output in an input-specific V316 folder; existing
# finished outputs are read back. A failed run can be restarted safely.
# Computational missing-data handling is the same gbmt 0.1.4 imputFun used by
# V295. No NA in the original trajectory table is changed or filled in.

options(stringsAsFactors=FALSE, scipen=999)
VERSION <- "V316"
EXPECTED <- c(
  "MIMIC_FINAL_MAKE30_ELECTROLYTEQC.csv"="f1eb2d7fac1303a682b53c679dde9abd",
  "SICdb_FINAL_MAKE30_ELECTROLYTEQC.csv"="d52a3ecdbb636bba2859eb9f102fdf45",
  "eICU_FINAL_RAW_ELECTROLYTEQC.csv"="a06fddbcbc3b93f2c83227201c1138cd",
  "V295_compact_fit_G3_S01.rds"="cf5881479ee206195aa80d59fa3bd17f",
  "V295_patient_class_G3_S01.csv"="dd12a3730af44c1035c79f4558c14db9",
  "V295_trajectory_G3_S01.csv"="3d47aa52a546266e937f50a26de2676a",
  "V295_task_summary_G3_S01.csv"="6b8eb067bb4013a54f6cb4dafe6d5f83"
)
Z <- c("msi_z","sodium_z","potassium_z","chloride_z","bicarbonate_z")
FROZEN_MEAN <- setNames(c(1.229422,138.931869,4.221718,
                          105.125243,22.905997),Z)
FROZEN_SD <- setNames(c(0.3219380,5.3785897,0.6426219,
                        6.5980074,4.7412984),Z)
G <- 3L

if (!requireNamespace("gbmt",quietly=TRUE))
  stop("gbmt 0.1.4 is required: install.packages('gbmt') on the laptop.")
if (as.character(utils::packageVersion("gbmt"))!="0.1.4")
  stop("Installed gbmt version differs from the frozen V295 engine (0.1.4).")

profile <- Sys.getenv("USERPROFILE", unset="")
if (!nzchar(profile) || !dir.exists(profile)) {
  profile <- path.expand("~")
  if (identical(basename(profile), "Documents"))
    stop("USERPROFILE is unavailable and ~ points to Documents.")
}
output_root <- file.path(profile,"Desktop","Research","MSI",
                         "提取轨迹文件","GBMTM_stage5_output")
if (!dir.exists(output_root)) stop("Exact STEP5 output root does not exist: ",output_root)
args <- commandArgs(trailingOnly=TRUE)
input_dir <- if(length(args)) args[1L] else file.path(output_root,"V307_STEP5_laptop_kit_G3")
input_dir <- normalizePath(input_dir,winslash="/",mustWork=TRUE)
files <- file.path(input_dir,names(EXPECTED))
if (any(!file.exists(files)))
  stop("Missing Step 5 inputs: ",paste(names(EXPECTED)[!file.exists(files)],collapse=", "))
actual <- unname(tools::md5sum(files))
names(actual) <- names(EXPECTED)
if (any(is.na(actual)) || any(tolower(actual)!=EXPECTED))
  stop("One or more Step 5 inputs differ from the uploaded frozen files: ",
       paste(names(EXPECTED)[tolower(actual)!=EXPECTED | is.na(actual)],collapse=", "))
input_key <- paste0(substr(actual["V295_compact_fit_G3_S01.rds"],1L,8L),"_",
                    substr(actual["MIMIC_FINAL_MAKE30_ELECTROLYTEQC.csv"],1L,8L),"_",
                    substr(actual["eICU_FINAL_RAW_ELECTROLYTEQC.csv"],1L,8L))
out_dir <- file.path(output_root,paste0(VERSION,"_STEP5_results_G3_",input_key))
dir.create(out_dir,recursive=TRUE,showWarnings=FALSE)
if (!dir.exists(out_dir)) stop("Cannot create output directory: ",out_dir)

write_csv_once <- function(d,name) {
  target <- file.path(out_dir,paste0(VERSION,"_",name,"_G3.csv"))
  if (file.exists(target)) return(invisible(target))
  tmp <- paste0(target,".tmp_",Sys.getpid())
  utils::write.csv(d,tmp,row.names=FALSE,na="",fileEncoding="UTF-8")
  if (!file.rename(tmp,target)) stop("Cannot save: ",target)
  invisible(target)
}
write_text_once <- function(lines,name) {
  target <- file.path(out_dir,paste0(VERSION,"_",name,"_G3.txt"))
  if (file.exists(target)) return(invisible(target))
  tmp <- paste0(target,".tmp_",Sys.getpid())
  writeLines(lines,tmp,useBytes=TRUE)
  if (!file.rename(tmp,target)) stop("Cannot save: ",target)
  invisible(target)
}
read_cohort <- function(filename,id_col) {
  utils::read.csv(file.path(input_dir,filename),check.names=FALSE,
                  colClasses=setNames("character",id_col),
                  na.strings=c("","NA","NULL"),stringsAsFactors=FALSE)
}
cohort <- list(
  MIMIC=read_cohort("MIMIC_FINAL_MAKE30_ELECTROLYTEQC.csv","stay_id"),
  SICdb=read_cohort("SICdb_FINAL_MAKE30_ELECTROLYTEQC.csv","case_id"),
  eICU=read_cohort("eICU_FINAL_RAW_ELECTROLYTEQC.csv","patientunitstayid")
)
id_col <- c(MIMIC="stay_id",SICdb="case_id",eICU="patientunitstayid")
expected_n <- c(MIMIC=9804L,SICdb=3831L,eICU=5839L)

fit <- readRDS(file.path(input_dir,"V295_compact_fit_G3_S01.rds"))
ref <- read_cohort("V295_patient_class_G3_S01.csv","unit_id")
traj <- utils::read.csv(file.path(input_dir,"V295_trajectory_G3_S01.csv"))
task <- utils::read.csv(file.path(input_dir,"V295_task_summary_G3_S01.csv"))
if (!inherits(fit,"gbmt") || fit$call$ng!=G || fit$call$d!=2L ||
    fit$call$scaling!=0 || !isTRUE(fit$call$pruning) ||
    !identical(fit$call$x.names,Z) || length(fit$beta)!=G ||
    length(fit$Sigma)!=G || length(fit$prior)!=G ||
    nrow(as.matrix(fit$Xmat))!=12L ||
    nrow(task)!=1L || task$status!="CONVERGED_VALID_FULL_G" ||
    task$requested_g!=G || task$start_id!=1L ||
    task$input_md5!="0361d8b59184b7a6b0a5fa9c946088a0")
  stop("V295 G3 S01 fitted model or result identity differs from frozen setup.")
if (nrow(ref)!=expected_n["MIMIC"] || anyDuplicated(ref$unit_id) ||
    !setequal(as.integer(ref$assigned_class),seq_len(G)) ||
    !all(paste0("p_class",1:G)%in%names(ref)))
  stop("MIMIC V295 patient classification is incomplete.")
times <- as.numeric(rownames(as.matrix(fit$Xmat)))
hours <- seq(1,23,by=2)
if (length(times)!=12L || anyNA(times) || anyDuplicated(times) ||
    is.unsorted(times,strictly=TRUE) ||
    max(abs(as.matrix(fit$Xmat)-cbind(1,(times-min(times)+1),
                         (times-min(times)+1)^2)))>1e-7 ||
    !all(vapply(hours,function(h) any(abs(traj$time_hour_center-h)<1e-8),logical(1))) ||
    max(abs(sort(unique(traj$time_model))-times))>1e-9)
  stop("V295 fitted design does not match the 12 fixed 0-24 h time anchors.")

for (db in names(cohort)) {
  d <- cohort[[db]]
  if(nrow(d)!=expected_n[db] ||
     anyNA(d[[id_col[db]]]) || anyDuplicated(d[[id_col[db]]]))
    stop(db,": cohort count/ID mismatch.")
  required <- c(paste0("msi_m",1:12),
                unlist(lapply(c("sodium","potassium","chloride","bicarbonate"),
                              function(v) paste0(v,"_t",1:4))),
                "make30","make30_evaluable_window","make30_death",
                "make30_new_krt","make30_persistent_kidney_dysfunction",
                "age","sex","sepsis_sofa_score","baseline_egfr_mdrd",
                "ckd_or_renal_disease","diabetes",
                "mechanical_ventilation_0_24h","vasopressor_0_24h")
  if (!all(required%in%names(d))) stop(db,": essential variables missing.")
  if (any(!is.na(d$make30) & !(d$make30%in%c(0,1))) ||
      any(is.na(d$make30_evaluable_window))) stop(db,": invalid MAKE30 coding.")
}
if (!setequal(cohort$MIMIC$stay_id,ref$unit_id))
  stop("MIMIC MAKE30 cohort IDs do not match frozen V295 training patients.")

# Recalculate development means/SDs from RAW observed windows only.
for (j in seq_along(Z)) {
  prefix <- c("msi_m","sodium_t","potassium_t","chloride_t","bicarbonate_t")[j]
  count <- if (j==1L) 12L else 4L
  observed <- unlist(cohort$MIMIC[paste0(prefix,1:count)],use.names=FALSE)
  observed <- as.numeric(observed[!is.na(observed)])
  if(abs(mean(observed)-FROZEN_MEAN[j])>0.00001 ||
     abs(stats::sd(observed)-FROZEN_SD[j])>0.00001)
    stop("MIMIC frozen standardization does not match raw export: ",prefix)
}

make_panel <- function(d,db) {
  ids <- as.character(d[[id_col[db]]])
  ii <- rep(seq_along(ids),each=12L)
  tt <- rep(seq_len(12L),times=length(ids))
  raw <- data.frame(unit_id=ids[ii],time_model=times[tt],
                    msi_z=as.numeric(as.matrix(d[paste0("msi_m",1:12)])[cbind(ii,tt)]),
                    sodium_z=NA_real_,potassium_z=NA_real_,chloride_z=NA_real_,
                    bicarbonate_z=NA_real_,check.names=FALSE)
  # t1..t4 were measured at ~3,9,15,21 h: MSI rows 2,5,8,11.
  match_t <- match(tt,c(2L,5L,8L,11L))
  for (j in 2:5) {
    prefix <- c("sodium","potassium","chloride","bicarbonate")[j-1L]
    sub <- as.matrix(d[paste0(prefix,"_t",1:4)])
    rr <- which(!is.na(match_t))
    raw[[Z[j]]][rr] <- as.numeric(sub[cbind(ii[rr],match_t[rr])])
  }
  for (v in Z) raw[[v]] <- (raw[[v]]-FROZEN_MEAN[v])/FROZEN_SD[v]
  raw <- raw[order(raw$unit_id,raw$time_model),,drop=FALSE]
  rownames(raw) <- NULL
  if (nrow(raw)!=12L*nrow(d) || anyDuplicated(raw[c("unit_id","time_model")]))
    stop(db,": unexpected panel dimension/keys.")
  # The V295 model requires at least one observation per indicator per stay.
  counts <- stats::aggregate(!is.na(raw[Z]),list(unit_id=raw$unit_id),sum)
  if (any(as.matrix(counts[Z])==0L))
    stop(db,": a patient has no observations for at least one trajectory marker.")
  raw
}

# Same V295 scoring pieces as the Step 4A EM; raw missing values stay NA.
imputFun <- getFromNamespace("imputFun","gbmt")
logdmN <- getFromNamespace("logdmN","gbmt")
pldBound <- getFromNamespace("pldBound","gbmt")
transport <- function(panel,db) {
  ids <- unique(panel$unit_id)
  tix <- match(panel$time_model,times)
  if (anyNA(tix)) stop(db,": unsupported time anchor.")
  score <- matrix(NA_real_,length(ids),G)
  for(k in seq_len(G)) {
    beta <- as.matrix(fit$beta[[k]])
    S <- as.matrix(fit$Sigma[[k]])
    if(!identical(colnames(beta),Z) || !identical(rownames(S),Z) ||
       !identical(colnames(S),Z) || nrow(beta)!=3L ||
       !is.finite(fit$prior[k]) || fit$prior[k]<=0)
      stop("Frozen G3 class parameters are invalid: ",k)
    completed <- imputFun(xdat=panel[Z],tval=panel$time_model-min(times)+1,
                          beta=beta,S=S,d=2L)
    if(nrow(completed)!=nrow(panel) || anyNA(completed[Z]))
      stop(db,": gbmt conditional missing-data step failed for class ",k)
    mu <- as.matrix(fit$Xmat)%*%beta
    D <- determinant(S,logarithm=TRUE)
    if(D$sign!=1) stop(db,": nonpositive class covariance.")
    diff <- as.matrix(completed[Z])-mu[tix,,drop=FALSE]
    density <- -0.5*(length(Z)*base::log(2*pi)+as.numeric(D$modulus)+
                       rowSums((diff%*%solve(S))*diff))
    if (any(!is.finite(density))) stop(db,": nonfinite class log density.")
    # Verify the fast matrix calculation against the installed gbmt version.
    for (i in unique(c(1L,length(ids)%/%2L,length(ids)))) {
      rr <- which(panel$unit_id==ids[i])
      check <- logdmN(as.matrix(completed[rr,Z,drop=FALSE]),mu,S)
      if (length(check)!=length(rr) || anyNA(check) || any(!is.finite(check)) ||
          max(abs(check-density[rr]))>1e-6)
        stop(db,": fast score is inconsistent with gbmt::logdmN.")
    }
    totals <- rowsum(matrix(density,ncol=1L),group=panel$unit_id,reorder=FALSE)
    if(!identical(rownames(totals),ids)) stop(db,": scorer ID order mismatch.")
    score[,k] <- pldBound(as.numeric(totals[,1L])+base::log(fit$prior[k]))
    rm(completed,diff,density,totals)
    invisible(gc(verbose=FALSE))
    message(db," class ",k,"/",G," scored")
  }
  score <- exp(score-apply(score,1L,max))
  score <- score/rowSums(score)
  if (anyNA(score) || any(!is.finite(score)) ||
      max(abs(rowSums(score)-1))>1e-10)
    stop(db,": class posterior is invalid.")
  out <- data.frame(unit_id=ids,assigned_class=max.col(score,ties.method="first"),
                    max_posterior=apply(score,1L,max),check.names=FALSE)
  for(k in seq_len(G)) out[[paste0("p_class",k)]] <- score[,k]
  out
}

ari <- function(x,y) {
  tab <- table(x,y)
  pairs <- function(z) sum(z*(z-1)/2)
  n <- length(x)
  rowp <- pairs(rowSums(tab)); colp <- pairs(colSums(tab))
  expval <- rowp*colp/(n*(n-1)/2)
  denom <- (rowp+colp)/2-expval
  if(denom==0) return(NA_real_)
  (pairs(tab)-expval)/denom
}

# MIMIC replay must run BEFORE looking at MAKE30. V295 stores posteriors from
# immediately before pruning and beta/Sigma from immediately after pruning.
# Record both and stop outcome modeling if the mismatch is material.
replay_file <- file.path(out_dir,"V316_MIMIC_final_parameter_replay_G3.csv")
if(file.exists(replay_file)) {
  mimic_score <- utils::read.csv(replay_file,colClasses=c(unit_id="character"))
} else {
  mimic_score <- transport(make_panel(cohort$MIMIC,"MIMIC"),"MIMIC")
  write_csv_once(mimic_score,"MIMIC_final_parameter_replay")
}
if(nrow(mimic_score)!=nrow(ref) || !setequal(mimic_score$unit_id,ref$unit_id))
  stop("Saved MIMIC replay file failed identity check.")
ref_ord <- ref[match(mimic_score$unit_id,ref$unit_id),,drop=FALSE]
training_agreement <- mean(mimic_score$assigned_class==ref_ord$assigned_class)
training_ari <- ari(mimic_score$assigned_class,ref_ord$assigned_class)
ref_p <- as.matrix(ref_ord[paste0("p_class",1:G)])
new_p <- as.matrix(mimic_score[paste0("p_class",1:G)])
audit <- data.frame(saved_v295_class_agreement=training_agreement,
                    saved_v295_ARI=training_ari,
                    posterior_mean_abs_diff=mean(abs(ref_p-new_p)),
                    posterior_max_abs_diff=max(abs(ref_p-new_p)),
                    model_input_md5=task$input_md5,
                    fit_md5=actual["V295_compact_fit_G3_S01.rds"],
                    gbmt_version=as.character(utils::packageVersion("gbmt")))
write_csv_once(audit,"MIMIC_training_replay_audit")
transition <- as.data.frame(table(
  saved_class=factor(ref_ord$assigned_class,levels=1:G),
  final_parameter_class=factor(mimic_score$assigned_class,levels=1:G)))
write_csv_once(transition,"MIMIC_saved_x_final_parameter_transition")
message("Training replay: agreement=",round(training_agreement,5),
        "; ARI=",round(training_ari,5),"; posterior MAE=",
        round(audit$posterior_mean_abs_diff,5))
if(!is.finite(training_ari) || training_agreement<0.95 || training_ari<0.90)
  stop("Model transport preflight: final pruned parameters materially change ",
       "the frozen MIMIC V295 class assignments. Audit CSV is saved. ",
       "Resolve parameter/assignment consistency BEFORE MAKE30 analysis.")

# In MIMIC preserve the official saved Step 4A classes; in external cohorts
# use only the final MIMIC G3 parameter posterior above. No outcomes enter it.
scored <- list(MIMIC=ref[,c("unit_id","assigned_class",paste0("p_class",1:G))])
for(db in c("SICdb","eICU")) {
  dest <- file.path(out_dir,paste0(VERSION,"_",db,"_parameter_transport_G3.csv"))
  if(file.exists(dest)) {
    scored[[db]] <- utils::read.csv(dest,colClasses=c(unit_id="character"))
  } else {
    scored[[db]] <- transport(make_panel(cohort[[db]],db),db)
    write_csv_once(scored[[db]],paste0(db,"_parameter_transport"))
  }
  if(nrow(scored[[db]])!=nrow(cohort[[db]]) ||
     !setequal(scored[[db]]$unit_id,cohort[[db]][[id_col[db]]]))
    stop(db,": saved external posterior IDs are invalid.")
}

check_outcome <- function(d,db) {
  death <- as.integer(d$make30_death)
  krt <- as.integer(d$make30_new_krt)
  renal <- as.integer(d$make30_persistent_kidney_dysfunction)
  derived <- ifelse(death==1L | krt==1L | renal==1L,1L,
                    ifelse(is.na(death)|is.na(krt)|is.na(renal),NA_integer_,0L))
  if(any((is.na(derived)!=is.na(d$make30)) |
         (!is.na(derived)&derived!=d$make30),na.rm=TRUE))
    stop(db,": MAKE30 conflicts with its components.")
  if(any(!is.na(d$make30)&d$make30_evaluable_window!=1L))
    stop(db,": a known outcome is outside the eligible landmark window.")
  good <- d$make30_evaluable_window==1L & !is.na(d$make30)
  if(db=="MIMIC") {
    beg <- as.POSIXct(d$sepsis_time,format="%d/%m/%Y %H:%M:%S",tz="UTC")
    end <- as.POSIXct(d$landmark_time,format="%d/%m/%Y %H:%M:%S",tz="UTC")
    if(anyNA(beg)||anyNA(end)||
       max(abs(as.numeric(difftime(end,beg,units="hours"))-24))>0.00001 ||
       any(d$make30_risk_start!=d$landmark_time))
      stop("MIMIC 24-hour landmark/risk time inconsistent.")
  } else if(db=="SICdb") {
    if(any(d$landmark_offset_sec-d$sepsis_offset_sec!=86400) ||
       any(d$make30_risk_start_offset_sec!=d$landmark_offset_sec))
      stop("SICdb 24-hour landmark/risk time inconsistent.")
  } else {
    if(any(d$landmark_offset-d$sepsis_offset!=1440))
      stop("eICU 24-hour landmark inconsistent.")
    delta <- (d$make30_censor_hosp_offset-
              d$make30_risk_start_hosp_offset)/1440
    if(anyNA(delta) || max(abs(delta-d$make30_followup_days_after_landmark))>1e-6)
      stop("eICU hospital time-base follow-up inconsistent.")
    shift <- d$make30_risk_start_hosp_offset-d$landmark_offset
    unified <- c("sepsis_offset_sec","landmark_offset_sec",
                 "make30_risk_start_offset_sec","make30_censor_offset_sec",
                 "first_new_krt_offset_sec","final_inpatient_creatinine_offset_sec")
    if(!all(unified %in% names(d))) stop("eICU converted second-based fields are missing.")
    if(anyNA(d[unified[1:4]]) || any(shift<0) ||
       any(d$sepsis_offset_sec!=60*d$sepsis_offset) ||
       any(d$landmark_offset_sec!=60*d$landmark_offset) ||
       any(d$make30_risk_start_offset_sec!=d$landmark_offset_sec) ||
       any(d$make30_censor_offset_sec!=60*(d$make30_censor_hosp_offset-shift)) ||
       any(d$make30_censor_offset_sec<=d$make30_risk_start_offset_sec) ||
       any(d$make30_censor_offset_sec>30*86400))
      stop("eICU unified ICU-origin time axis is inconsistent.")
    for(pair in list(c("first_new_krt_hosp_offset","first_new_krt_offset_sec"),
                    c("final_inpatient_creatinine_hosp_offset","final_inpatient_creatinine_offset_sec"))) {
      old <- d[[pair[1L]]]; new <- d[[pair[2L]]]
      if(any(is.na(old)!=is.na(new)) ||
         any(new[!is.na(new)]!=60*(old[!is.na(old)]-shift[!is.na(old)])))
        stop("eICU event time conversion is inconsistent: ",pair[2L])
    }
    write_csv_once(data.frame(original_n=5867L,excluded_time_order_n=27L,
                             excluded_nonpositive_window_n=1L,retained_n=nrow(d),
                             unit="seconds",origin="index ICU admission",
                             source_SQL_conversion_verified=TRUE),"eICU_timebase_audit")
  }
  data.frame(dataset=db,export_n=nrow(d),outside_landmark_n=sum(d$make30_evaluable_window!=1L),
             make30_unknown_n=sum(is.na(d$make30)),
             analyzed_n=sum(good),events_n=sum(d$make30[good]==1L),
             event_pct=100*mean(d$make30[good]==1L))
}

quality <- list(); risks <- list(); tables <- list(); analyses <- list()
missing_outcomes <- list()
covars <- c("age","sex","race_ethnicity","weight_kg","bmi",
            "sepsis_sofa_score","baseline_creatinine","baseline_egfr_mdrd",
            "ckd_or_renal_disease","diabetes","congestive_heart_failure",
            "liver_disease","malignancy","mechanical_ventilation_0_24h",
            "vasopressor_0_24h")
binary_covars <- c("ckd_or_renal_disease","diabetes","congestive_heart_failure",
                   "liver_disease","malignancy","mechanical_ventilation_0_24h",
                   "vasopressor_0_24h")
core <- c("age","sex","sepsis_sofa_score","baseline_egfr_mdrd",
          "ckd_or_renal_disease","diabetes",
          "mechanical_ventilation_0_24h","vasopressor_0_24h")

robust_components <- function(d,db,adjusted=FALSE) {
  rhs <- c("phenotype",if(adjusted) core else character(0))
  f <- stats::as.formula(paste("make30 ~",paste(rhs,collapse=" + ")))
  mod <- stats::glm(f,data=d,family=stats::poisson(link="log"))
  if (anyNA(stats::coef(mod)) || nobs(mod)!=nrow(d))
    stop(db,": singular or incomplete Poisson model; inspect common covariates.")
  X <- stats::model.matrix(mod)
  U <- X*as.numeric(mod$y-stats::fitted(mod))
  invH <- stats::vcov(mod)
  sandwich <- invH%*%crossprod(U)%*%invH*nrow(X)/(nrow(X)-ncol(X))
  se <- sqrt(diag(sandwich))
  idx <- grep("^phenotypeC[23]$",names(stats::coef(mod)))
  if(length(idx)!=2L || any(!is.finite(se[idx])))
    stop(db,": robust phenotype contrast unavailable.")
  list(est=stats::coef(mod)[idx],se=se[idx],npar=ncol(X))
}
robust_rr <- function(d,db,adjusted=FALSE) {
  z <- robust_components(d,db,adjusted)
  est <- z$est
  data.frame(dataset=db,model=if(adjusted) "adjusted_common_covariates" else "crude",
             comparison=paste0(sub("^phenotype","",names(est))," vs C1"),
             RR=exp(est),CI95_low=exp(est-1.96*z$se),
             CI95_high=exp(est+1.96*z$se),
             robust_SE_logRR=z$se,p_value=2*stats::pnorm(-abs(est/z$se)),
             n=nrow(d),events=sum(d$make30),row.names=NULL)
}
robust_rr_mice <- function(d,db) {
  if (!requireNamespace("mice",quietly=TRUE))
    stop(db,": two sex records need MICE for adjusted analysis. ",
         "Install the laptop R package: install.packages('mice').")
  variables <- unique(c("sex",setdiff(core,"sex"),"phenotype","make30"))
  input <- d[variables]
  methods <- stats::setNames(rep("",length(variables)),variables)
  methods["sex"] <- "logreg"
  predictors <- matrix(0L,length(variables),length(variables),
                       dimnames=list(variables,variables))
  predictors["sex",setdiff(variables,"sex")] <- 1L
  # MAKE30 is a covariate-imputation predictor only; phenotype G and its
  # external posterior were fixed long before this analysis.
  imp <- mice::mice(input,m=20L,maxit=5L,method=methods,
                    predictorMatrix=predictors,seed=20260927L,printFlag=FALSE)
  na_index <- which(is.na(d$sex))
  fits <- vector("list",20L)
  record <- vector("list",20L)
  for (j in seq_len(20L)) {
    completed <- mice::complete(imp,action=j)
    z <- d
    z$sex <- factor(as.character(completed$sex),levels=c("Female","Male"))
    if(anyNA(z$sex)) stop(db,": incomplete MICE sex imputation.")
    fits[[j]] <- robust_components(z,db,adjusted=TRUE)
    record[[j]] <- data.frame(imputation=j,unit_id=as.character(
      z[[id_col[db]]][na_index]),reported_sex=z$sex_reported[na_index],
      imputed_binary_sex=as.character(z$sex[na_index]),row.names=NULL)
  }
  write_csv_once(do.call(rbind,record),paste0(db,"_sex_MICE_20x"))
  do.call(rbind,lapply(seq_len(2L),function(k) {
    q <- vapply(fits,function(z) unname(z$est[k]),numeric(1))
    u <- vapply(fits,function(z) unname(z$se[k]^2),numeric(1))
    pooled <- mice::pool.scalar(Q=q,U=u,n=nrow(d),k=fits[[1]]$npar,
                                rule="rubin1987")
    SE <- sqrt(pooled$t)
    crit <- stats::qt(0.975,df=pooled$df)
    data.frame(dataset=db,model="adjusted_common_covariates_MICE_20",
               comparison=paste0("C",k+1L," vs C1"),
               RR=exp(pooled$qbar),CI95_low=exp(pooled$qbar-crit*SE),
               CI95_high=exp(pooled$qbar+crit*SE),
               robust_SE_logRR=SE,
               p_value=2*stats::pt(-abs(pooled$qbar/SE),df=pooled$df),
               n=nrow(d),events=sum(d$make30),row.names=NULL)
  }))
}

for(db in names(cohort)) {
  d <- cohort[[db]]
  risks[[db]] <- check_outcome(d,db)
  s <- scored[[db]]
  d$assigned_class <- s$assigned_class[match(d[[id_col[db]]],s$unit_id)]
  d$max_posterior <- apply(as.matrix(s[match(d[[id_col[db]]],s$unit_id),
                                       paste0("p_class",1:G),drop=FALSE]),1L,max)
  if(anyNA(d$assigned_class) || !setequal(d$assigned_class,1:G))
    stop(db,": invalid final G3 class assignment.")
  pp <- as.matrix(s[paste0("p_class",1:G)])
  missing_outcomes[[db]] <- do.call(rbind,lapply(1:G,function(k) {
    z <- d[d$assigned_class==k,,drop=FALSE]
    data.frame(dataset=db,class=k,total_n=nrow(z),
               make30_unknown_n=sum(is.na(z$make30)),
               outside_landmark_n=sum(z$make30_evaluable_window!=1L),
               outcome_analysis_n=sum(z$make30_evaluable_window==1L &
                                      !is.na(z$make30)))
  }))
  ent <- 1+sum(pp[pp>0]*log(pp[pp>0]))/(nrow(pp)*log(G))
  quality[[db]] <- do.call(rbind,lapply(1:G,function(k) {
    mask <- s$assigned_class==k
    appa <- mean(pp[mask,k]); share <- mean(mask)
    data.frame(dataset=db,class=k,n=sum(mask),pct=100*share,
               APPA=appa,OCC=if(appa<1 && share<1)
                 (appa/(1-appa))/(share/(1-share)) else NA_real_,
               mean_max_posterior=mean(apply(pp[mask,,drop=FALSE],1,max)),
               pct_posterior_ge_090=100*mean(apply(pp[mask,,drop=FALSE],1,max)>=0.9),
               entropy=ent)
  }))
  keep <- d$make30_evaluable_window==1L & !is.na(d$make30)
  a <- d[keep,,drop=FALSE]
  a$phenotype <- factor(paste0("C",a$assigned_class),levels=paste0("C",1:G))
  a$sex_reported <- as.character(a$sex)
  a$sex <- tolower(trimws(as.character(a$sex)))
  a$sex <- ifelse(a$sex %in% c("f","female"),"Female",
                  ifelse(a$sex %in% c("m","male"),"Male",NA_character_))
  a$sex <- factor(a$sex,levels=c("Female","Male"))
  for(v in setdiff(core,"sex")) a[[v]] <- as.numeric(a[[v]])
  if(anyNA(a[setdiff(core,"sex")]))
    stop(db,": unexpectedly missing a frozen common adjustment variable.")
  if(db!="eICU" && anyNA(a$sex)) stop(db,": unexpectedly unknown sex.")
  if(db=="eICU" && sum(is.na(a$sex))!=2L)
    stop(db,": eICU binary sex unknown/other count changed from the audit.")
  analyses[[db]] <- a
  write_csv_once(data.frame(dataset=db,unit_id=a[[id_col[db]]],
                            phenotype=a$phenotype,posterior_max=a$max_posterior,
                            make30=a$make30),paste0(db,"analysis_cohort"))
  tables[[db]] <- do.call(rbind,lapply(1:G,function(k) {
    sub <- a[a$assigned_class==k,,drop=FALSE]
    do.call(rbind,lapply(covars,function(v) {
      x <- if(v=="sex") sub$sex_reported else sub[[v]]
      all_missing <- all(is.na(x))
      data.frame(dataset=db,phenotype=paste0("C",k),variable=v,
                 total_n=nrow(sub),missing_n=sum(is.na(x)),
                 status=if(all_missing) "unavailable" else "observed",
                 summary=if(all_missing) "not available"
                 else if(v%in%binary_covars)
                   sprintf("%d/%d (%.1f%%)",sum(as.numeric(x)==1L,na.rm=TRUE),
                           sum(!is.na(x)),100*mean(as.numeric(x)==1L,na.rm=TRUE))
                 else if(is.numeric(x))
                   paste0(format(round(stats::median(x,na.rm=TRUE),2),trim=TRUE),
                          " [",format(round(stats::quantile(x,0.25,na.rm=TRUE),2),trim=TRUE),
                          ", ",format(round(stats::quantile(x,0.75,na.rm=TRUE),2),trim=TRUE),"]")
                 else paste(names(sort(table(x),decreasing=TRUE)),
                            as.integer(sort(table(x),decreasing=TRUE)),
                            collapse="; "),row.names=NULL)
    }))
  }))
}
write_csv_once(do.call(rbind,risks),"landmark_outcome_flow")
write_csv_once(do.call(rbind,missing_outcomes),"MAKE30_missing_by_class")
write_csv_once(do.call(rbind,quality),"posterior_quality")
write_csv_once(do.call(rbind,tables),"table1_by_dataset_class")
event_rates <- do.call(rbind,lapply(names(analyses),function(db) {
  d <- analyses[[db]]
  do.call(rbind,lapply(1:G,function(k) {
    z <- d[d$assigned_class==k,,drop=FALSE]
    data.frame(dataset=db,class=k,n=nrow(z),make30_events=sum(z$make30),
               make30_event_pct=100*mean(z$make30))
  }))
}))
write_csv_once(event_rates,"MAKE30_rates_after_landmark")
associations <- do.call(rbind,lapply(names(analyses),function(db) {
  rbind(robust_rr(analyses[[db]],db,FALSE),
        if(db=="eICU") {
          cc <- robust_rr(analyses[[db]][!is.na(analyses[[db]]$sex),,drop=FALSE],
                          db,TRUE)
          cc$model <- "adjusted_binary_sex_complete_case_sensitivity"
          rbind(robust_rr_mice(analyses[[db]],db),cc)
        }
        else robust_rr(analyses[[db]],db,TRUE))
}))
write_csv_once(associations,"MAKE30_modified_Poisson_robust_RR")

notes <- c(
  "V316 Step 5 uses frozen G3 S01; V305/V306 Step 4B runs separately.",
  paste("Input dir:",input_dir),paste("Output dir:",out_dir),
  paste("gbmt:",as.character(utils::packageVersion("gbmt"))),
  paste("MIMIC final-parameter replay agreement:",training_agreement),
  paste("MIMIC final-parameter replay ARI:",training_ari),
  "Na/K/Cl/HCO3 t1..t4 align to MSI rows 2,5,8,11 (~3,9,15,21 h).",
  "V295 EM fitting retains source NA and uses internal conditional missing-value calculations.",
  "External posteriors use exactly V295 final pruned beta/Sigma/prior and gbmt 0.1.4 score.",
  "MIMIC analyses use V295 saved development cohort assignments.",
  "Only patients with evaluable 24 h landmark and known MAKE30 enter outcome analysis.",
  "SICdb weight/BMI/race/CHF/liver disease/malignancy have no observed records.",
  "The common adjusted model includes age, sex, SOFA, baseline eGFR, CKD, diabetes,",
  "mechanical ventilation and vasopressors. Two eICU recorded sex values (Unknown/Other)",
  "require MICE (20 logistic-chain imputations) for the binary adjusted sex covariate.",
  "No outcome or clinical prognosis was used to choose G or assign a phenotype.",
  "eICU landmark_offset uses ICU minutes; make30_risk_start_hosp_offset uses hospital minutes.",
  "eICU source SQL conversion was inspected; unified second-based ICU-origin fields are included.",
  "Outcome risk begins at sepsis onset +24h and ends at ICU admission +30 days or hospital discharge, whichever is earlier.",
  "eICU time-order exclusions: 27 reverse-origin records and one nonpositive observation window; original model and outcomes were not re-estimated."
)
write_text_once(notes,"run_notes")
message("Step 5 finished. Results saved in: ",out_dir)
