# V316 STEP5: MIMIC G3 S01 fixed-parameter replay ONLY.
# Derived from the inspected V307 scoring functions; no GBMTM refit.
# Run in RStudio: source(file.choose(), encoding="UTF-8")
# Reads the exact unpacked V307 kit under the current Windows profile.
# Or: Rscript V316_STEP5_replay_G3.R <unpacked_V307_input_directory>
# Requires exactly gbmt 0.1.4. No outcome models and no external transport.
# Five frozen input files are checked by exact names, sizes/MD5, IDs and settings.
# All outputs are new; official V295 files are never modified.
# Audit folder and ZIP are saved under GBMTM_stage5_output.

main <- function() {
old_options <- options(stringsAsFactors=FALSE, scipen=999, error=NULL)
on.exit(options(old_options), add=TRUE)
VERSION <- "V316"
EXPECTED <- c(
  "MIMIC_FINAL_MAKE30_ELECTROLYTEQC.csv"="f1eb2d7fac1303a682b53c679dde9abd",
  "V295_compact_fit_G3_S01.rds"="cf5881479ee206195aa80d59fa3bd17f",
  "V295_patient_class_G3_S01.csv"="dd12a3730af44c1035c79f4558c14db9",
  "V295_trajectory_G3_S01.csv"="3d47aa52a546266e937f50a26de2676a",
  "V295_task_summary_G3_S01.csv"="6b8eb067bb4013a54f6cb4dafe6d5f83"
)
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
run_name <- paste0(VERSION,"_STEP5_replay_G3_",format(Sys.time(),"%Y%m%d_%H%M%S"),"_",Sys.getpid())
out_dir <- file.path(output_root, run_name)
if (dir.exists(out_dir)) stop("Output already exists; rerun to create a new audit.")
if (!dir.create(out_dir, recursive=FALSE)) stop("Cannot create audit output: ", out_dir)
log_path <- file.path(out_dir, paste0(VERSION,"_replay_log_G3.txt"))
progress_log <- function(...) {
  line <- paste0(format(Sys.time(),"%Y-%m-%d %H:%M:%S")," | ",paste0(...,collapse=""))
  cat(line,"\n",sep="")
  cat(line,"\n",file=log_path,append=TRUE,sep="")
}
save_csv <- function(d, content) {
  target <- file.path(out_dir, paste0(VERSION,"_",content,"_G3.csv"))
  utils::write.csv(d,target,row.names=FALSE,na="",fileEncoding="UTF-8")
  invisible(target)
}
save_text <- function(lines, content) {
  target <- file.path(out_dir, paste0(VERSION,"_",content,"_G3.txt"))
  writeLines(enc2utf8(lines),target,useBytes=TRUE)
  invisible(target)
}
progress_log("Replay started; no refitting, outcome analysis or external scoring.")
progress_log("Input directory: ",input_dir)
status <- tryCatch({
files <- file.path(input_dir,names(EXPECTED))
exists <- file.exists(files)
actual <- rep(NA_character_,length(files)); names(actual) <- names(EXPECTED)
if (any(exists)) actual[exists] <- unname(tools::md5sum(files[exists]))
save_csv(data.frame(filename=names(EXPECTED),exists=exists,
                   expected_md5=unname(EXPECTED),actual_md5=unname(actual),
                   matches=!is.na(actual)&actual==EXPECTED),"input_audit")
if (any(!exists)) stop("Missing exact input files: ",paste(names(EXPECTED)[!exists],collapse=", "))
if (anyNA(actual) || any(actual!=EXPECTED)) stop("Frozen input MD5 mismatch; no inference or substitution allowed.")
if (!requireNamespace("gbmt",quietly=TRUE)) stop("gbmt 0.1.4 is required; install.packages('gbmt').")
if (as.character(utils::packageVersion("gbmt"))!="0.1.4") stop("Requires exactly gbmt 0.1.4, matching V295.")
progress_log("All five input hashes match. gbmt version: 0.1.4.")
Z <- c("msi_z","sodium_z","potassium_z","chloride_z","bicarbonate_z")
FROZEN_MEAN <- setNames(c(1.229422,138.931869,4.221718,
                          105.125243,22.905997),Z)
FROZEN_SD <- setNames(c(0.3219380,5.3785897,0.6426219,
                        6.5980074,4.7412984),Z)
G <- 3L


read_cohort <- function(filename,id_col) {
  utils::read.csv(file.path(input_dir,filename),check.names=FALSE,
                 colClasses=setNames("character",id_col),
                 na.strings=c("","NA","NULL"),stringsAsFactors=FALSE)
}
cohort <- list(MIMIC=read_cohort("MIMIC_FINAL_MAKE30_ELECTROLYTEQC.csv","stay_id"))
id_col <- c(MIMIC="stay_id")
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
  stop("Formal V295 G3 S01 model identity/settings do not match.")
if (nrow(ref)!=9804L || anyNA(ref$unit_id) || anyDuplicated(ref$unit_id) ||
    anyNA(ref$assigned_class) || !setequal(as.integer(ref$assigned_class),seq_len(G)) ||
    !all(paste0("p_class",1:G)%in%names(ref))) stop("Formal patient classification is invalid.")
d <- cohort$MIMIC
raw_cols <- c(paste0("msi_m",1:12),
              unlist(lapply(c("sodium","potassium","chloride","bicarbonate"),
                             function(v) paste0(v,"_t",1:4))))
if(nrow(d)!=9804L || anyNA(d$stay_id) || anyDuplicated(d$stay_id) ||
   !all(raw_cols%in%names(d)) || !setequal(d$stay_id,ref$unit_id))
  stop("MIMIC count/IDs/observed trajectory fields do not match the frozen training cohort.")
times <- as.numeric(rownames(as.matrix(fit$Xmat)))
hours <- seq(1,23,by=2)
if(length(times)!=12L || anyNA(times) || anyDuplicated(times) ||
   is.unsorted(times,strictly=TRUE) ||
   max(abs(as.matrix(fit$Xmat)-cbind(1,(times-min(times)+1),
                                (times-min(times)+1)^2)))>1e-7 ||
   !all(vapply(hours,function(h) any(abs(traj$time_hour_center-h)<1e-8),logical(1))) ||
   max(abs(sort(unique(traj$time_model))-times))>1e-9)
  stop("Frozen polynomial design and true time anchors differ.")
standardization <- vector("list",length(Z))
for (j in seq_along(Z)) {
  prefix <- c("msi_m","sodium_t","potassium_t","chloride_t","bicarbonate_t")[j]
  count <- if (j==1L) 12L else 4L
  observed <- unlist(d[paste0(prefix,1:count)],use.names=FALSE)
  observed <- as.numeric(observed[!is.na(observed)])
  standardization[[j]] <- data.frame(variable=Z[j],n_observed=length(observed),
      mean_raw=mean(observed),sd_raw=stats::sd(observed),
      frozen_mean=FROZEN_MEAN[j],frozen_sd=FROZEN_SD[j],
      mean_abs_diff=abs(mean(observed)-FROZEN_MEAN[j]),
      sd_abs_diff=abs(stats::sd(observed)-FROZEN_SD[j]))
}
standardization <- do.call(rbind,standardization)
save_csv(standardization,"standardization_audit")
if(any(standardization$mean_abs_diff>0.00001 | standardization$sd_abs_diff>0.00001))
  stop("Frozen MIMIC standardization does not match the raw export.")
# Retain only IDs and original observed trajectory columns; no outcome is used.
cohort$MIMIC <- d[c("stay_id",raw_cols)]
rm(d)
progress_log("9804 patient IDs, polynomial design, and raw standardization verified.")
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


panel <- make_panel(cohort$MIMIC,"MIMIC")
missing_before <- colSums(is.na(panel[Z]))
mimic_score <- transport(panel,"MIMIC")
if(!identical(missing_before,colSums(is.na(panel[Z])))) stop("Original panel missingness changed.")
save_csv(data.frame(variable=Z,n_missing_original=as.integer(missing_before)),"panel_missingness")
save_csv(mimic_score,"MIMIC_final_parameter_replay")
if(nrow(mimic_score)!=nrow(ref) || !setequal(mimic_score$unit_id,ref$unit_id))
  stop("Replay output patient IDs differ from the formal reference.")
ref_ord <- ref[match(mimic_score$unit_id,ref$unit_id),,drop=FALSE]
training_agreement <- mean(mimic_score$assigned_class==ref_ord$assigned_class)
training_ari <- ari(mimic_score$assigned_class,ref_ord$assigned_class)
ref_p <- as.matrix(ref_ord[paste0("p_class",1:G)])
new_p <- as.matrix(mimic_score[paste0("p_class",1:G)])
if(anyNA(ref_p) || any(!is.finite(ref_p))) stop("Formal posterior values are invalid.")
diff_p <- abs(ref_p-new_p)
legacy_preflight <- is.finite(training_ari) && training_agreement>=0.95 && training_ari>=0.90
same_class <- mimic_score$assigned_class==ref_ord$assigned_class
audit <- data.frame(n_patients=nrow(ref),n_class_changed=sum(!same_class),
     saved_v295_class_agreement=training_agreement,saved_v295_ARI=training_ari,
     posterior_mean_abs_diff=mean(diff_p),posterior_max_abs_diff=max(diff_p),
     posterior_mean_abs_diff_after_rounding_4dp=mean(abs(ref_p-round(new_p,4))),
     posterior_all_identical_after_rounding_4dp=all(abs(ref_p-round(new_p,4))<1e-12),
     exact_class_reproduction=all(same_class),
     V307_original_preflight_pass=legacy_preflight,
     model_input_md5=task$input_md5,
     fit_md5=actual["V295_compact_fit_G3_S01.rds"],gbmt_version="0.1.4")
save_csv(audit,"MIMIC_training_replay_audit")
comparison <- data.frame(unit_id=mimic_score$unit_id,saved_class=ref_ord$assigned_class,
                         replay_class=mimic_score$assigned_class,class_changed=!same_class,
                         posterior_max_abs_diff=apply(diff_p,1L,max))
for(k in seq_len(G)) {
  comparison[[paste0("saved_p_class",k)]] <- ref_p[,k]
  comparison[[paste0("replay_p_class",k)]] <- new_p[,k]
}
save_csv(comparison,"MIMIC_patient_comparison")
save_csv(comparison[!same_class,,drop=FALSE],"MIMIC_changed_patients")
transition <- as.data.frame(table(saved_class=factor(ref_ord$assigned_class,levels=1:G),
                                 final_parameter_class=factor(mimic_score$assigned_class,levels=1:G)))
save_csv(transition,"MIMIC_saved_x_final_parameter_transition")
class_audit <- do.call(rbind,lapply(seq_len(G),function(k) {
  old <- ref_ord$assigned_class==k; new <- mimic_score$assigned_class==k
  data.frame(class=k,n_saved=sum(old),n_replay=sum(new),n_overlap=sum(old&new),
             recall=sum(old&new)/sum(old),jaccard=sum(old&new)/sum(old|new),
             proportion_drift_pp=100*(mean(new)-mean(old)))
}))
save_csv(class_audit,"MIMIC_class_replay_audit")
review_status <- if(all(same_class)) "EXACT_CLASS_REPRODUCTION" else "CLASS_DIFFERENCES_REQUIRE_REVIEW"
save_text(c(
  paste0("Status: ",review_status),
  paste0("Patients: ",nrow(ref)),
  paste0("Changed patient classes: ",sum(!same_class)),
  paste0("Class agreement: ",format(training_agreement,digits=12)),
  paste0("ARI: ",format(training_ari,digits=12)),
  paste0("Posterior mean absolute difference: ",format(mean(diff_p),digits=12)),
  paste0("Posterior maximum absolute difference: ",format(max(diff_p),digits=12)),
  paste0("Identical probabilities after 4dp rounding: ",audit$posterior_all_identical_after_rounding_4dp),
  paste0("V307 original preflight (agreement >=0.95 and ARI >=0.90): ",legacy_preflight),
  "The V307 preflight is an existing stop rule, NOT proof of exact reproduction.",
  "V295 saves pre-pruning posterior/assignments but post-pruning beta/Sigma.",
  "Only final frozen parameters are used here. No parameter is estimated or replaced.",
  "Reference class numbers are kept unchanged; no relabeling to conceal discrepancies.",
  "All results require review before proceeding to outcome analysis.",
  "No MAKE30 model, external classification, or GBMTM refitting was performed."
),"replay_summary")
progress_log("Agreement=",format(training_agreement,digits=8),"; ARI=",format(training_ari,digits=8),
    "; changed classes=",sum(!same_class),"; posterior MAE=",format(mean(diff_p),digits=8))
list(ok=TRUE,status=review_status)
}, error=function(e) {
  progress_log("ERROR: ",conditionMessage(e))
  save_text(c("Status: REPLAY_ERROR",conditionMessage(e)),"replay_summary")
  list(ok=FALSE,error=conditionMessage(e))
})
save_text(capture.output(utils::sessionInfo()),"session_info")
zip_path <- file.path(output_root,paste0(run_name,".zip"))
ps_quote <- function(x) paste0("'",gsub("'","''",x,fixed=TRUE),"'")
zip_error <- NULL
if(.Platform$OS.type=="windows") {
  ps <- Sys.which("powershell.exe")
  if(!nzchar(ps)) {
    zip_error <- "powershell.exe unavailable; the output folder is complete but ZIP was not created."
  } else {
    ps_file <- tempfile("V316_zip_",fileext=".ps1")
    on.exit(unlink(ps_file),add=TRUE)
    writeLines(c("$ErrorActionPreference = 'Stop'",
      paste0("$src = ",ps_quote(normalizePath(out_dir,winslash="/",mustWork=TRUE))),
      paste0("$dst = ",ps_quote(zip_path)),
      "$items = @(Get-ChildItem -LiteralPath $src -File | Select-Object -ExpandProperty FullName)",
      "Compress-Archive -LiteralPath $items -DestinationPath $dst -CompressionLevel Optimal"),
      ps_file,useBytes=TRUE)
    code <- tryCatch(system2(ps,c("-NoProfile","-NonInteractive","-ExecutionPolicy","Bypass",
                                "-File",shQuote(ps_file))),error=function(e) 1L)
    if(!identical(as.integer(code),0L) || !file.exists(zip_path))
      zip_error <- "ZIP creation failed; the output folder is complete and can be compressed manually."
  }
} else {
  zip_bin <- Sys.which("zip")
  if(!nzchar(zip_bin)) {
    zip_error <- "zip executable unavailable; compress the output folder manually."
  } else {
    oldwd <- getwd()
    on.exit(setwd(oldwd),add=TRUE)
    setwd(out_dir)
    code <- tryCatch(utils::zip(zip_path,list.files(out_dir,full.names=FALSE),zip=zip_bin),
                     error=function(e) 1L)
    setwd(oldwd)
    if(!file.exists(zip_path)) zip_error <- "ZIP creation failed; compress the output folder manually."
  }
}
cat("\nOutput folder: ",out_dir,"\n",sep="")
if(is.null(zip_error)) cat("Send this ZIP for review: ",zip_path,"\n",sep="") else warning(zip_error,call.=FALSE)
if(!isTRUE(status$ok)) {
  message("Replay failed; the error was recorded. Send the audit ZIP for review.")
  return(invisible(list(status=status,output_directory=out_dir,
                       zip=if(file.exists(zip_path)) zip_path else NULL)))
}
cat("Replay completed; review the results before any outcome analysis.\n")
invisible(list(status=status,output_directory=out_dir,zip=if(file.exists(zip_path)) zip_path else NULL))
}
main()
