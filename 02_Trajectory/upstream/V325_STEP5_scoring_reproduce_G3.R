# V325: paths only; frozen V318 classifier retained without refitting.
# V318: recover the original V295 G3 S01 classification rule for fixed scoring.
# This executable only reads the delivered frozen scoring state. No EM fitting,
# pruning, G selection, outcome models or replacement of any old file occurs.
# Inputs are read from the archive at Z:/Research/MSI-1.
# RStudio: open this complete file, Ctrl+A, Ctrl+Enter.
# CLI: Rscript V325_STEP5_scoring_reproduce_G3.R

main_alignment <- function() {
  old_options <- options(stringsAsFactors=FALSE,scipen=999,error=NULL)
  on.exit(options(old_options),add=TRUE)
  VERSION <- "V318"
  Z <- c("msi_z","sodium_z","potassium_z","chloride_z","bicarbonate_z")
  G <- 3L
  EXPECTED_STATE <- "89f684d00680fafc30577d677fee6a30"
  EXPECTED <- c(
    "MIMIC_FINAL_MAKE30_ELECTROLYTEQC.csv"="f1eb2d7fac1303a682b53c679dde9abd",
    "SICdb_FINAL_MAKE30_ELECTROLYTEQC.csv"="d52a3ecdbb636bba2859eb9f102fdf45",
    "eICU_FINAL_RAW_ELECTROLYTEQC.csv"="a06fddbcbc3b93f2c83227201c1138cd",
    "V295_compact_fit_G3_S01.rds"="cf5881479ee206195aa80d59fa3bd17f",
    "V295_patient_class_G3_S01.csv"="dd12a3730af44c1035c79f4558c14db9",
    "V295_trajectory_G3_S01.csv"="3d47aa52a546266e937f50a26de2676a",
    "V295_task_summary_G3_S01.csv"="6b8eb067bb4013a54f6cb4dafe6d5f83"
  )
  archive_root <- "Z:/Research/MSI-1"
  output_root <- file.path(archive_root,"STEP5_冻结分型与预后外部验证","独立复现")
  dir.create(output_root,recursive=TRUE,showWarnings=FALSE)
  if(!dir.exists(output_root)) stop("Exact STEP5 root does not exist: ",output_root)
  args <- character(0L)
  input_dir <- file.path(archive_root,"STEP5_冻结分型与预后外部验证","冻结分类")
  state_path <- file.path(archive_root,"RDS","V318_frozen_scoring_state_G3.rds")
  input_dir <- normalizePath(input_dir,winslash="/",mustWork=TRUE)
  if(!file.exists(state_path)) stop("Exact delivered scoring state is missing: ",state_path)
  state_path <- normalizePath(state_path,winslash="/",mustWork=TRUE)
  if(!identical(unname(tools::md5sum(state_path)),EXPECTED_STATE))
    stop("Delivered V318 scoring state MD5 differs; no substitution allowed.")
  input_file <- function(name) {
    if(name %in% c("MIMIC_FINAL_MAKE30_ELECTROLYTEQC.csv","eICU_FINAL_RAW_ELECTROLYTEQC.csv","SICdb_FINAL_MAKE30_ELECTROLYTEQC.csv")) file.path(archive_root,name)
    else if(identical(name,"V295_compact_fit_G3_S01.rds")) file.path(archive_root,"RDS",name)
    else if(name %in% c("V295_patient_class_G3_S01.csv","V295_trajectory_G3_S01.csv","V295_task_summary_G3_S01.csv")) file.path(archive_root,"STEP4_多起点精修与稳定性验证","多起点结果",name)
    else file.path(input_dir,name)
  }
  input_paths <- vapply(names(EXPECTED),input_file,character(1L))
  if(any(!file.exists(input_paths))) stop("Missing exact input files: ",paste(names(EXPECTED)[!file.exists(input_paths)],collapse=", "))
  actual <- unname(tools::md5sum(input_paths))
  if(!identical(actual,unname(EXPECTED)))
    stop("Frozen input identity differs: ",paste(names(EXPECTED)[is.na(actual)|actual!=EXPECTED],collapse=", "))
  if(!requireNamespace("gbmt",quietly=TRUE) || as.character(utils::packageVersion("gbmt"))!="0.1.4")
    stop("Requires exactly gbmt 0.1.4, matching the original V295 computation.")
  state <- readRDS(state_path)
  required_state <- c("schema","G","degree","variables","source_long_input_md5","source_fit_md5",
    "source_patient_class_md5","source_task_md5","source_MIMIC_md5","original_seed",
    "completion_parameters_iteration","score_parameters_iteration","original_stop_reason","gbmt_version",
    "raw_center","raw_scale","time_model","Xmat","prior","completion_beta","completion_Sigma",
    "score_beta","score_Sigma","archived_class_agreement","archived_posterior_4dp_max_abs_diff")
  if(!all(required_state%in%names(state)) || state$schema!="V318_recovered_V295_G3_S01" ||
     state$G!=3L || state$degree!=2L || !identical(state$variables,Z) ||
     state$source_long_input_md5!="0361d8b59184b7a6b0a5fa9c946088a0" ||
     state$source_fit_md5!=EXPECTED["V295_compact_fit_G3_S01.rds"] ||
     state$source_patient_class_md5!=EXPECTED["V295_patient_class_G3_S01.csv"] ||
     state$source_task_md5!=EXPECTED["V295_task_summary_G3_S01.csv"] ||
     state$source_MIMIC_md5!=EXPECTED["MIMIC_FINAL_MAKE30_ELECTROLYTEQC.csv"] ||
     state$completion_parameters_iteration!=162L || state$score_parameters_iteration!=163L ||
     state$gbmt_version!="0.1.4" || state$archived_class_agreement!=1 ||
     state$archived_posterior_4dp_max_abs_diff!=0)
    stop("Delivered frozen-state identity/settings differ from the verified G3 S01 rule.")
  times <- state$time_model
  if(length(times)!=12L || anyNA(times) || anyDuplicated(times) || is.unsorted(times,strictly=TRUE) ||
     !identical(names(state$raw_center),Z) || !identical(names(state$raw_scale),Z) ||
     any(!is.finite(state$raw_center)) || any(!is.finite(state$raw_scale)|state$raw_scale<=0) ||
     !identical(dim(state$Xmat),c(12L,3L)) ||
     max(abs(state$Xmat-cbind(1,times-min(times)+1,(times-min(times)+1)^2)))>1e-12 ||
     length(state$prior)!=3L || any(!is.finite(state$prior)|state$prior<=0) ||
     abs(sum(state$prior)-1)>1e-10)
    stop("Frozen time design, standardization or class prior is invalid.")
  for(k in seq_len(G)) {
    for(field in c("completion_beta","score_beta")) {
      b <- state[[field]][[k]]
      if(!identical(dim(b),c(3L,5L)) || !identical(colnames(b),Z) || any(!is.finite(b)))
        stop("Invalid frozen beta: ",field," class ",k)
    }
    for(field in c("completion_Sigma","score_Sigma")) {
      S <- state[[field]][[k]]
      if(!identical(dim(S),c(5L,5L)) || !identical(rownames(S),Z) ||
         !identical(colnames(S),Z) || any(!is.finite(S)) || det(S)<=0)
        stop("Invalid frozen covariance: ",field," class ",k)
    }
  }
  key <- paste0(substr(EXPECTED_STATE,1L,10L),"_",substr(EXPECTED["MIMIC_FINAL_MAKE30_ELECTROLYTEQC.csv"],1L,8L))
  out_dir <- file.path(output_root,paste0(VERSION,"_STEP5_scoring_G3_",key))
  dir.create(out_dir,showWarnings=FALSE)
  if(!dir.exists(out_dir)) stop("Cannot create output directory: ",out_dir)
  log_file <- file.path(out_dir,paste0(VERSION,"_scoring_log_G3.txt"))
  progress_log <- function(...) {
    line <- paste0(format(Sys.time(),"%Y-%m-%d %H:%M:%S")," | ",paste0(...,collapse=""))
    cat(line,"\n",sep=""); cat(line,"\n",file=log_file,append=TRUE,sep="")
  }
  write_csv <- function(d,content) {
    dest <- file.path(out_dir,paste0(VERSION,"_",content,"_G3.csv"))
    tmp <- paste0(dest,".tmp_",Sys.getpid())
    utils::write.csv(d,tmp,row.names=FALSE,na="",fileEncoding="UTF-8")
    if(file.exists(dest)) unlink(dest)
    if(!file.rename(tmp,dest)) stop("Cannot save generated report: ",dest)
    invisible(dest)
  }
  write_text <- function(lines,content) {
    dest <- file.path(out_dir,paste0(VERSION,"_",content,"_G3.txt"))
    writeLines(enc2utf8(lines),dest,useBytes=TRUE)
    invisible(dest)
  }
  progress_log("Frozen G3 S01 scoring started; no fitting, pruning or outcome analysis.")
  status <- tryCatch({
    write_csv(data.frame(filename=c(names(EXPECTED),basename(state_path)),
      expected_md5=c(unname(EXPECTED),EXPECTED_STATE),actual_md5=c(actual,EXPECTED_STATE),matches=TRUE),"input_audit")
    read_cohort <- function(filename,id) utils::read.csv(input_file(filename),
      check.names=FALSE,colClasses=setNames("character",id),na.strings=c("","NA","NULL"))
    datasets <- c("MIMIC","SICdb","eICU")
    filenames <- c(MIMIC="MIMIC_FINAL_MAKE30_ELECTROLYTEQC.csv",SICdb="SICdb_FINAL_MAKE30_ELECTROLYTEQC.csv",eICU="eICU_FINAL_RAW_ELECTROLYTEQC.csv")
    id_cols <- c(MIMIC="stay_id",SICdb="case_id",eICU="patientunitstayid")
    expected_n <- c(MIMIC=9804L,SICdb=3831L,eICU=5839L)
    ref <- read_cohort("V295_patient_class_G3_S01.csv","unit_id")
    archived <- readRDS(input_file("V295_compact_fit_G3_S01.rds"))
    task <- utils::read.csv(input_file("V295_task_summary_G3_S01.csv"),check.names=FALSE)
    if(nrow(task)!=1L || task$requested_g!=G || task$start_id!=1L || task$seed!=state$original_seed ||
       task$em_iterations!=163L || task$stop_reason!=state$original_stop_reason ||
       task$input_md5!=state$source_long_input_md5 || nrow(ref)!=9804L ||
       anyNA(ref$unit_id) || anyDuplicated(ref$unit_id) ||
       !all(c("assigned_class",paste0("p_class",1:G))%in%names(ref)) ||
       !inherits(archived,"gbmt") || !identical(archived$call$x.names,Z) ||
       max(abs(state$prior-archived$prior))>1e-12)
      stop("Formal G3 S01 archive/task identity differs.")
    arch_class <- as.integer(archived$assign[match(ref$unit_id,names(archived$assign))])
    arch_p <- archived$posterior[match(ref$unit_id,rownames(archived$posterior)),,drop=FALSE]
    if(anyNA(arch_class) || !all(arch_class==ref$assigned_class) ||
       max(abs(arch_p-as.matrix(ref[paste0("p_class",1:G)])))!=0)
      stop("Original archived assignments/posteriors differ from the formal CSV.")
    raw_cols <- c(paste0("msi_m",1:12),unlist(lapply(c("sodium","potassium","chloride","bicarbonate"),function(v)paste0(v,"_t",1:4))))
    make_panel <- function(d,db) {
      id <- id_cols[db]
      if(!all(c(id,raw_cols)%in%names(d)) || nrow(d)!=expected_n[db] ||
         anyNA(d[[id]]) || anyDuplicated(d[[id]])) stop(db,": exact cohort schema/count/IDs differ.")
      ids <- as.character(d[[id]])
      ii <- rep(seq_along(ids),each=12L); tt <- rep(1:12,times=length(ids))
      panel <- data.frame(unit_id=ids[ii],time_model=times[tt],
        msi_z=as.numeric(as.matrix(d[paste0("msi_m",1:12)])[cbind(ii,tt)]),
        sodium_z=NA_real_,potassium_z=NA_real_,chloride_z=NA_real_,bicarbonate_z=NA_real_)
      anchor <- match(tt,c(2L,5L,8L,11L))
      rr <- which(!is.na(anchor))
      for(j in 2:5) {
        prefix <- c("sodium","potassium","chloride","bicarbonate")[j-1L]
        panel[[Z[j]]][rr] <- as.numeric(as.matrix(d[paste0(prefix,"_t",1:4)])[cbind(ii[rr],anchor[rr])])
      }
      for(v in Z) panel[[v]] <- (panel[[v]]-state$raw_center[v])/state$raw_scale[v]
      panel <- panel[order(panel$unit_id,panel$time_model),,drop=FALSE]; rownames(panel)<-NULL
      if(anyDuplicated(panel[c("unit_id","time_model")]) || nrow(panel)!=12L*nrow(d) ||
         any(is.infinite(as.matrix(panel[Z])))) stop(db,": invalid observed panel.")
      counts <- stats::aggregate(!is.na(panel[Z]),list(unit_id=panel$unit_id),sum)
      if(any(as.matrix(counts[Z])==0L)) stop(db,": a patient has no observation for at least one marker.")
      panel
    }
    imputFun <- getFromNamespace("imputFun","gbmt")
    logdmN <- getFromNamespace("logdmN","gbmt")
    pldBound <- getFromNamespace("pldBound","gbmt")
    score_fixed <- function(panel,db) {
      ids <- unique(panel$unit_id); ix <- match(panel$time_model,times)
      stopifnot(!anyNA(ix))
      scores <- matrix(NA_real_,length(ids),G)
      observed <- as.matrix(panel[Z])
      original_mask <- is.na(observed)
      for(k in 1:G) {
        # Exactly reproduce the final original EM scoring state: its
        # conditional-computation parameters are from iteration162,
        # while its likelihood parameters are from iteration163.
        completed <- imputFun(observed,panel$time_model-min(times)+1,
          state$completion_beta[[k]],state$completion_Sigma[[k]],2L)
        if(anyNA(completed) || any(!is.finite(completed))) stop(db,": conditional calculation failed.")
        b <- state$score_beta[[k]]; S <- state$score_Sigma[[k]]
        mu <- state$Xmat%*%b
        diff <- completed-mu[ix,,drop=FALSE]
        density <- -ncol(S)/2*log(2*pi)-.5*log(det(S))-.5*rowSums((diff%*%solve(S))*diff)
        for(i in unique(c(1L,max(1L,length(ids)%/%2L),length(ids)))) {
          rows <- which(panel$unit_id==ids[i])
          check <- logdmN(completed[rows,,drop=FALSE],mu[ix[rows],,drop=FALSE],S)
          if(length(check)!=length(rows) || anyNA(check) || max(abs(check-density[rows]))>1e-8)
            stop(db,": density calculation differs from installed gbmt0.1.4.")
        }
        summed <- rowsum(matrix(density,ncol=1L),panel$unit_id,reorder=FALSE)
        stopifnot(identical(rownames(summed),ids))
        scores[,k] <- vapply(as.numeric(summed[,1L])+log(state$prior[k]),pldBound,numeric(1L),USE.NAMES=FALSE)
        progress_log(db," class ",k,"/3 scored")
      }
      stopifnot(identical(is.na(observed),original_mask),identical(is.na(as.matrix(panel[Z])),original_mask))
      prob <- exp(scores-apply(scores,1L,max)); prob <- prob/rowSums(prob)
      if(anyNA(prob) || any(!is.finite(prob)) || max(abs(rowSums(prob)-1))>1e-12)
        stop(db,": invalid frozen class posterior.")
      ans <- data.frame(unit_id=ids,assigned_class=max.col(prob,ties.method="first"),max_posterior=apply(prob,1L,max))
      for(k in 1:G) ans[[paste0("p_class",k)]] <- prob[,k]
      ans
    }
    score_or_resume <- function(db) {
      d <- read_cohort(filenames[db],id_cols[db])
      panel <- make_panel(d,db)
      dest <- file.path(out_dir,paste0(VERSION,"_",db,"_fixed_scoring_G3.csv"))
      if(file.exists(dest)) {
        scored <- utils::read.csv(dest,check.names=FALSE,colClasses=c(unit_id="character"))
        expected_columns <- c("unit_id","assigned_class","max_posterior",paste0("p_class",1:G))
        if(!identical(names(scored),expected_columns) || nrow(scored)!=nrow(d) || anyNA(scored) ||
           anyDuplicated(scored$unit_id) || !identical(scored$unit_id,unique(panel$unit_id)) ||
           any(!scored$assigned_class%in%1:G) ||
           any(as.matrix(scored[paste0("p_class",1:G)])<0|as.matrix(scored[paste0("p_class",1:G)])>1) ||
           max(abs(rowSums(scored[paste0("p_class",1:G)])-1))>1e-10)
          stop("Existing V318 scoring output is incompatible: ",dest)
        progress_log(db," existing complete score file reused")
      } else {
        scored <- score_fixed(panel,db)
        write_csv(scored,paste0(db,"_fixed_scoring"))
      }
      list(scored=scored,panel_missing=colSums(is.na(panel[Z])))
    }
    mimic <- score_or_resume("MIMIC")
    scored <- mimic$scored
    ref_order <- ref[match(scored$unit_id,ref$unit_id),,drop=FALSE]
    if(anyNA(ref_order$unit_id)) stop("Replayed MIMIC IDs differ from the formal classification.")
    ref_prob <- as.matrix(ref_order[paste0("p_class",1:G)])
    new_prob <- as.matrix(scored[paste0("p_class",1:G)])
    audit <- data.frame(n_patients=nrow(scored),n_class_changed=sum(scored$assigned_class!=ref_order$assigned_class),
      class_agreement=mean(scored$assigned_class==ref_order$assigned_class),
      posterior_max_abs_diff=max(abs(new_prob-ref_prob)),
      posterior_max_abs_diff_4dp=max(abs(round(new_prob,4)-ref_prob)),
      exact_classes=all(scored$assigned_class==ref_order$assigned_class),
      exact_saved_posterior_4dp=all(round(new_prob,4)==ref_prob),
      fit_md5=EXPECTED["V295_compact_fit_G3_S01.rds"],scoring_state_md5=EXPECTED_STATE,gbmt_version="0.1.4")
    write_csv(audit,"MIMIC_replay_audit")
    if(!isTRUE(audit$exact_classes) || !isTRUE(audit$exact_saved_posterior_4dp))
      stop("Original G3 classification rule is not exactly reproduced; external scoring was not started. See audit.")
    progress_log("MIMIC PASS: 9804/9804 original classes; original four-decimal posterior reproduced.")
    # Preserve official MIMIC classes/probabilities as the analysis input.
    write_csv(ref[,c("unit_id","assigned_class",paste0("p_class",1:G))],"MIMIC_formal_classification")
    all_scored <- list(MIMIC=mimic)
    for(db in c("SICdb","eICU")) all_scored[[db]] <- score_or_resume(db)
    quality <- do.call(rbind,lapply(datasets,function(db) {
      s <- all_scored[[db]]$scored
      do.call(rbind,lapply(1:G,function(k) {
        q <- s$assigned_class==k
        data.frame(dataset=db,class=k,n=sum(q),pct=100*mean(q),
          APPA=if(any(q))mean(s[[paste0("p_class",k)]][q]) else NA_real_,
          median_max_posterior=if(any(q))median(s$max_posterior[q]) else NA_real_,
          pct_max_posterior_ge_070=if(any(q))100*mean(s$max_posterior[q]>=.7) else NA_real_,
          pct_max_posterior_ge_090=if(any(q))100*mean(s$max_posterior[q]>=.9) else NA_real_)
      }))
    }))
    write_csv(quality,"classification_quality")
    write_csv(do.call(rbind,lapply(datasets,function(db)data.frame(dataset=db,variable=Z,n_missing_original=all_scored[[db]]$panel_missing))),"observed_panel_missingness")
    write_csv(data.frame(variable=Z,mean_original_precision=state$raw_center,sd_original_precision=state$raw_scale),"frozen_standardization")
    if(!identical(unname(tools::md5sum(input_paths)),unname(EXPECTED)) ||
       !identical(unname(tools::md5sum(state_path)),EXPECTED_STATE)) stop("An original input/state file changed during scoring.")
    write_text(c("Status: EXACT_ORIGINAL_CLASSIFICATION_REPRODUCED",
      "MIMIC original classes: 9804/9804 (100%).",
      "MIMIC original posterior after4dp rounding: exact.",
      "SICdb and eICU: scored using exactly the same frozen original-classification rule.",
      "Original V295 model, trajectories, assignments and all cohort files remain unchanged.",
      "No EM fitting, G selection, pruning or outcome analysis was run by this executable.",
      "V318 uses the recovered original scoring state; V317 post-pruning replay is retained separately.",
      "This result establishes technical scoring consistency; external clinical validation still requires review."),"scoring_summary")
    write_text(capture.output(utils::sessionInfo()),"session_info")
    progress_log("Fixed scoring finished; all original input hashes remain unchanged.")
    list(ok=TRUE,audit=audit)
  },error=function(e) {
    write_text(c("Status: SCORING_ERROR",conditionMessage(e)),"scoring_summary")
    progress_log("ERROR: ",conditionMessage(e))
    list(ok=FALSE,error=conditionMessage(e))
  })
  zip_path <- paste0(out_dir,".zip")
  zip_ok <- FALSE
  if(.Platform$OS.type=="windows") {
    ps <- Sys.which("powershell.exe")
    if(nzchar(ps)) {
      ps_quote <- function(x)paste0("'",gsub("'","''",x,fixed=TRUE),"'")
      ps_file <- tempfile("V318_zip_",fileext=".ps1")
      on.exit(unlink(ps_file),add=TRUE)
      lines <- c("$ErrorActionPreference = 'Stop'",
        paste0("$src = ",ps_quote(normalizePath(out_dir,winslash="/",mustWork=TRUE))),
        paste0("$dst = ",ps_quote(zip_path)),
        "$items = @(Get-ChildItem -LiteralPath $src -File | Select-Object -ExpandProperty FullName)",
        "Compress-Archive -LiteralPath $items -DestinationPath $dst -CompressionLevel Optimal -Force")
      con <- file(ps_file,"wb")
      writeBin(c(as.raw(c(239,187,191)),charToRaw(enc2utf8(paste(lines,collapse="\r\n")))),con)
      close(con)
      code <- tryCatch(system2(ps,c("-NoProfile","-NonInteractive","-ExecutionPolicy","Bypass","-File",shQuote(ps_file))),error=function(e)1L)
      zip_ok <- identical(as.integer(code),0L)&&file.exists(zip_path)
    }
  } else if(nzchar(Sys.which("zip"))) {
    oldwd <- getwd(); on.exit(setwd(oldwd),add=TRUE); setwd(out_dir)
    tryCatch(utils::zip(zip_path,list.files(out_dir,full.names=FALSE),flags="-q",zip=Sys.which("zip")),error=function(e)NULL)
    setwd(oldwd); zip_ok <- file.exists(zip_path)
  }
  cat("\nOutput directory: ",out_dir,"\n",sep="")
  if(zip_ok) cat("Send this ZIP: ",zip_path,"\n",sep="") else message("ZIP unavailable; compress the generated output directory manually.")
  invisible(list(status=status,output_directory=out_dir,zip=if(zip_ok)zip_path else NULL))
}
main_alignment()
