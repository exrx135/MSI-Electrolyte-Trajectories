# V316: install the supplied eICU conversion, then replay frozen MIMIC G3 S01.
# Extract V316_STEP5_update_G3.zip directly under GBMTM_stage5_output.
# In RStudio run source(file.choose(), encoding="UTF-8") and select this file.
# Existing V295 models, MIMIC/SICdb data and finished results are untouched.

main_update <- function() {
  old_options <- options(stringsAsFactors=FALSE, error=NULL)
  on.exit(options(old_options), add=TRUE)
  profile <- Sys.getenv("USERPROFILE", unset="")
  if(!nzchar(profile) || !dir.exists(profile)) {
    profile <- path.expand("~")
    if(identical(basename(profile),"Documents"))
      stop("USERPROFILE is unavailable and ~ points to Documents.")
  }
  output_root <- file.path(profile,"Desktop","Research","MSI",
                           "提取轨迹文件","GBMTM_stage5_output")
  kit_dir <- file.path(output_root,"V307_STEP5_laptop_kit_G3")
  update_dir <- file.path(output_root,"V316_STEP5_update_G3")
  if(!dir.exists(kit_dir)) stop("Exact input kit does not exist: ",kit_dir)
  if(!dir.exists(update_dir))
    stop("Extract V316_STEP5_update_G3.zip under this directory first: ",output_root)

  frozen <- c(
    "MIMIC_FINAL_MAKE30_ELECTROLYTEQC.csv"="f1eb2d7fac1303a682b53c679dde9abd",
    "SICdb_FINAL_MAKE30_ELECTROLYTEQC.csv"="d52a3ecdbb636bba2859eb9f102fdf45",
    "V295_compact_fit_G3_S01.rds"="cf5881479ee206195aa80d59fa3bd17f",
    "V295_patient_class_G3_S01.csv"="dd12a3730af44c1035c79f4558c14db9",
    "V295_trajectory_G3_S01.csv"="3d47aa52a546266e937f50a26de2676a",
    "V295_task_summary_G3_S01.csv"="6b8eb067bb4013a54f6cb4dafe6d5f83"
  )
  fp <- file.path(kit_dir,names(frozen))
  if(any(!file.exists(fp))) stop("An exact frozen model/input file is missing.")
  if(!identical(unname(tools::md5sum(fp)),unname(frozen)))
    stop("Frozen MIMIC/SICdb/model identity mismatch; no files have been replaced.")
  old_hash <- "66650775841cfd26c78a368adcd9c805"
  new_hash <- "a06fddbcbc3b93f2c83227201c1138cd"
  payload <- file.path(update_dir,"V316_eICU_validation_G3.csv")
  target <- file.path(kit_dir,"eICU_FINAL_RAW_ELECTROLYTEQC.csv")
  if(!file.exists(payload) || !identical(unname(tools::md5sum(payload)),new_hash))
    stop("Supplied V316 eICU conversion file is missing or has changed.")
  if(!file.exists(target)) stop("The exact eICU input file does not exist: ",target)
  current_hash <- unname(tools::md5sum(target))
  if(!current_hash %in% c(old_hash,new_hash))
    stop("Existing eICU input is neither the original V307 file nor this V316 conversion.")

  audit_dir <- file.path(output_root,"V316_eICU_audit_G3")
  dir.create(audit_dir,recursive=TRUE,showWarnings=FALSE)
  if(!dir.exists(audit_dir)) stop("Cannot create eICU audit directory.")
  backup <- file.path(audit_dir,"V316_eICU_original_G3.csv")
  if(identical(current_hash,old_hash)) {
    if(file.exists(backup)) {
      if(!identical(unname(tools::md5sum(backup)),old_hash))
        stop("Existing original-data backup has a different hash.")
    } else {
      if(!file.copy(target,backup,overwrite=FALSE) ||
         !identical(unname(tools::md5sum(backup)),old_hash))
        stop("Could not preserve the original eICU input before replacement.")
    }
    stage <- tempfile("V316_eICU_",tmpdir=kit_dir,fileext=".csv")
    on.exit(unlink(stage),add=TRUE)
    if(!file.copy(payload,stage,overwrite=FALSE) ||
       !identical(unname(tools::md5sum(stage)),new_hash))
      stop("Could not stage and validate the replacement eICU file.")
    if(!file.copy(stage,target,overwrite=TRUE) ||
       !identical(unname(tools::md5sum(target)),new_hash)) {
      restored <- file.copy(backup,target,overwrite=TRUE)
      stop("eICU replacement failed. Original backup: ",backup,
           "; original restored: ",restored)
    }
    message("eICU replaced: 5867 original -> 5839 retained; 28 time-rule exclusions.")
  } else {
    message("V316 eICU data already installed; replacement skipped.")
  }

  for(filename in c("V316_eICU_excluded_G3.csv","V316_eICU_flow_G3.csv")) {
    src <- file.path(update_dir,filename)
    dst <- file.path(audit_dir,filename)
    if(!file.exists(src)) stop("Missing exact exclusion audit file: ",filename)
    if(!file.copy(src,dst,overwrite=TRUE)) stop("Cannot save exclusion audit: ",dst)
  }
  for(filename in c("V316_STEP5_replay_G3.R","V316_STEP5_transport_make30_G3.R")) {
    src <- file.path(update_dir,filename)
    dst <- file.path(kit_dir,filename)
    if(!file.exists(src)) stop("Missing exact complete R script: ",filename)
    if(!file.copy(src,dst,overwrite=TRUE) ||
       !identical(unname(tools::md5sum(src)),unname(tools::md5sum(dst))))
      stop("Cannot install complete R script: ",filename)
  }
  if(!identical(unname(tools::md5sum(fp)),unname(frozen)))
    stop("Frozen model/input identity changed during installation.")
  message("G3 model and MIMIC/SICdb hashes verified unchanged.")
  message("Starting corrected replay. Outputs: ",output_root)
  result <- source(file.path(kit_dir,"V316_STEP5_replay_G3.R"),encoding="UTF-8")
  invisible(result$value)
}

tryCatch(main_update(),error=function(e) {
  message("V316 update stopped: ",conditionMessage(e))
  invisible(list(ok=FALSE,error=conditionMessage(e)))
})
