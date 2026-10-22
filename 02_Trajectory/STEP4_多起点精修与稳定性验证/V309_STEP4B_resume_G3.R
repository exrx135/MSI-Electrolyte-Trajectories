# ============================================================
# V309 STEP 4B: resume V305 G3 resample refits after signature audit
# ============================================================
# Run on Windows desktop with Rscript V309_STEP4B_resume_G3.R.
# Requirements: data.table and gbmt 0.1.4 (same as V295).
# This standalone file contains the exact V295 data preparation and EM class
# kernels below, followed by a resampling controller. No outcome is read.
# 30 independent 80% subsamples + 10 ordinary bootstrap resamples.
# Up to four independent fits run concurrently on a 9800X3D (8C/16T),
# subject to free-memory headroom; each fit has three class-specific workers.
# Run again after interruption: the SAME V305 A/B checkpoints resume.
# The version of the per-draw results remains V305 by design. V309 permits
# exactly two text renderings of the same 1e-4 checkpoint tolerance value:
# "1e-04" and "0.0001". Every other signature field and structural check
# remains exact. Stop the old controller before starting this file.
# Edit the two task counts below before the first run if required.
# The original V295 result directory is read only.

options(stringsAsFactors = FALSE, scipen = 999)
VERSION_TAG <- "V305"
MODEL_DEGREE <- 2L
REFINE_TOL <- 1e-4
NONBINDING_MAXIT <- .Machine$integer.max
G <- 3L
N_SUBSAMPLES <- 30L
N_BOOTSTRAPS <- 10L
SUBSAMPLE_FRACTION <- 0.80
MASTER_SEED <- 20260927L
INNER_WORKERS <- 3L
MAX_OUTER_FITS <- 4L   # 4 x 3 busy EM workers, plus 4 fit coordinators
ESTIMATED_PEAK_MB_PER_FIT <- 3200  # V295 measured ~2.5 GB; margin included
RAM_HEADROOM_MB <- 4096             # keep OS and other applications responsive
EXPECTED_PATIENTS <- 9804L
EXPECTED_LONG_ROWS <- 117648L
EXPECTED_INPUT_MD5 <- "0361d8b59184b7a6b0a5fa9c946088a0"
Z_MARKERS <- c("msi_z", "sodium_z", "potassium_z", "chloride_z", "bicarbonate_z")
REQUIRED_COLS <- c("unit_id", "time_model", Z_MARKERS)
SCALING_TABLE <- data.frame(
  variable = Z_MARKERS,
  indicator = c("MSI", "Na", "K", "Cl", "HCO3"),
  mean_native = c(1.229422, 138.931869, 4.221718, 105.125243, 22.905997),
  sd_native = c(0.3219380, 5.3785897, 0.6426219, 6.5980074, 4.7412984)
)

if (.Platform$OS.type != "windows") stop("This script runs on Windows desktop.")
for (pkg in c("gbmt", "data.table")) {
  if (!requireNamespace(pkg, quietly = TRUE)) stop("Required R package missing: ", pkg)
}
if (as.character(packageVersion("gbmt")) != "0.1.4") {
  stop("gbmt must be version 0.1.4; installed: ", as.character(packageVersion("gbmt")))
}
suppressPackageStartupMessages(library(data.table))
data.table::setDTthreads(1L)
Sys.setenv(OMP_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1", MKL_NUM_THREADS="1",
           BLAS_NUM_THREADS="1", VECLIB_MAXIMUM_THREADS="1", NUMEXPR_NUM_THREADS="1")
profile_root <- Sys.getenv("USERPROFILE", unset="")
if (!nzchar(profile_root) || !dir.exists(profile_root)) stop("USERPROFILE is unavailable.")
project_dir <- file.path(profile_root, "Desktop", "Research", "MSI", "提取轨迹文件")
input_rds <- file.path(project_dir, "GBMTM_stage2_output", "STEP2 SAVE",
                       "04_RUN_CONFIGURATION", "V286_STEP2_model_input_source_V251_G2G20.rds")
reference_dir <- file.path(project_dir, "GBMTM_stage4_output", "desktop")
output_dir <- file.path(project_dir, "GBMTM_stage4_output", "V305_STEP4B_G3")
reference_patient_path <- file.path(reference_dir, "V295_patient_class_G3_S01.csv")
reference_traj_path <- file.path(reference_dir, "V295_trajectory_G3_S01.csv")
reference_summary_path <- file.path(reference_dir, "V295_task_summary_G3_S01.csv")
for (p in c(input_rds, reference_patient_path, reference_traj_path, reference_summary_path)) {
  if (!file.exists(p)) stop("Required exact input is missing: ", p)
}
input_md5 <- unname(tools::md5sum(input_rds))
if (!identical(tolower(input_md5), EXPECTED_INPUT_MD5)) stop("Input RDS MD5 mismatch; stop.")
reference_md5 <- unname(tools::md5sum(reference_patient_path))
dir.create(output_dir, recursive=TRUE, showWarnings=FALSE)
if (!dir.exists(output_dir)) stop("Cannot create Step 4B output folder: ", output_dir)
timestamp_now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_path <- file.path(output_dir,
                      sprintf("V306_progress_pid%d_G3.log", Sys.getpid()))
log_line <- function(...) {
  line <- paste0(timestamp_now(), " | ", paste0(..., collapse=""))
  cat(line, "\n"); cat(line, "\n", file=log_path, append=TRUE)
}
safe_rds <- function(obj, path, compress=FALSE) {
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  saveRDS(obj, tmp, compress=compress)
  if (!file.exists(tmp) || is.na(file.info(tmp)$size) || file.info(tmp)$size < 1L) stop("RDS write failed: ",tmp)
  if (file.exists(path)) file.remove(path)
  if (!file.rename(tmp,path)) stop("RDS rename failed: ", path)
  invisible(path)
}
safe_csv <- function(obj,path) {
  tmp <- paste0(path,".tmp_",Sys.getpid())
  data.table::fwrite(obj,tmp)
  if (!file.exists(tmp) || is.na(file.info(tmp)$size) || file.info(tmp)$size < 1L) stop("CSV write failed: ",tmp)
  if (file.exists(path)) file.remove(path)
  if (!file.rename(tmp,path)) stop("CSV rename failed: ",path)
  invisible(path)
}
task_paths <- function(task_id) {
  f <- function(content,ext) file.path(output_dir,paste0("V305_",content,"_",task_id,"_G3.",ext))
  list(map=f("patient_draw","csv"),
       checkpoint_a=f("checkpoint_A","rds"),
       checkpoint_b=f("checkpoint_B","rds"),
       history=f("EM_history","csv"),
       fit=f("compact_fit","rds"),
       quality=f("class_quality","csv"),
       patient=f("patient_class","csv"),
       traj=f("trajectory","csv"),
       alignment=f("alignment","csv"),
       comparison=f("comparison","csv"),
       class_metrics=f("class_metrics","csv"),
       done=f("DONE","txt"),
       error=f("ERROR","txt"))
}
ref_patient <- data.table::fread(reference_patient_path, colClasses=list(character="unit_id"))
ref_traj <- data.table::fread(reference_traj_path)
ref_summary <- data.table::fread(reference_summary_path)
if (!all(c("unit_id","assigned_class", "p_class1","p_class2","p_class3") %in% names(ref_patient)) ||
    !all(c("class","variable","time_model","fitted_z") %in% names(ref_traj)) ||
    !all(c("status","requested_g","input_md5") %in% names(ref_summary))) stop("V295 reference column mismatch.")
if (nrow(ref_summary)!=1L || ref_summary$status[1L]!="CONVERGED_VALID_FULL_G" ||
    ref_summary$requested_g[1L]!=G ||
    tolower(ref_summary$input_md5[1L])!=EXPECTED_INPUT_MD5) stop("V295 G3 S01 reference validation failed.")
if (nrow(ref_patient)!=EXPECTED_PATIENTS || anyDuplicated(ref_patient$unit_id) ||
    !setequal(ref_patient$assigned_class,seq_len(G))) stop("V295 reference patient table invalid.")
if (!all(c("indicator","time_hour_center") %in% names(ref_traj)) ||
    anyDuplicated(ref_traj[, .(class,variable,time_model)]) ||
    !setequal(ref_traj$class,seq_len(G)) ||
    !setequal(ref_traj$variable,Z_MARKERS)) stop("V295 reference trajectory invalid.")
input_dt <- as.data.table(readRDS(input_rds))
if (!setequal(names(input_dt),REQUIRED_COLS) || length(names(input_dt))!=length(REQUIRED_COLS)) stop("Frozen input columns mismatch.")
input_dt[, unit_id := as.character(unit_id)]
input_dt[, time_model := as.numeric(time_model)]
for (v in Z_MARKERS) input_dt[, (v) := as.numeric(get(v))]
if (nrow(input_dt)!=EXPECTED_LONG_ROWS || uniqueN(input_dt$unit_id)!=EXPECTED_PATIENTS ||
    anyNA(input_dt$unit_id) || anyNA(input_dt$time_model) ||
    input_dt[, anyDuplicated(.SD), .SDcols=c("unit_id","time_model")] != 0L) stop("Frozen input counts or keys invalid.")
all_ids <- sort(unique(input_dt$unit_id))
if (!setequal(all_ids, ref_patient$unit_id)) stop("V295 patient IDs differ from frozen Stage 2 input.")
input_row_index <- split(seq_len(nrow(input_dt)), input_dt$unit_id)
if (any(lengths(input_row_index)!=12L)) stop("Expected exactly 12 longitudinal rows per patient.")
setkey(ref_patient,unit_id)

plan_path <- file.path(output_dir,"V305_task_plan_G3.csv")
make_task_order <- function(n_sub,n_boot) {
  out <- character(0); s <- 1L; b <- 1L
  while(s<=n_sub || b<=n_boot) {
    if(s<=n_sub) {
      upto <- min(n_sub,s+2L)
      out <- c(out,sprintf("SUB%03d",seq.int(s,upto)))
      s <- upto+1L
    }
    if(b<=n_boot) {
      out <- c(out,sprintf("BOOT%03d",b))
      b <- b+1L
    }
  }
  out
}
expected_tasks <- make_task_order(N_SUBSAMPLES,N_BOOTSTRAPS)
expected_designs <- ifelse(startsWith(expected_tasks,"SUB"),"subsample","bootstrap")
if (!file.exists(plan_path)) {
  task_plan <- data.table(task_id=expected_tasks,design=expected_designs,
    seed=MASTER_SEED+seq_len(N_SUBSAMPLES+N_BOOTSTRAPS)*104729L,
    input_md5=input_md5,reference_patient_md5=reference_md5)
  safe_csv(task_plan,plan_path)
}
task_plan <- data.table::fread(plan_path)
if (!identical(as.character(task_plan$task_id),expected_tasks) ||
    !identical(as.character(task_plan$design),expected_designs) ||
    any(as.character(task_plan$input_md5)!=input_md5) ||
    any(as.character(task_plan$reference_patient_md5)!=reference_md5)) {
  stop("Existing V305 task plan differs from current inputs/configuration. No result was overwritten.")
}

# BEGIN EXACT V295 DATA PREPARATION / CLASS-SPECIFIC EM KERNELS.
prepare_static_data <- function(path) {
  raw_data <- readRDS(path)
  missing_cols <- setdiff(REQUIRED_COLS, names(raw_data))
  if (length(missing_cols) > 0L) {
    stop(
      paste0(
        "Model input is missing required column(s): ",
        paste(missing_cols, collapse = ", ")
      )
    )
  }
  data <- as.data.frame(raw_data)[, REQUIRED_COLS, drop = FALSE]
  rm(raw_data)
  data$unit_id <- as.character(data$unit_id)
  data$time_model <- as.numeric(data$time_model)
  for (v in Z_MARKERS) data[[v]] <- as.numeric(data[[v]])

  data[, "unit_id"] <- factor(data[, "unit_id"])
  cou <- levels(data[, "unit_id"])

  for (i in seq_along(cou)) {
    idx <- which(data[, "unit_id"] == cou[i])
    idat <- data[idx, , drop = FALSE]
    if (sum(duplicated(idat[, "time_model"])) > 0L) {
      stop(paste0("Unit '", cou[i], "' has duplicated time points."))
    }
    data[idx, ] <- idat[order(idat[, "time_model"]), , drop = FALSE]
  }

  tnam <- levels(factor(data[, "time_model"]))
  tval <- as.numeric(data[, "time_model"] - min(data[, "time_model"]) + 1)
  nt <- length(tnam)

  Xmat <- matrix(nrow = nt, ncol = MODEL_DEGREE + 1L)
  for (degree in 0:MODEL_DEGREE) {
    Xmat[, degree + 1L] <- sort(unique(tval))^degree
  }
  rownames(Xmat) <- tnam
  colnames(Xmat) <- paste0("I(t^", 0:MODEL_DEGREE, ")")

  dataI <- data
  ismiss <- sum(is.na(data[, Z_MARKERS])) > 0L
  if (ismiss) {
    for (v in Z_MARKERS) {
      idx <- which(is.na(data[, v]))
      if (length(idx) > 0L) {
        dataI[idx, v] <- mean(data[, v], na.rm = TRUE)
      }
    }
  }

  unit_index <- as.integer(data[, "unit_id"])
  row_indices <- split(seq_len(nrow(data)), unit_index)

  list(
    data = data,
    dataI = dataI,
    cou = cou,
    n = length(cou),
    np = length(Z_MARKERS),
    nt = nt,
    tnam = tnam,
    tval = tval,
    Xmat = Xmat,
    ismiss = ismiss,
    unit_index = unit_index,
    row_indices = row_indices,
    x_names = Z_MARKERS,
    unit_name = "unit_id",
    time_name = "time_model",
    degree = MODEL_DEGREE
  )
}


inner_init_chunk <- function(class_ids, Z, p0) {
  wlsFit <- getFromNamespace("wlsFit", "gbmt")
  logdmN <- getFromNamespace("logdmN", "gbmt")
  pldBound <- getFromNamespace("pldBound", "gbmt")
  s <- V295_STATIC
  out <- vector("list", length(class_ids))

  for (pos in seq_along(class_ids)) {
    class_id <- as.integer(class_ids[pos])
    iwei <- as.numeric(Z[s$unit_index] == class_id)
    imod <- wlsFit(
      x.names = s$x_names,
      data = s$dataI,
      t = s$tval,
      weights = iwei,
      d = s$degree
    )

    postd_col <- numeric(s$n)
    for (w in seq_len(s$n)) {
      idx <- s$row_indices[[w]]
      wdat <- as.matrix(s$dataI[idx, s$x_names, drop = FALSE])
      if (s$np > 1L) {
        ijD <- sum(logdmN(wdat, s$Xmat %*% imod$beta, imod$S))
      } else {
        ijD <- sum(
          dnorm(
            s$dataI[idx, s$x_names],
            s$Xmat %*% imod$beta,
            sqrt(imod$S),
            log = TRUE
          )
        )
      }
      postd_col[w] <- pldBound(ijD + log(p0[class_id]))
    }

    out[[pos]] <- list(
      class = class_id,
      beta = imod$beta,
      Sigma = imod$S,
      postd = postd_col
    )
  }
  out
}

inner_mstep_chunk <- function(class_ids, postp, beta, Sigma) {
  wlsFit <- getFromNamespace("wlsFit", "gbmt")
  imputFun <- getFromNamespace("imputFun", "gbmt")
  logdmN <- getFromNamespace("logdmN", "gbmt")
  pldBound <- getFromNamespace("pldBound", "gbmt")
  s <- V295_STATIC

  if (!exists("V295_DATAI_CACHE", envir = .GlobalEnv, inherits = FALSE)) {
    assign("V295_DATAI_CACHE", list(), envir = .GlobalEnv)
  }
  cache <- get("V295_DATAI_CACHE", envir = .GlobalEnv, inherits = FALSE)
  out <- vector("list", length(class_ids))

  for (pos in seq_along(class_ids)) {
    class_id <- as.integer(class_ids[pos])
    iwei <- postp[, class_id][s$unit_index]

    if (s$ismiss) {
      idat <- imputFun(
        xdat = s$data[, s$x_names, drop = FALSE],
        tval = s$tval,
        beta = beta[[class_id]],
        S = Sigma[[class_id]],
        d = s$degree
      )
    } else {
      idat <- s$data[, s$x_names, drop = FALSE]
    }

    imod <- wlsFit(
      x.names = s$x_names,
      data = idat,
      t = s$tval,
      weights = iwei,
      d = s$degree
    )

    postd_col <- numeric(s$n)
    for (w in seq_len(s$n)) {
      idx <- s$row_indices[[w]]
      wdat <- as.matrix(idat[idx, s$x_names, drop = FALSE])
      if (s$np > 1L) {
        ijD <- sum(logdmN(wdat, s$Xmat %*% imod$beta, imod$S))
      } else {
        ijD <- sum(
          dnorm(
            wdat,
            s$Xmat %*% imod$beta,
            sqrt(imod$S),
            log = TRUE
          )
        )
      }
      postd_col[w] <- pldBound(ijD + log(mean(postp[, class_id])))
    }

    cache[[as.character(class_id)]] <- idat
    out[[pos]] <- list(
      class = class_id,
      beta = imod$beta,
      Sigma = imod$S,
      postd = postd_col
    )
  }

  assign("V295_DATAI_CACHE", cache, envir = .GlobalEnv)
  out
}

flatten_class_results <- function(result_chunks, g, n) {
  flat <- unlist(result_chunks, recursive = FALSE)
  if (length(flat) != g) {
    stop(paste0("Parallel class result count ", length(flat), " differs from G=", g, "."))
  }
  beta <- vector("list", g)
  Sigma <- vector("list", g)
  postd <- matrix(NA_real_, nrow = n, ncol = g)
  seen <- integer(0)
  for (item in flat) {
    k <- as.integer(item$class)
    if (k < 1L || k > g || k %in% seen) stop("Invalid/duplicate class result.")
    seen <- c(seen, k)
    beta[[k]] <- item$beta
    Sigma[[k]] <- item$Sigma
    postd[, k] <- as.numeric(item$postd)
  }
  if (!identical(sort(seen), seq_len(g))) stop("One or more class results are missing.")
  list(beta = beta, Sigma = Sigma, postd = postd)
}


# END EXACT V295 DATA PREPARATION / CLASS-SPECIFIC EM KERNELS.

hungarian_min_cost <- function(cost) {
  cost <- as.matrix(cost)
  if (
    nrow(cost) == 0L ||
    nrow(cost) != ncol(cost) ||
    any(!is.finite(cost))
  ) {
    stop("Hungarian alignment requires a finite, nonempty square cost matrix.")
  }

  n <- nrow(cost)
  m <- ncol(cost)
  u <- numeric(n + 1L)
  v <- numeric(m + 1L)
  p <- integer(m + 1L)
  way <- integer(m + 1L)

  for (i in seq_len(n)) {
    p[1L] <- i
    j0 <- 1L
    minv <- rep(Inf, m + 1L)
    used <- rep(FALSE, m + 1L)

    repeat {
      used[j0] <- TRUE
      i0 <- p[j0]
      delta <- Inf
      j1 <- NA_integer_

      for (j in 2:(m + 1L)) {
        if (!used[j]) {
          cur <- cost[i0, j - 1L] - u[i0 + 1L] - v[j]
          if (cur < minv[j]) {
            minv[j] <- cur
            way[j] <- j0
          }
          if (minv[j] < delta) {
            delta <- minv[j]
            j1 <- j
          }
        }
      }

      if (!is.finite(delta) || is.na(j1)) {
        stop("Hungarian alignment failed to find an augmenting path.")
      }

      for (j in seq_len(m + 1L)) {
        if (used[j]) {
          u[p[j] + 1L] <- u[p[j] + 1L] + delta
          v[j] <- v[j] - delta
        } else {
          minv[j] <- minv[j] - delta
        }
      }

      j0 <- j1
      if (p[j0] == 0L) break
    }

    repeat {
      j1 <- way[j0]
      p[j0] <- p[j1]
      j0 <- j1
      if (j0 == 1L) break
    }
  }

  assignment <- integer(n)
  for (j in 2:(m + 1L)) {
    if (p[j] > 0L) assignment[p[j]] <- j - 1L
  }
  if (!identical(sort(assignment), seq_len(n))) {
    stop("Hungarian alignment returned an invalid permutation.")
  }
  assignment
}

adjusted_rand_index <- function(labels_a, labels_b) {
  if (length(labels_a) != length(labels_b) || length(labels_a) == 0L) {
    return(NA_real_)
  }
  tab <- table(labels_a, labels_b)
  choose2 <- function(x) x * (x - 1) / 2
  total_pairs <- choose2(sum(tab))
  if (total_pairs <= 0) return(NA_real_)

  sum_cells <- sum(choose2(tab))
  sum_rows <- sum(choose2(rowSums(tab)))
  sum_cols <- sum(choose2(colSums(tab)))
  expected <- sum_rows * sum_cols / total_pairs
  upper <- 0.5 * (sum_rows + sum_cols)
  denominator <- upper - expected

  if (abs(denominator) <= .Machine$double.eps) {
    return(if (all(labels_a == labels_b)) 1 else 0)
  }
  (sum_cells - expected) / denominator
}

# A saved task consists of a deterministic draw map, a newly estimated G3,
# and comparisons of the patients present in that draw with V295 G3 S01.
make_draw <- function(task) {
  set.seed(as.integer(task$seed))
  if (task$design == "subsample") {
    ids <- sample(all_ids, size=floor(length(all_ids)*SUBSAMPLE_FRACTION), replace=FALSE)
  } else if (task$design == "bootstrap") {
    ids <- sample(all_ids, size=length(all_ids), replace=TRUE)
  } else stop("Unknown draw design.")
  data.table(unit_id=sprintf("R%05d",seq_along(ids)), source_unit_id=ids,
             draw_position=seq_along(ids))
}
make_signature <- function(task, draw_md5) {
  paste(VERSION_TAG, task$task_id,task$design,task$seed,
        input_md5,reference_md5,draw_md5,MODEL_DEGREE,REFINE_TOL,
        "scaling=0", "pruning=TRUE", "delete.empty=FALSE",
        "nstart=1", "engine=V295_class_EM",sep="|")
}
save_cp <- function(state,paths) {
  p <- if (state$iteration %% 2L == 0L) paths$checkpoint_a else paths$checkpoint_b
  safe_rds(state,p,compress=FALSE)
  check <- readRDS(p)
  if (!identical(check$signature,state$signature) ||
      !identical(check$iteration,state$iteration)) stop("Checkpoint verification failed.")
}
load_cp <- function(paths,signature,n) {
  slots <- list(paths$checkpoint_a,paths$checkpoint_b)
  present <- vapply(slots,file.exists,logical(1))
  if (!any(present)) return(NULL)
  loaded <- lapply(slots[present],function(p) tryCatch(readRDS(p),error=function(e) NULL))
  signature_matches <- function(saved) {
    if (!is.character(saved) || length(saved)!=1L || is.na(saved)) return(FALSE)
    if (identical(saved,signature)) return(TRUE)
    current_fields <- strsplit(signature,"|",fixed=TRUE)[[1L]]
    saved_fields <- strsplit(saved,"|",fixed=TRUE)[[1L]]
    length(current_fields)==14L && length(saved_fields)==14L &&
      identical(current_fields[9L],"0.0001") &&
      identical(saved_fields[9L],"1e-04") &&
      identical(current_fields[-9L],saved_fields[-9L])
  }
  good <- vapply(loaded,function(z) {
    is.list(z) && signature_matches(z$signature) &&
    length(z$Z)==n && nrow(as.matrix(z$postd))==n &&
    ncol(as.matrix(z$postd))==G && length(z$beta)==G &&
    length(z$Sigma)==G && is.finite(z$ll_old) && is.finite(z$iteration)
  },logical(1))
  if (!any(good)) stop("Checkpoints exist but do not match this task. No reset performed.")
  loaded <- loaded[good]
  selected <- loaded[[which.max(vapply(loaded,function(z) as.integer(z$iteration),integer(1)))]]
  if (!identical(selected$signature,signature))
    log_line("Resuming checkpoint with equivalent 1e-04 tolerance spelling; all other signature fields match exactly.")
  selected
}
run_refit <- function(static,paths,seed,signature) {
  # Three workers correspond to the three class-specific V295 EM chunks.
  cl <- parallel::makePSOCKcluster(INNER_WORKERS, outfile="")
  on.exit(parallel::stopCluster(cl),add=TRUE)
  V295_STATIC <- static
  parallel::clusterExport(cl,"V295_STATIC",envir=environment())
  parallel::clusterExport(cl,c("inner_init_chunk","inner_mstep_chunk",
                               "Z_MARKERS","MODEL_DEGREE"),envir=.GlobalEnv)
  parallel::clusterEvalQ(cl,{
    Sys.setenv(OMP_NUM_THREADS="1",OPENBLAS_NUM_THREADS="1",MKL_NUM_THREADS="1",
               BLAS_NUM_THREADS="1",VECLIB_MAXIMUM_THREADS="1")
    suppressPackageStartupMessages(library(gbmt));
    if (exists("V295_DATAI_CACHE",envir=.GlobalEnv,inherits=FALSE))
      rm("V295_DATAI_CACHE",envir=.GlobalEnv)
    NULL
  })
  chunks <- as.list(seq_len(G))
  state <- load_cp(paths,signature,static$n)
  if (is.null(state)) {
    set.seed(seed)
    Z <- sample(rep(seq_len(G),length=static$n))
    names(Z) <- static$cou
    p0 <- prop.table(table(factor(Z,levels=seq_len(G))))
    init <- flatten_class_results(parallel::clusterMap(
      cl,fun=inner_init_chunk,class_ids=chunks,
      MoreArgs=list(Z=Z,p0=p0),SIMPLIFY=FALSE,USE.NAMES=FALSE,RECYCLE=FALSE),G,static$n)
    beta <- init$beta; Sigma <- init$Sigma; postd <- init$postd
    rownames(postd) <- static$cou
    ll_old <- sum(log(apply(exp(postd),1L,sum)))
    count <- 0L; count2 <- 0L
    em_history <- data.table(iteration=0L,logLik=ll_old,delta_logLik=NA_real_,
                             seconds=NA_real_,stop_reason="INITIALIZATION")
    save_cp(list(signature=signature,iteration=count,count2=count2,
                 ll_old=ll_old,p0=p0,beta=beta,Sigma=Sigma,postd=postd,
                 Z=Z,em_history=em_history),paths)
    log_line("fresh EM initialization; logLik=",format(ll_old,digits=16))
  } else {
    count <- as.integer(state$iteration); count2 <- as.integer(state$count2)
    ll_old <- as.numeric(state$ll_old); p0 <- state$p0
    beta <- state$beta; Sigma <- state$Sigma; postd <- as.matrix(state$postd)
    Z <- state$Z; names(Z) <- static$cou; rownames(postd) <- static$cou
    em_history <- as.data.table(state$em_history)
    log_line("resumed EM iteration ",count,"; logLik=",format(ll_old,digits=16))
  }
  safe_csv(em_history,paths$history)
  if (file.exists(file.path(output_dir,"V305_PAUSE_G3.txt")))
    return(list(paused=TRUE,iteration=count))
  fine <- 0L; stop_reason <- ""
  p0_old <- p0; beta_old <- beta; Sigma_old <- Sigma
  while (fine == 0L) {
    start_time <- Sys.time(); count <- count + 1L
    postp <- exp(postd); postp <- postp / rowSums(postp)
    if (any(!is.finite(postp))) stop("Nonfinite posterior during EM.")
    Z_old <- Z
    Z <- apply(postp,1L,which.max)
    p0 <- colMeans(postp)
    result <- flatten_class_results(parallel::clusterMap(
      cl,fun=inner_mstep_chunk,class_ids=chunks,
      MoreArgs=list(postp=postp,beta=beta,Sigma=Sigma),
      SIMPLIFY=FALSE,USE.NAMES=FALSE,RECYCLE=FALSE),G,static$n)
    beta <- result$beta; Sigma <- result$Sigma; postd <- result$postd
    rownames(postd) <- static$cou
    ll_new <- sum(log(apply(exp(postd),1L,sum)))
    if (!is.finite(ll_new)) stop("Nonfinite log likelihood during EM.")
    reasons <- character(0)
    if (ll_new >= ll_old) {
      if (ll_new-ll_old < REFINE_TOL) {
        fine <- 1L; reasons <- c(reasons,"DELTA_LL_BELOW_TOL")
      }
    } else {
      p0 <- p0_old; beta <- beta_old; Sigma <- Sigma_old
      ll_new <- ll_old; count <- count-1L; fine <- 1L
      reasons <- c(reasons,"LL_DECREASE_ROLLBACK")
    }
    if (count >= NONBINDING_MAXIT) {
      fine <- 1L; reasons <- c(reasons,"MAXIT")
    }
    # Replicates the exact V295 stable-assignment rule (including sum check).
    if (sum(Z_old-Z)==0L) {
      count2 <- count2+1L
      if (ll_new-ll_old < REFINE_TOL*sqrt(count2)) {
        fine <- 1L; reasons <- c(reasons,"STABLE_ASSIGNMENT_RULE")
      }
    } else count2 <- 0L
    stop_iter <- if (length(reasons)) paste(unique(reasons),collapse="+") else ""
    em_history <- rbind(em_history,data.table(iteration=count,logLik=ll_new,
         delta_logLik=ll_new-ll_old,
         seconds=as.numeric(difftime(Sys.time(),start_time,units="secs")),
         stop_reason=stop_iter))
    safe_csv(em_history,paths$history)
    if (fine==0L) {
      ll_old <- ll_new; p0_old <- p0; beta_old <- beta; Sigma_old <- Sigma
      save_cp(list(signature=signature,iteration=count,count2=count2,
                   ll_old=ll_old,p0=p0,beta=beta,Sigma=Sigma,postd=postd,
                   Z=Z,em_history=em_history),paths)
      if (count %% 10L == 0L) log_line("accepted EM ",count," logLik=",format(ll_old,digits=16))
      if (file.exists(file.path(output_dir,"V305_PAUSE_G3.txt")))
        return(list(paused=TRUE,iteration=count))
    } else stop_reason <- stop_iter
  }
  postp <- exp(postd); postp <- postp/rowSums(postp)
  Z <- apply(postp,1L,which.max)
  names(Z) <- static$cou
  cache_chunks <- parallel::clusterEvalQ(cl,{
    if (exists("V295_DATAI_CACHE",envir=.GlobalEnv,inherits=FALSE))
      get("V295_DATAI_CACHE",envir=.GlobalEnv) else list()
  })
  dataI_list <- vector("list",G)
  for (chunk in cache_chunks) for (nm in names(chunk))
    dataI_list[[as.integer(nm)]] <- chunk[[nm]]
  if (static$ismiss && any(vapply(dataI_list,is.null,logical(1))))
    stop("V295 class-specific EM cache incomplete at pruning.")
  if (static$ismiss) {
    auxdat <- dataI_list[[1L]][,static$x_names,drop=FALSE]*p0[1L]
    for (k in 2:G)
      auxdat <- auxdat + dataI_list[[k]][,static$x_names,drop=FALSE]*p0[k]
    newdat <- data.frame(static$data[,c(static$unit_name,static$time_name)],auxdat)
  } else {
    newdat <- static$data[,c(static$unit_name,static$time_name,static$x_names)]
  }
  pruneFun <- getFromNamespace("pruneFun","gbmt")
  nparCalc <- getFromNamespace("nparCalc","gbmt")
  icCalc <- getFromNamespace("icCalc","gbmt")
  appaCalc <- getFromNamespace("appaCalc","gbmt")
  pruned <- pruneFun(beta=beta,posterior=postp,data=newdat)
  beta <- pruned$beta; Sigma <- pruned$S
  ll_final <- pruned$logLik
  npar <- nparCalc(beta)
  ic <- icCalc(ll_final,npar=npar,ss=nrow(newdat))
  fitted <- vector("list",G)
  for (k in seq_len(G)) {
    fitted[[k]] <- data.frame(static$Xmat %*% beta[[k]])
    rownames(fitted[[k]]) <- static$tnam
    colnames(fitted[[k]]) <- static$x_names
  }
  names(p0) <- colnames(postp) <- names(beta) <- names(Sigma) <- names(fitted) <- seq_len(G)
  rownames(postp) <- static$cou
  fit <- list(call=list(x.names=static$x_names,unit=static$unit_name,time=static$time_name,
             ng=G,d=MODEL_DEGREE,scaling=0,pruning=TRUE,nstart=1L,
             tol=REFINE_TOL,maxit=NONBINDING_MAXIT,delete.empty=FALSE),
             prior=p0,beta=beta,Sigma=Sigma,posterior=round(postp,4),
             Xmat=static$Xmat,fitted=fitted,assign=Z,logLik=ll_final,
             npar=npar,ic=ic,appa=appaCalc(postp),
             em=c(logLik=ll_final,n.iter=count,converged=ifelse(count>=NONBINDING_MAXIT,0,1)))
  class(fit) <- "gbmt"
  list(paused=FALSE,fit=fit,iteration=count,stop_reason=stop_reason,
       em_logLik_before_pruning=ll_new)
}

trajectory_table <- function(fit,task_id) {
  x <- vector("list",G*length(Z_MARKERS))
  j <- 0L
  for (k in seq_len(G)) for (v in Z_MARKERS) {
    j <- j+1L
    m <- as.matrix(fit$fitted[[k]])
    row <- SCALING_TABLE[SCALING_TABLE$variable==v,,drop=FALSE]
    z <- as.numeric(m[,v])
    x[[j]] <- data.table(task_id=task_id,class=k,variable=v,
                         indicator=as.character(row$indicator),
                         time_model=as.numeric(rownames(m)),
                         time_hour_center=2*seq_along(z)-1,
                         fitted_z=z,
                         fitted_native=as.numeric(row$mean_native+row$sd_native*z))
  }
  rbindlist(x)
}
patient_and_quality <- function(fit,task_id) {
  post <- as.matrix(fit$posterior)
  lab <- as.integer(fit$assign[rownames(post)])
  if (length(lab)!=nrow(post) || anyNA(lab) || any(!is.finite(post)))
    stop("Refit posterior or assignment failed validation.")
  maxp <- apply(post,1L,max)
  patient <- data.table(task_id=task_id,unit_id=rownames(post),
                        assigned_class=lab,max_posterior=maxp,
                        p_class1=post[,1L],p_class2=post[,2L],p_class3=post[,3L])
  quality <- rbindlist(lapply(seq_len(G),function(k) {
    idx <- which(lab==k)
    n <- length(idx)
    ap <- if(n) mean(post[idx,k]) else NA_real_
    prior <- as.numeric(fit$prior[k])
    occ <- if(is.finite(ap) && ap>0 && ap<1 && is.finite(prior) && prior>0 && prior<1)
             ap/(1-ap)/(prior/(1-prior)) else NA_real_
    data.table(task_id=task_id,class=k,n=n,pct=100*n/length(lab),
               prior=prior,APPA=ap,OCC=occ,
               mean_max_posterior=if(n) mean(maxp[idx]) else NA_real_,
               pct_max_posterior_ge_090=if(n) 100*mean(maxp[idx]>=0.9) else NA_real_)
  }))
  list(patient=patient,quality=quality)
}
join_trajectories <- function(x,y) {
  # Fwrite/fread may round a floating-point time_model by one last digit.
  # Match its ordinal time slot and verify that the underlying hours agree.
  x <- copy(x); y <- copy(y)
  x[, time_slot := match(time_model,sort(unique(time_model)))]
  y[, time_slot := match(time_model,sort(unique(time_model)))]
  xy <- merge(x,y,by=c("variable","time_slot"))
  if (nrow(xy) &&
      (any(abs(xy$time_model.x-xy$time_model.y)>1e-8) ||
       any(xy$time_hour_center.x != xy$time_hour_center.y)))
    stop("Refit time grid does not match frozen V295 reference.")
  xy
}
observed_anchor_rows <- function(tab) {
  tab[variable=="msi_z" | time_hour_center %in% c(3,9,15,21)]
}
balanced_rmse <- function(differences,variable) {
  sqrt(mean(vapply(split(differences^2,variable),mean,numeric(1))))
}
compare_one <- function(fit,map,ref_patient,ref_traj,task_id,paths) {
  refit_traj <- trajectory_table(fit,task_id)
  extra <- patient_and_quality(fit,task_id)
  safe_csv(refit_traj,paths$traj)
  safe_csv(extra$patient,paths$patient)
  safe_csv(extra$quality,paths$quality)
  if (nrow(extra$patient)!=nrow(map) ||
      !setequal(extra$patient$unit_id,map$unit_id)) stop("Refit/draw patient mismatch.")
  pair <- extra$patient[match(map$unit_id,extra$patient$unit_id)]
  pair[, source_unit_id := map$source_unit_id]
  # A bootstrap draw can contain the same original person several times.
  # Report classification metrics once per unique original patient.
  unique_pair <- pair[, .(p_class1=mean(p_class1),p_class2=mean(p_class2),
                           p_class3=mean(p_class3)),by=source_unit_id]
  unique_pair[, fit_class := max.col(as.matrix(.SD),ties.method="first"),
              .SDcols=c("p_class1","p_class2","p_class3")]
  ref_rows <- ref_patient[match(unique_pair$source_unit_id,ref_patient$unit_id)]
  if (anyNA(ref_rows$assigned_class)) stop("Missing original patient class on draw.")
  unique_pair[, ref_class := as.integer(ref_rows$assigned_class)]
  n_unique <- nrow(unique_pair)
  cost <- matrix(NA_real_,G,G)
  for (s in seq_len(G)) for (r in seq_len(G)) {
    a <- unique_pair$fit_class==s; b <- unique_pair$ref_class==r
    union <- sum(a|b)
    overlap_cost <- if(union>0L) 1-sum(a&b)/union else 1
    xs <- observed_anchor_rows(refit_traj[class==s,
            .(variable,time_model,time_hour_center,fit_z=fitted_z)])
    yr <- observed_anchor_rows(ref_traj[class==r,
            .(variable,time_model,time_hour_center,ref_z=fitted_z)])
    xy <- join_trajectories(xs,yr)
    if(nrow(xy)!=28L) stop("Observed anchor grid does not match V295 reference.")
    trajectory_tiebreak <- balanced_rmse(xy$fit_z-xy$ref_z,xy$variable)
    cost[s,r] <- overlap_cost+0.001*trajectory_tiebreak
  }
  matching <- hungarian_min_cost(cost) # sample class -> reference class
  align <- data.table(task_id=task_id,fit_class=seq_len(G),
                      reference_class=matching,
                      matching_cost=cost[cbind(seq_len(G),matching)])
  safe_csv(align,paths$alignment)
  unique_pair[, aligned_class := matching[fit_class]]
  ari <- adjusted_rand_index(unique_pair$ref_class,unique_pair$aligned_class)
  all_dist <- list(); class_rows <- list(); offset <- 0L
  for(s in seq_len(G)) {
    r <- matching[s]
    original <- unique_pair$ref_class==r
    fitted <- unique_pair$fit_class==s
    denom <- sum(original|fitted)
    jac <- if(denom>0L) sum(original&fitted)/denom else NA_real_
    rpct <- 100*mean(ref_patient$assigned_class==r)
    fpct <- 100*mean(as.integer(fit$assign)==s)
    appa <- extra$quality$APPA[extra$quality$class==s]
    for(v in c("ALL",Z_MARKERS)) {
      x <- observed_anchor_rows(refit_traj[class==s & (v=="ALL" | variable==v),
                .(variable,time_model,time_hour_center,z=fitted_z,native=fitted_native)])
      y <- observed_anchor_rows(ref_traj[class==r & (v=="ALL" | variable==v),
                .(variable,time_model,time_hour_center,z0=fitted_z,native0=fitted_native)])
      joined <- join_trajectories(x,y)
      expected <- if(v=="ALL") 28L else if(v=="msi_z") 12L else 4L
      if (nrow(joined)!=expected) stop("Trajectory matching has missing time points.")
      joined[, `:=`(cz=z-z[1L], cz0=z0-z0[1L]),by=variable]
      # RMSE in native units is meaningful for one marker only.
      native_rmse <- if(v=="ALL") NA_real_ else sqrt(mean((joined$native-joined$native0)^2))
      offset <- offset+1L
      class_rows[[offset]] <- data.table(task_id=task_id,fit_class=s,reference_class=r,
        variable=v,jaccard=jac,
        n_reference_shared=sum(original),n_fit_shared=sum(fitted),
        n_overlap=sum(original&fitted),
        reference_pct_full=rpct,fit_pct_draw=fpct,proportion_drift_pp=fpct-rpct,
        APPA=appa,rmse_z=balanced_rmse(joined$z-joined$z0,joined$variable),
        rmse_centered_z=balanced_rmse(joined$cz-joined$cz0,joined$variable),
        rmse_native=native_rmse,
        correlation_native=if(sd(joined$z)>0 && sd(joined$z0)>0)
          cor(joined$z,joined$z0) else NA_real_,
        correlation_centered=if(sd(joined$cz)>0 && sd(joined$cz0)>0)
          cor(joined$cz,joined$cz0) else NA_real_)
    }
  }
  classes <- rbindlist(class_rows)
  safe_csv(classes,paths$class_metrics)
  maxp <- extra$patient$max_posterior
  post <- as.matrix(fit$posterior)
  positive <- post[post>0]
  relative_entropy <- 1+sum(positive*log(positive))/(nrow(post)*log(G))
  stat <- data.table(task_id=task_id,status=if(all(extra$quality$n>0L))
                       "CONVERGED_VALID_FULL_G" else "CONVERGED_EMPTY_CLASS",
    requested_g=G,n_drawn=nrow(map),n_unique=n_unique,
    final_logLik=as.numeric(fit$logLik),em_iterations=as.integer(fit$em["n.iter"]),
    ARI=ari,mean_class_Jaccard=mean(classes[variable=="ALL",jaccard]),
    min_class_Jaccard=min(classes[variable=="ALL",jaccard]),
    mean_class_rmse_z=mean(classes[variable=="ALL",rmse_z]),
    mean_class_rmse_centered_z=mean(classes[variable=="ALL",rmse_centered_z]),
    max_abs_proportion_drift_pp=max(abs(classes[variable=="ALL",proportion_drift_pp])),
    min_class_pct=min(extra$quality$pct),min_APPA=min(extra$quality$APPA),
    min_OCC=if(all(is.na(extra$quality$OCC))) NA_real_ else min(extra$quality$OCC,na.rm=TRUE),
    entropy=relative_entropy,mean_max_posterior=mean(maxp),
    median_max_posterior=median(maxp),
    pct_max_posterior_ge_070=100*mean(maxp>=0.7),
    pct_max_posterior_ge_090=100*mean(maxp>=0.9))
  list(summary=stat,classes=classes)
}
refresh_aggregate <- function() {
  summary_files <- Sys.glob(file.path(output_dir,"V305_comparison_*_G3.csv"))
  if (length(summary_files)) {
    report <- rbindlist(lapply(summary_files,data.table::fread),fill=TRUE)
    setorder(report,task_id)
    safe_csv(report,file.path(output_dir,"V305_summary_G3.csv"))
  }
  class_files <- Sys.glob(file.path(output_dir,"V305_class_metrics_*_G3.csv"))
  if (length(class_files)) {
    all_classes <- rbindlist(lapply(class_files,data.table::fread),fill=TRUE)
    setorder(all_classes,task_id,reference_class,variable)
    safe_csv(all_classes,file.path(output_dir,"V305_class_summary_G3.csv"))
  }
}
run_one_draw <- function(task) {
  task_id <- as.character(task$task_id)
  paths <- task_paths(task_id)
  required_outputs <- c(paths$done,paths$fit,paths$comparison,paths$class_metrics,
                        paths$alignment,paths$quality)
  if (file.exists(paths$done)) {
    if (!all(file.exists(required_outputs))) stop("DONE marker has incomplete outputs: ",task_id)
    log_line(task_id," already complete: skip")
    return(invisible(NULL))
  }
  draw <- make_draw(task)
  if (file.exists(paths$map)) {
    old <- data.table::fread(paths$map,colClasses=list(character=c("unit_id","source_unit_id")))
    if (!identical(old$unit_id,draw$unit_id) ||
        !identical(old$source_unit_id,draw$source_unit_id))
      stop("Saved draw map does not match deterministic seed: ",task_id)
  } else safe_csv(draw,paths$map)
  draw_md5 <- unname(tools::md5sum(paths$map))
  signature <- make_signature(task,draw_md5)
  sample_path <- file.path(output_dir,paste0("V305_model_input_",task_id,"_G3.rds"))
  if (file.exists(paths$fit)) {
    fit <- readRDS(paths$fit)
    info <- fit$V305_metadata
    if (is.null(info) || !identical(info$task_id,task_id) ||
        !identical(info$draw_md5,draw_md5) ||
        !identical(info$source_input_md5,input_md5) ||
        !identical(info$reference_patient_md5,reference_md5))
      stop("Saved fit has a mismatched task identity: ",task_id)
    result <- list(stop_reason=info$stop_reason,
                   em_logLik_before_pruning=info$em_logLik_before_pruning)
    log_line("reusing already saved refit for downstream results: ",task_id)
  } else {
    row_lists <- input_row_index[draw$source_unit_id]
    if (any(lengths(row_lists)!=12L)) stop("Malformed draw rows: ",task_id)
    row_index <- unlist(row_lists,use.names=FALSE)
    sample_data <- as.data.frame(input_dt[row_index,..REQUIRED_COLS])
    sample_data$unit_id <- rep(draw$unit_id,times=lengths(row_lists))
    # The input can be reconstructed from the draw map; it is removed after DONE.
    safe_rds(sample_data,sample_path,compress=FALSE)
    rm(sample_data,row_lists,row_index)
    static <- prepare_static_data(sample_path)
    if (static$n!=nrow(draw) || static$nt!=12L) stop("Refit design matrix invalid: ",task_id)
    log_line("starting ",task_id," design=",task$design," draw=",nrow(draw),
             " unique=",uniqueN(draw$source_unit_id))
    result <- run_refit(static,paths,as.integer(task$seed),signature)
    rm(static); gc(verbose=FALSE)
    if (isTRUE(result$paused)) {
      log_line("PAUSED after checkpoint at ",task_id," iteration=",result$iteration)
      return(invisible(FALSE))
    }
    fit <- result$fit
    if (fit$call$ng!=G || ncol(fit$posterior)!=G || length(fit$beta)!=G)
      stop("Refit lost its requested full-G parameter structure: ",task_id)
    fit$V305_metadata <- list(version=VERSION_TAG,task_id=task_id,
            design=as.character(task$design),draw_md5=draw_md5,
            source_input_md5=input_md5,reference_patient_md5=reference_md5,
            stop_reason=result$stop_reason,
            em_logLik_before_pruning=result$em_logLik_before_pruning)
    safe_rds(fit,paths$fit,compress=TRUE)
  }
  if (fit$call$ng!=G || ncol(fit$posterior)!=G || length(fit$beta)!=G)
    stop("Saved fit has invalid full-G structure: ",task_id)
  result_metrics <- compare_one(fit,draw,ref_patient,ref_traj,task_id,paths)
  result_metrics$summary[, `:=`(design=as.character(task$design),seed=as.integer(task$seed),
         stop_reason=result$stop_reason,
         em_logLik_before_pruning=result$em_logLik_before_pruning,
         input_md5=input_md5,reference_patient_md5=reference_md5,draw_md5=draw_md5,
         completed_at=timestamp_now())]
  safe_csv(result_metrics$summary,paths$comparison)
  writeLines(c(paste0("Task: ",task_id),
               paste0("Status: ",result_metrics$summary$status),
               paste0("Completed: ",timestamp_now()),
               paste0("ARI: ",result_metrics$summary$ARI)),paths$done)
  if (file.exists(paths$error)) file.remove(paths$error)
  if (file.exists(paths$checkpoint_a)) file.remove(paths$checkpoint_a)
  if (file.exists(paths$checkpoint_b)) file.remove(paths$checkpoint_b)
  if (file.exists(sample_path)) file.remove(sample_path)
  log_line("completed ",task_id,"; ARI=",round(result_metrics$summary$ARI,4),
           " minJ=",round(result_metrics$summary$min_class_Jaccard,4),
           " minAPPA=",round(result_metrics$summary$min_APPA,4))
  # Only the parent V306 controller writes the shared aggregate CSVs.
  invisible(TRUE)
}

v306_script_path <- function() {
  args <- commandArgs(trailingOnly=FALSE)
  ff <- grep("^--file=", args, value=TRUE)
  if (length(ff)) {
    f <- sub("^--file=", "", ff[1L])
    if (file.exists(f)) return(normalizePath(f, winslash="/", mustWork=TRUE))
  }
  frames <- sys.frames()
  for (j in rev(seq_along(frames))) {
    f <- frames[[j]]$ofile
    if (!is.null(f) && length(f)==1L && file.exists(f) &&
        identical(basename(f), "V309_STEP4B_resume_G3.R"))
      return(normalizePath(f, winslash="/", mustWork=TRUE))
  }
  stop("Run with Rscript <full_path_to_V309_STEP4B_resume_G3.R> ",
       "or source(file.choose()) after selecting the complete V309 file.")
}

v306_worker_task <- function(task) {
  task_id <- as.character(task$task_id)
  if (file.exists(file.path(output_dir,"V305_PAUSE_G3.txt")))
    return(list(task_id=task_id, status="PAUSED_BEFORE_START",
                message="Remove V305_PAUSE_G3.txt before restarting."))
  tryCatch({
    result <- run_one_draw(task)
    list(task_id=task_id,
         status=if (identical(result,FALSE)) "PAUSED_AT_CHECKPOINT" else "OK",
         message="")
  }, error=function(e) {
    msg <- paste0(timestamp_now(), " | ",task_id," | ",conditionMessage(e))
    writeLines(msg, task_paths(task_id)$error)
    log_line("ERROR: ",msg)
    list(task_id=task_id, status="ERROR", message=msg)
  })
}

v306_available_ram_mb <- function() {
  shell <- Sys.which("powershell")
  if (!nzchar(shell)) shell <- Sys.which("pwsh")
  if (!nzchar(shell)) return(NA_real_)
  out <- tryCatch(suppressWarnings(system2(shell,
    args=c("-NoProfile","-Command",
           "\"(Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory\""),
    stdout=TRUE,stderr=FALSE)),error=function(e) character(0))
  values <- suppressWarnings(as.numeric(trimws(out)))
  values <- values[is.finite(values)&values>0]
  if (!length(values)) NA_real_ else tail(values,1L)/1024
}

v306_main <- function() {
  log_line("V306 controller started; existing V305 fit/checkpoint format preserved; input MD5=",
           input_md5,"; G=",G,"; subsamples=",N_SUBSAMPLES,
           "; bootstrap=",N_BOOTSTRAPS)
  if (file.exists(file.path(output_dir,"V305_PAUSE_G3.txt")))
    stop("V305_PAUSE_G3.txt is present. Remove it only after the old V305 process exits.")
  refresh_aggregate()
  completed <- vapply(task_plan$task_id,function(id) {
    p <- task_paths(as.character(id))
    all(file.exists(c(p$done,p$fit,p$comparison,p$class_metrics,
                      p$alignment,p$quality)))
  },logical(1))
  remaining <- task_plan[!completed]
  if (!nrow(remaining)) {
    log_line("All V305 resample fits already complete; nothing to rerun.")
    return(invisible(NULL))
  }
  logical_cpus <- suppressWarnings(as.integer(parallel::detectCores(logical=TRUE)))
  if (length(logical_cpus)!=1L || is.na(logical_cpus) || logical_cpus<1L)
    logical_cpus <- 16L  # expected desktop; edit MAX_OUTER_FITS for another CPU
  free_ram_mb <- v306_available_ram_mb()
  ram_fit_limit <- if(is.finite(free_ram_mb))
    max(1L,as.integer(floor(max(0,free_ram_mb-RAM_HEADROOM_MB)/
                             ESTIMATED_PEAK_MB_PER_FIT))) else 3L
  outer_n <- min(MAX_OUTER_FITS,ram_fit_limit,
                 max(1L,as.integer(floor(logical_cpus/(INNER_WORKERS+1L)))),
                 nrow(remaining))
  log_line("Logical CPUs=",logical_cpus,"; outer simultaneous fits=",outer_n,
           "; class workers per fit=",INNER_WORKERS,
           "; currently available RAM MB=",round(free_ram_mb),
           "; RAM-based fit limit=",ram_fit_limit,
           "; fits remaining=",nrow(remaining),
           ". One coordinator per fit is mostly idle while EM classes run.")
  script <- v306_script_path()
  outer <- parallel::makePSOCKcluster(outer_n, outfile="")
  on.exit(parallel::stopCluster(outer),add=TRUE)
  worker_ready <- parallel::clusterCall(outer,function(path) {
    options(V306_IS_WORKER=TRUE)
    sys.source(path,envir=.GlobalEnv)
    list(pid=Sys.getpid(),gbmt=as.character(utils::packageVersion("gbmt")))
  }, script)
  if (length(worker_ready)!=outer_n ||
      any(vapply(worker_ready,function(x) x$gbmt!="0.1.4",logical(1))))
    stop("Outer workers did not load the correct gbmt engine.")
  tasks <- lapply(seq_len(nrow(remaining)),function(i) remaining[i])
  # Dynamic scheduling keeps healthy workers busy when tasks have very
  # different run times. Every task has its own fixed seed and file identity.
  result <- parallel::parLapplyLB(outer,tasks,v306_worker_task,chunk.size=1L)
  refresh_aggregate()
  statuses <- vapply(result,function(x) x$status,character(1))
  log_line("V306 controller ended. Fits OK=",sum(statuses=="OK"),
           "; paused=",sum(grepl("PAUSED",statuses)),
           "; errors=",sum(statuses=="ERROR"),
           ". Review V305_summary_G3.csv and per-task ERROR files.")
  if (any(statuses=="ERROR"))
    warning("Some draws failed. Their checkpoints were retained; fix the cause and rerun V306.")
  invisible(result)
}

if (!isTRUE(getOption("V306_IS_WORKER", FALSE))) v306_main()
