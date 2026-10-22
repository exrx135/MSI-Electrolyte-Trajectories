# ============================================================
# V295 - STEP 4A DESKTOP TWO-LANE G3-G5 MULTIPLE-START REFINEMENT
# ============================================================
#
# This is the complete desktop continuation copy for low-G STEP 4A refinement.
# It runs G3, G4, and G5 only, with five independent starts per G, and stops
# unless Windows reports the frozen Ryzen 7 9800X3D desktop specification.
# Existing V295 G3-G5 DONE files are skipped and valid V295 A/B checkpoints are
# resumed without changing their signatures, seeds, filenames, or file types.
#
# A task means one complete fit for one G x one independent start_id.
# Each fit uses nstart=1 internally. Multiple independent fits are compared
# after completion.
#
# Scientific settings:
#   - five indicators: MSI + Na + K + Cl + HCO3
#   - d = 2
#   - scaling = 0
#   - pruning = TRUE
#   - delete.empty = FALSE
#   - tol = 1e-4
#   - maxit = .Machine$integer.max (non-binding technical sentinel)
#   - input NA values are preserved; gbmt 0.1.4 internal missing-value
#     handling is retained exactly as in the V258 checkpoint engine
#   - no outcome file is read or used
#
# Execution settings:
#   - Windows PSOCK multiprocessing
#   - the already completed V295 G6-S01 is read only as the measured memory
#     reference; this continuation script never schedules another G6 fit
#   - pending G3-G5 starts use exactly two simultaneous outer lanes
#   - each fit parallelizes the class-specific part of every EM iteration
#   - desktop CPU order: one thread on every physical core before SMT
#   - BLAS/OpenMP nesting is disabled inside every worker
#   - accepted EM state is saved after every accepted non-final iteration
#   - alternating A/B checkpoints permit iteration-level restart
#   - a pause request is honored only after a safe checkpoint boundary
#
# GPU note:
#   gbmt 0.1.4 has no CUDA/GPU implementation. The RTX 5070 Ti is therefore
#   not used by this statistically identical CPU refinement engine.
#
# Required archived local file:
#   %USERPROFILE%/Desktop/Research/MSI/提取轨迹文件/
#     GBMTM_stage2_output/STEP2 SAVE/04_RUN_CONFIGURATION/
#       V286_STEP2_model_input_source_V251_G2G20.rds
#
# Output:
#   %USERPROFILE%/Desktop/Research/MSI/提取轨迹文件/
#     GBMTM_stage4_output/desktop/
#
# ============================================================

options(stringsAsFactors = FALSE, scipen = 999)


# ============================================================
# 0. Frozen configuration
# ============================================================

VERSION_TAG <- "V295"
ENGINE_TAG <- "V295_parallel_EM_checkpoint"
EXPECTED_MACHINE_ROLE <- "desktop"

MODEL_DEGREE <- 2L
REFINE_TOL <- 1e-4
NONBINDING_MAXIT <- .Machine$integer.max

Z_MARKERS <- c(
  "msi_z",
  "sodium_z",
  "potassium_z",
  "chloride_z",
  "bicarbonate_z"
)

REQUIRED_COLS <- c(
  "unit_id",
  "time_model",
  Z_MARKERS
)

EXPECTED_PATIENTS <- 9804L
EXPECTED_LONG_ROWS <- 117648L
EXPECTED_TIME_POINTS <- 12L
EXPECTED_INPUT_MD5 <- "0361d8b59184b7a6b0a5fa9c946088a0"

STAGE2_MASTER_SEED <- 20260904L
START_SEED_STEP <- 1009L

CHECKPOINT_STATE_VERSION <- 2L
CHECKPOINT_EVERY_ACCEPTED_ITER <- 1L

MAX_SUPPORTED_SIMULTANEOUS_FITS <- 2L
CALIBRATION_G <- 6L
CALIBRATION_START_ID <- 1L
START_STAGGER_SECONDS <- 20
FORCE_TWO_LANES <- TRUE

# New fits are held back when system memory is already constrained. Once at
# least one completed task has supplied an observed process-tree peak, that
# observed peak is used by the gate. No guessed per-fit GB value is used.
MEMORY_MIN_AVAILABLE_GB <- 6
MEMORY_MIN_AVAILABLE_FRACTION <- 0.20
MEMORY_MAX_COMMIT_FRACTION <- 0.86
MEMORY_PEAK_SAFETY_FACTOR <- 1.25
MEMORY_POLL_SECONDS <- 15
MEMORY_POST_CALIBRATION_SETTLE_SECONDS <- 5

PROCESS_PRIORITY_CLASS <- "AboveNormal"

PREFERRED_MIN_CLASS_PCT <- 5
PREFERRED_MIN_APPA <- 0.70

SCALING_TABLE <- data.frame(
  variable = Z_MARKERS,
  indicator = c("MSI", "Na", "K", "Cl", "HCO3"),
  mean_native = c(
    1.229422,
    138.931869,
    4.221718,
    105.125243,
    22.905997
  ),
  sd_native = c(
    0.3219380,
    5.3785897,
    0.6426219,
    6.5980074,
    4.7412984
  ),
  stringsAsFactors = FALSE
)

# Frozen Stage-2 nstart=1 reference values extracted by V270 from the
# completed V251 fits supplied for this project. They are audit benchmarks
# only and never initialize or influence a V295 fit.
STAGE2_BENCHMARK <- data.frame(
  requested_g = c(3L, 4L, 5L, 6L),
  stage2_tol = rep(NA_real_, 4L),
  stage2_em_iterations = c(59L, 65L, 103L, 139L),
  stage2_logLik = c(
    -324711.964417534,
    -314159.304266617,
    -308197.159933641,
    -300819.052704939
  ),
  stage2_BIC = c(
    650474.719550568,
    629731.338272962,
    618145.637726447,
    603751.362293270
  ),
  stringsAsFactors = FALSE
)


# ============================================================
# 1. Packages and single-thread nesting guard
# ============================================================

required_packages <- c("data.table", "gbmt")

missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]

if (length(missing_packages) > 0L) {
  stop(
    paste0(
      "Missing required R package(s): ",
      paste(missing_packages, collapse = ", "),
      "\nInstall the exact project environment before running V295."
    )
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(gbmt)
  library(parallel)
})

if (!identical(as.character(packageVersion("gbmt")), "0.1.4")) {
  stop(
    paste0(
      "V295 is source-matched to gbmt 0.1.4.\n",
      "Installed gbmt version: ",
      as.character(packageVersion("gbmt")),
      "\nStop here rather than silently changing the fitting engine."
    )
  )
}

required_gbmt_internals <- c(
  "wlsFit",
  "imputFun",
  "logdmN",
  "pldBound",
  "pruneFun",
  "nparCalc",
  "icCalc",
  "appaCalc"
)

missing_gbmt_internals <- required_gbmt_internals[
  !vapply(
    required_gbmt_internals,
    function(nm) {
      !inherits(try(getFromNamespace(nm, "gbmt"), silent = TRUE), "try-error")
    },
    logical(1)
  )
]

if (length(missing_gbmt_internals) > 0L) {
  stop(
    paste0(
      "Required gbmt 0.1.4 internal function(s) unavailable: ",
      paste(missing_gbmt_internals, collapse = ", ")
    )
  )
}

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

data.table::setDTthreads(1L)

if (.Platform$OS.type != "windows") {
  stop(
    paste0(
      "V295 is written for the frozen Windows desktop computer.\n",
      "Detected OS type: ",
      .Platform$OS.type
    )
  )
}


# ============================================================
# 2. Exact relative project paths
# ============================================================

profile_root <- Sys.getenv("USERPROFILE", unset = "")

if (!nzchar(profile_root) || !dir.exists(profile_root)) {
  stop(
    paste0(
      "USERPROFILE did not resolve to an existing Windows profile directory:\n",
      normalizePath(profile_root, winslash = "/", mustWork = FALSE)
    )
  )
}

project_dir <- file.path(
  profile_root,
  "Desktop",
  "Research",
  "MSI",
  "提取轨迹文件"
)

stage2_dir <- file.path(
  project_dir,
  "GBMTM_stage2_output",
  "STEP2 SAVE"
)

input_rds <- file.path(
  stage2_dir,
  "04_RUN_CONFIGURATION",
  "V286_STEP2_model_input_source_V251_G2G20.rds"
)

output_root <- file.path(project_dir, "GBMTM_stage4_output")

if (!dir.exists(project_dir)) {
  stop(
    paste0(
      "Required project directory does not exist. Checked exactly:\n",
      normalizePath(project_dir, winslash = "/", mustWork = FALSE),
      "\nV295 will not guess OneDrive, another drive, or another hierarchy."
    )
  )
}

if (!dir.exists(stage2_dir)) {
  stop(
    paste0(
      "Required STEP2 SAVE archive directory does not exist. Checked exactly:\n",
      normalizePath(stage2_dir, winslash = "/", mustWork = FALSE),
      "\nV295 will not guess another folder hierarchy."
    )
  )
}

if (!file.exists(input_rds)) {
  stop(
    paste0(
      "Required model input does not exist. Checked exactly:\n",
      normalizePath(input_rds, winslash = "/", mustWork = FALSE)
    )
  )
}

dir.create(output_root, recursive = TRUE, showWarnings = FALSE)

if (!dir.exists(output_root)) {
  stop(
    paste0(
      "Could not create output directory:\n",
      normalizePath(output_root, winslash = "/", mustWork = FALSE)
    )
  )
}


# ============================================================
# 3. Windows CPU/memory helper
# ============================================================

powershell_exe <- Sys.which("powershell")
if (!nzchar(powershell_exe)) powershell_exe <- Sys.which("pwsh")
if (!nzchar(powershell_exe)) stop("PowerShell/pwsh was not found.")

windows_helper <- file.path(output_root, "V295_windows_helper.ps1")

helper_lines <- c(
  'param(',
  '  [Parameter(Mandatory=$true)][ValidateSet("topology","affinity","memory","process_memory")][string]$Mode,',
  '  [int]$TargetPid = 0,',
  '  [int]$LogicalIndex = -1,',
  '  [string]$PriorityClass = "AboveNormal",',
  '  [string]$PidList = ""',
  ')',
  '$ErrorActionPreference = "Stop"',
  '',
  'if ($Mode -eq "memory") {',
  '  $os = Get-CimInstance Win32_OperatingSystem',
  '  [PSCustomObject]@{',
  '    TotalPhysicalMB = [math]::Round($os.TotalVisibleMemorySize / 1024.0, 3)',
  '    AvailablePhysicalMB = [math]::Round($os.FreePhysicalMemory / 1024.0, 3)',
  '    CommitLimitMB = [math]::Round($os.TotalVirtualMemorySize / 1024.0, 3)',
  '    CommitUsedMB = [math]::Round(($os.TotalVirtualMemorySize - $os.FreeVirtualMemory) / 1024.0, 3)',
  '  } | ConvertTo-Csv -NoTypeInformation',
  '  exit 0',
  '}',
  '',
  'if ($Mode -eq "process_memory") {',
  '  $ids = @()',
  '  if ($PidList) {',
  '    $ids = $PidList.Split(",") | ForEach-Object { [int]$_.Trim() }',
  '  }',
  '  $procs = @()',
  '  foreach ($id in $ids) {',
  '    $p = Get-Process -Id $id -ErrorAction SilentlyContinue',
  '    if ($null -ne $p) { $procs += $p }',
  '  }',
  '  $ws = ($procs | Measure-Object -Property WorkingSet64 -Sum).Sum',
  '  $pm = ($procs | Measure-Object -Property PrivateMemorySize64 -Sum).Sum',
  '  if ($null -eq $ws) { $ws = 0 }',
  '  if ($null -eq $pm) { $pm = 0 }',
  '  [PSCustomObject]@{',
  '    ProcessCount = $procs.Count',
  '    WorkingSetMB = [math]::Round($ws / 1MB, 3)',
  '    PrivateMemoryMB = [math]::Round($pm / 1MB, 3)',
  '  } | ConvertTo-Csv -NoTypeInformation',
  '  exit 0',
  '}',
  '',
  'if ($Mode -eq "affinity") {',
  '  if ($TargetPid -le 0) { throw "TargetPid must be positive." }',
  '  if ($LogicalIndex -lt 0 -or $LogicalIndex -gt 62) { throw "LogicalIndex must be 0-62." }',
  '  $p = Get-Process -Id $TargetPid',
  '  $mask = ([Int64]1 -shl $LogicalIndex)',
  '  $p.ProcessorAffinity = [IntPtr]$mask',
  '  if ($PriorityClass -and $PriorityClass -ne "Normal") {',
  '    $p.PriorityClass = $PriorityClass',
  '  }',
  '  [PSCustomObject]@{',
  '    PID = $TargetPid',
  '    LogicalIndex = $LogicalIndex',
  '    PriorityClass = $p.PriorityClass.ToString()',
  '    Status = "OK"',
  '  } | ConvertTo-Csv -NoTypeInformation',
  '  exit 0',
  '}',
  '',
  '$src = @"',
  'using System;',
  'using System.Runtime.InteropServices;',
  'public static class V295CpuSetNative {',
  '  [DllImport("kernel32.dll", SetLastError=true)]',
  '  public static extern bool GetSystemCpuSetInformation(',
  '    IntPtr Information, uint BufferLength, out uint ReturnedLength,',
  '    IntPtr Process, uint Flags);',
  '}',
  '"@',
  'Add-Type -TypeDefinition $src',
  '$needed = [uint32]0',
  '[void][V295CpuSetNative]::GetSystemCpuSetInformation([IntPtr]::Zero, 0, [ref]$needed, [IntPtr]::Zero, 0)',
  'if ($needed -le 0) { throw "GetSystemCpuSetInformation returned zero bytes." }',
  '$buf = [Runtime.InteropServices.Marshal]::AllocHGlobal([int]$needed)',
  'try {',
  '  $returned = [uint32]0',
  '  $ok = [V295CpuSetNative]::GetSystemCpuSetInformation($buf, $needed, [ref]$returned, [IntPtr]::Zero, 0)',
  '  if (-not $ok) { throw "GetSystemCpuSetInformation failed." }',
  '  $cpuName = (Get-CimInstance Win32_Processor | Select-Object -First 1 -ExpandProperty Name)',
  '  $offset = 0',
  '  $rows = @()',
  '  while ($offset -lt $returned) {',
  '    $ptr = [IntPtr]::Add($buf, $offset)',
  '    $size = [Runtime.InteropServices.Marshal]::ReadInt32($ptr, 0)',
  '    $type = [Runtime.InteropServices.Marshal]::ReadInt32($ptr, 4)',
  '    if ($size -le 0) { throw "Invalid CPU-set structure size." }',
  '    if ($type -eq 0) {',
  '      $rows += [PSCustomObject]@{',
  '        ProcessorName = $cpuName',
  '        Group = [int][Runtime.InteropServices.Marshal]::ReadInt16($ptr, 12)',
  '        LogicalProcessorIndex = [int][Runtime.InteropServices.Marshal]::ReadByte($ptr, 14)',
  '        CoreIndex = [int][Runtime.InteropServices.Marshal]::ReadByte($ptr, 15)',
  '        LastLevelCacheIndex = [int][Runtime.InteropServices.Marshal]::ReadByte($ptr, 16)',
  '        EfficiencyClass = [int][Runtime.InteropServices.Marshal]::ReadByte($ptr, 18)',
  '        SchedulingClass = [int][Runtime.InteropServices.Marshal]::ReadByte($ptr, 20)',
  '      }',
  '    }',
  '    $offset += $size',
  '  }',
  '  $rows | Sort-Object Group,LogicalProcessorIndex | ConvertTo-Csv -NoTypeInformation',
  '} finally {',
  '  [Runtime.InteropServices.Marshal]::FreeHGlobal($buf)',
  '}'
)

writeLines(helper_lines, windows_helper, useBytes = TRUE)

read_ps_csv <- function(lines) {
  if (length(lines) == 0L) return(data.table())
  suppressWarnings(
    data.table::fread(
      text = paste(lines, collapse = "\n"),
      showProgress = FALSE
    )
  )
}

run_ps_helper <- function(
  mode,
  target_pid = 0L,
  logical_index = -1L,
  priority_class = PROCESS_PRIORITY_CLASS,
  pid_list = integer(0),
  helper_path = windows_helper,
  powershell_path = powershell_exe
) {
  args <- c(
    "-NoProfile",
    "-ExecutionPolicy", "Bypass",
    "-File", shQuote(helper_path),
    "-Mode", mode
  )
  if (identical(mode, "affinity")) {
    args <- c(
      args,
      "-TargetPid", as.character(as.integer(target_pid)),
      "-LogicalIndex", as.character(as.integer(logical_index)),
      "-PriorityClass", priority_class
    )
  }
  if (identical(mode, "process_memory")) {
    args <- c(
      args,
      "-PidList", paste(as.integer(pid_list), collapse = ",")
    )
  }
  out <- system2(
    powershell_path,
    args = args,
    stdout = TRUE,
    stderr = TRUE
  )
  status <- attr(out, "status")
  if (!is.null(status) && status != 0L) {
    stop(
      paste0(
        "PowerShell helper failed in mode '", mode, "':\n",
        paste(out, collapse = "\n")
      )
    )
  }
  out
}

query_memory <- function() {
  x <- read_ps_csv(run_ps_helper("memory"))
  required <- c(
    "TotalPhysicalMB",
    "AvailablePhysicalMB",
    "CommitLimitMB",
    "CommitUsedMB"
  )
  if (nrow(x) != 1L || !all(required %in% names(x))) {
    stop("Windows memory query returned an unexpected schema.")
  }
  for (nm in required) x[[nm]] <- as.numeric(x[[nm]])
  x
}

query_process_memory <- function(pids) {
  x <- read_ps_csv(run_ps_helper("process_memory", pid_list = pids))
  required <- c("ProcessCount", "WorkingSetMB", "PrivateMemoryMB")
  if (nrow(x) != 1L || !all(required %in% names(x))) {
    stop("Windows process-memory query returned an unexpected schema.")
  }
  for (nm in required) x[[nm]] <- as.numeric(x[[nm]])
  x
}

apply_affinity <- function(pid, logical_index) {
  x <- read_ps_csv(
    run_ps_helper(
      "affinity",
      target_pid = pid,
      logical_index = logical_index
    )
  )
  if (
    nrow(x) != 1L ||
    !"Status" %in% names(x) ||
    !identical(as.character(x$Status[1L]), "OK")
  ) {
    stop(
      paste0(
        "Failed to pin PID ", pid,
        " to logical CPU ", logical_index, "."
      )
    )
  }
  invisible(TRUE)
}


# ============================================================
# 4. Detect exact machine and build CPU lanes
# ============================================================

topology_raw <- read_ps_csv(run_ps_helper("topology"))

topology_required <- c(
  "ProcessorName",
  "Group",
  "LogicalProcessorIndex",
  "CoreIndex",
  "LastLevelCacheIndex",
  "EfficiencyClass",
  "SchedulingClass"
)

if (nrow(topology_raw) == 0L || !all(topology_required %in% names(topology_raw))) {
  stop("Windows CPU topology query returned an unexpected schema.")
}

for (nm in setdiff(topology_required, "ProcessorName")) {
  topology_raw[[nm]] <- as.integer(topology_raw[[nm]])
}

cpu_name <- trimws(as.character(topology_raw$ProcessorName[1L]))

if (!grepl("AMD Ryzen 7 9800X3D", cpu_name, fixed = TRUE)) {
  stop(
    paste0(
      "This desktop script requires the frozen AMD Ryzen 7 9800X3D.\n",
      "Windows reported exactly:\n", cpu_name,
      "\nStop here rather than run the desktop allocation on another CPU."
    )
  )
}

MACHINE_ROLE <- "desktop"
EXPECTED_PHYSICAL_CORES <- 8L
EXPECTED_LOGICAL_CPUS <- 16L
STEP4_G_VALUES <- c(3L, 4L, 5L)
STEP4_PRIORITY_G <- c(5L, 4L, 3L)

task_table <- rbindlist(
  lapply(
    STEP4_PRIORITY_G,
    function(g_value) {
      data.table(
        requested_g = rep(as.integer(g_value), 5L),
        start_id = 1:5
      )
    }
  ),
  use.names = TRUE
)

role_dir <- file.path(output_root, MACHINE_ROLE)
dir.create(role_dir, recursive = TRUE, showWarnings = FALSE)
if (!dir.exists(role_dir)) stop("Could not create machine-role output folder.")

pause_request_path <- file.path(role_dir, "V295_PAUSE_REQUEST.txt")
if (file.exists(pause_request_path)) {
  stop(
    paste0(
      "Pause request already exists:\n",
      normalizePath(pause_request_path, winslash = "/", mustWork = TRUE),
      "\nDelete it only when you intend to resume."
    )
  )
}

topology <- copy(topology_raw)
topology[, core_key := paste(Group, CoreIndex, sep = "_")]

cores <- topology[
  ,
  .(
    ProcessorName = ProcessorName[1L],
    Group = Group[1L],
    CoreIndex = CoreIndex[1L],
    LastLevelCacheIndex = LastLevelCacheIndex[1L],
    EfficiencyClass = max(EfficiencyClass),
    SchedulingClass = max(SchedulingClass),
    n_logical = .N,
    logical_indices = list(sort(LogicalProcessorIndex))
  ),
  by = core_key
]

if (nrow(cores) != EXPECTED_PHYSICAL_CORES || nrow(topology) != EXPECTED_LOGICAL_CPUS) {
  fwrite(topology, file.path(role_dir, "V295_cpu_topology_unexpected.csv"))
  stop(
    paste0(
      "CPU core/thread count does not match the frozen machine specification.\n",
      "Detected physical cores: ", nrow(cores),
      "; expected: ", EXPECTED_PHYSICAL_CORES, "\n",
      "Detected logical CPUs: ", nrow(topology),
      "; expected: ", EXPECTED_LOGICAL_CPUS
    )
  )
}

if (!all(cores$n_logical == 2L)) {
  stop("9800X3D topology did not expose exactly two logical CPUs per core.")
}
cores[, core_type := "P"]

slot_list <- list()
slot_index <- 0L

for (i in seq_len(nrow(cores))) {
  logicals <- sort(as.integer(unlist(cores$logical_indices[[i]])))
  primary <- logicals[1L]
  core_type <- cores$core_type[i]
  tier_rank <- switch(core_type, "P" = 1L, "E" = 2L, "LP-E" = 3L)

  slot_index <- slot_index + 1L
  slot_list[[slot_index]] <- data.table(
    Group = cores$Group[i],
    CoreIndex = cores$CoreIndex[i],
    LogicalProcessorIndex = primary,
    tier_name = core_type,
    tier_rank = tier_rank,
    is_smt = 0L
  )

  if (length(logicals) > 1L) {
    for (sibling in logicals[-1L]) {
      slot_index <- slot_index + 1L
      slot_list[[slot_index]] <- data.table(
        Group = cores$Group[i],
        CoreIndex = cores$CoreIndex[i],
        LogicalProcessorIndex = sibling,
        tier_name = "SMT",
        tier_rank = 4L,
        is_smt = 1L
      )
    }
  }
}

cpu_slots <- rbindlist(slot_list)
setorder(cpu_slots, tier_rank, CoreIndex, LogicalProcessorIndex)
cpu_slots[, slot_order := .I]

if (nrow(cpu_slots) != EXPECTED_LOGICAL_CPUS) {
  stop("Constructed CPU-slot count differs from detected logical CPU count.")
}

# Keep the original V295 six-worker lane layout so resumed checkpoints continue
# under the same process allocation. G3-G5 use only their actual number of
# class jobs; unused workers remain idle. Two lanes reserve 12 logical CPU slots
# for inner pools plus two distinct controller slots.
MAX_CLASS_WORKERS_PER_FIT <- 6L

single_lane_worker_slots <- cpu_slots[
  seq_len(MAX_CLASS_WORKERS_PER_FIT)
]

single_lane_controller_slot <- cpu_slots[
  MAX_CLASS_WORKERS_PER_FIT + 1L
]

heavy_slots <- cpu_slots[
  seq_len(2L * MAX_CLASS_WORKERS_PER_FIT)
]

controller_slots <- cpu_slots[
  (2L * MAX_CLASS_WORKERS_PER_FIT + 1L):
    (2L * MAX_CLASS_WORKERS_PER_FIT + 2L)
]

calibration_lane_slots <- list(
  as.integer(single_lane_worker_slots$LogicalProcessorIndex)
)

calibration_controller_slots <- as.integer(
  single_lane_controller_slot$LogicalProcessorIndex[1L]
)

fwrite(topology, file.path(role_dir, "V295_cpu_topology.csv"))
fwrite(cpu_slots, file.path(role_dir, "V295_cpu_slots.csv"))


# ============================================================
# 5. Validate the exact Stage-2 model input
# ============================================================

cat("Reading and validating archived V286 STEP2 model input...\n")

input_md5 <- unname(tools::md5sum(input_rds))
if (
  length(input_md5) != 1L ||
  is.na(input_md5) ||
  !identical(tolower(input_md5), EXPECTED_INPUT_MD5)
) {
  stop(
    paste0(
      "Archived V286 STEP2 model input MD5 mismatch.\n",
      "Expected: ", EXPECTED_INPUT_MD5, "\n",
      "Observed: ", ifelse(length(input_md5) == 1L, input_md5, "unavailable"),
      "\nThe script will not read or fit an unverified input file."
    )
  )
}

input_check <- as.data.table(readRDS(input_rds))

missing_cols <- setdiff(REQUIRED_COLS, names(input_check))
if (length(missing_cols) > 0L) {
  stop(
    paste0(
      "Archived V286 STEP2 model input is missing exact required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

extra_cols <- setdiff(names(input_check), REQUIRED_COLS)
if (length(extra_cols) > 0L || length(names(input_check)) != length(REQUIRED_COLS)) {
  stop(
    paste0(
      "Archived V286 STEP2 model input contains columns outside the frozen seven-column input: ",
      ifelse(length(extra_cols) > 0L, paste(extra_cols, collapse = ", "), "duplicated names")
    )
  )
}

input_check[, unit_id := as.character(unit_id)]
input_check[, time_model := as.numeric(time_model)]
for (v in Z_MARKERS) input_check[, (v) := as.numeric(get(v))]

if (uniqueN(input_check$unit_id) != EXPECTED_PATIENTS) {
  stop(
    paste0(
      "Unexpected patient count: ", uniqueN(input_check$unit_id),
      "; expected ", EXPECTED_PATIENTS, "."
    )
  )
}

if (nrow(input_check) != EXPECTED_LONG_ROWS) {
  stop(
    paste0(
      "Unexpected long-row count: ", nrow(input_check),
      "; expected ", EXPECTED_LONG_ROWS, "."
    )
  )
}

if (anyNA(input_check$unit_id)) stop("unit_id contains missing values.")
if (anyNA(input_check$time_model)) stop("time_model contains missing values.")

duplicates <- input_check[, .N, by = .(unit_id, time_model)][N > 1L]
if (nrow(duplicates) > 0L) stop("Duplicate unit_id + time_model rows detected.")

rows_per_patient <- input_check[, .N, by = unit_id]
if (!all(rows_per_patient$N == EXPECTED_TIME_POINTS)) {
  stop("Not every patient has exactly 12 long-format time rows.")
}

input_audit <- data.table(
  version = VERSION_TAG,
  input_file = basename(input_rds),
  input_md5 = input_md5,
  n_patients = uniqueN(input_check$unit_id),
  n_rows = nrow(input_check),
  n_time_points = uniqueN(input_check$time_model),
  missing_msi = sum(is.na(input_check$msi_z)),
  missing_sodium = sum(is.na(input_check$sodium_z)),
  missing_potassium = sum(is.na(input_check$potassium_z)),
  missing_chloride = sum(is.na(input_check$chloride_z)),
  missing_bicarbonate = sum(is.na(input_check$bicarbonate_z))
)

fwrite(input_audit, file.path(role_dir, "V295_input_audit.csv"))

rm(input_check, duplicates, rows_per_patient)
gc(full = TRUE)


# ============================================================
# 6. General safe-write and audit helpers
# ============================================================

timestamp_now <- function() format(Sys.time(), "%Y-%m-%d %H:%M:%S")

append_log <- function(path, text) {
  cat(
    paste0(timestamp_now(), " | ", text, "\n"),
    file = path,
    append = TRUE
  )
}

safe_save_rds <- function(object, path, compress = FALSE) {
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  if (file.exists(tmp)) file.remove(tmp)
  saveRDS(object, tmp, compress = compress)
  if (!file.exists(tmp) || is.na(file.info(tmp)$size) || file.info(tmp)$size <= 0) {
    stop(paste0("Temporary RDS write failed: ", tmp))
  }
  if (file.exists(path)) file.remove(path)
  if (!file.rename(tmp, path)) {
    if (!file.copy(tmp, path, overwrite = TRUE)) {
      stop(paste0("Could not finalize RDS: ", path))
    }
    file.remove(tmp)
  }
  invisible(path)
}

safe_fwrite <- function(x, path) {
  tmp <- paste0(path, ".tmp_", Sys.getpid())
  if (file.exists(tmp)) file.remove(tmp)
  data.table::fwrite(x, tmp)
  if (!file.exists(tmp) || is.na(file.info(tmp)$size) || file.info(tmp)$size <= 0) {
    stop(paste0("Temporary CSV write failed: ", tmp))
  }
  if (file.exists(path)) file.remove(path)
  if (!file.rename(tmp, path)) {
    if (!file.copy(tmp, path, overwrite = TRUE)) {
      stop(paste0("Could not finalize CSV: ", path))
    }
    file.remove(tmp)
  }
  invisible(path)
}

relative_output_path <- function(path, root = output_root) {
  root_norm <- normalizePath(root, winslash = "/", mustWork = TRUE)
  path_norm <- normalizePath(path, winslash = "/", mustWork = TRUE)
  prefix <- paste0(root_norm, "/")
  if (!startsWith(path_norm, prefix)) {
    stop("Output path is not inside the V295 output root.")
  }
  substring(path_norm, nchar(prefix) + 1L)
}

seed_for_task <- function(g, start_id) {
  as.integer(
    STAGE2_MASTER_SEED +
      as.integer(g) * 100003L +
      (as.integer(start_id) - 1L) * START_SEED_STEP
  )
}

task_suffix <- function(g, start_id) {
  paste0("G", as.integer(g), "_S", sprintf("%02d", as.integer(start_id)))
}

task_paths <- function(g, start_id, output_dir = role_dir) {
  suffix <- task_suffix(g, start_id)
  list(
    full_fit = file.path(output_dir, paste0(VERSION_TAG, "_refine_fit_", suffix, ".rds")),
    compact_fit = file.path(output_dir, paste0(VERSION_TAG, "_compact_fit_", suffix, ".rds")),
    summary = file.path(output_dir, paste0(VERSION_TAG, "_task_summary_", suffix, ".csv")),
    class_quality = file.path(output_dir, paste0(VERSION_TAG, "_class_quality_", suffix, ".csv")),
    patient_class = file.path(output_dir, paste0(VERSION_TAG, "_patient_class_", suffix, ".csv")),
    trajectory = file.path(output_dir, paste0(VERSION_TAG, "_trajectory_", suffix, ".csv")),
    beta = file.path(output_dir, paste0(VERSION_TAG, "_beta_", suffix, ".csv")),
    sigma = file.path(output_dir, paste0(VERSION_TAG, "_sigma_", suffix, ".csv")),
    em_history = file.path(output_dir, paste0(VERSION_TAG, "_EM_history_", suffix, ".csv")),
    identity = file.path(output_dir, paste0(VERSION_TAG, "_task_identity_", suffix, ".rds")),
    checkpoint_a = file.path(output_dir, paste0(VERSION_TAG, "_checkpoint_A_", suffix, ".rds")),
    checkpoint_b = file.path(output_dir, paste0(VERSION_TAG, "_checkpoint_B_", suffix, ".rds")),
    running = file.path(output_dir, paste0(VERSION_TAG, "_RUNNING_", suffix, ".txt")),
    done = file.path(output_dir, paste0(VERSION_TAG, "_DONE_", suffix, ".txt")),
    error = file.path(output_dir, paste0(VERSION_TAG, "_ERROR_", suffix, ".txt")),
    stage_log = file.path(output_dir, paste0(VERSION_TAG, "_stage_", suffix, ".log")),
    resource_log = file.path(output_dir, paste0(VERSION_TAG, "_resource_", suffix, ".csv")),
    resume_audit = file.path(output_dir, paste0(VERSION_TAG, "_resume_audit_", suffix, ".csv"))
  )
}

task_is_complete <- function(g, start_id, output_dir = role_dir) {
  p <- task_paths(g, start_id, output_dir)
  if (!file.exists(p$done) || !file.exists(p$summary) || !file.exists(p$full_fit)) {
    return(FALSE)
  }
  s <- tryCatch(fread(p$summary, showProgress = FALSE), error = function(e) NULL)
  if (is.null(s) || nrow(s) != 1L || !"status" %in% names(s)) return(FALSE)
  as.character(s$status[1L]) %in% c(
    "CONVERGED_VALID_FULL_G",
    "CONVERGED_EMPTY_CLASS",
    "CONVERGED_STRUCTURE_MISMATCH",
    "NONCONVERGED_FULL_G",
    "NONCONVERGED_EMPTY_CLASS",
    "NONCONVERGED_STRUCTURE_MISMATCH"
  )
}


# ============================================================
# 7. Checkpoint helpers
# ============================================================

make_checkpoint_signature <- function(g, start_id, seed_used) {
  paste(
    VERSION_TAG,
    paste0("G=", as.integer(g)),
    paste0("start_id=", as.integer(start_id)),
    paste0("seed=", as.integer(seed_used)),
    paste0("input_md5=", input_md5),
    paste0("x=", paste(Z_MARKERS, collapse = ",")),
    paste0("d=", MODEL_DEGREE),
    paste0("tol=", format(REFINE_TOL, scientific = TRUE)),
    paste0("maxit=", NONBINDING_MAXIT),
    "scaling=0",
    "pruning=TRUE",
    "delete.empty=FALSE",
    "nstart=1",
    sep = "|"
  )
}

make_checkpoint_state <- function(
  signature,
  requested_g,
  start_id,
  seed_used,
  iteration,
  count2,
  ll_old,
  p0,
  beta,
  Sigma,
  postd,
  Z,
  em_history,
  peak_working_set_mb,
  peak_private_mb
) {
  list(
    state_version = CHECKPOINT_STATE_VERSION,
    signature = as.character(signature),
    requested_g = as.integer(requested_g),
    start_id = as.integer(start_id),
    seed = as.integer(seed_used),
    iteration = as.integer(iteration),
    count2 = as.integer(count2),
    ll_old = as.numeric(ll_old),
    p0 = p0,
    beta = beta,
    Sigma = Sigma,
    postd = postd,
    Z = Z,
    em_history = as.data.frame(em_history),
    peak_working_set_mb = as.numeric(peak_working_set_mb),
    peak_private_mb = as.numeric(peak_private_mb),
    saved_at = timestamp_now()
  )
}

validate_checkpoint_state <- function(state, signature, g, start_id) {
  required <- c(
    "state_version", "signature", "requested_g", "start_id", "seed",
    "iteration", "count2", "ll_old", "p0", "beta", "Sigma",
    "postd", "Z", "em_history", "peak_working_set_mb",
    "peak_private_mb", "saved_at"
  )
  if (!is.list(state)) return(list(valid = FALSE, reason = "not a list"))
  missing <- setdiff(required, names(state))
  if (length(missing) > 0L) {
    return(list(valid = FALSE, reason = paste("missing", paste(missing, collapse = ","))))
  }
  if (as.integer(state$state_version) != CHECKPOINT_STATE_VERSION) {
    return(list(valid = FALSE, reason = "state version mismatch"))
  }
  if (!identical(as.character(state$signature), as.character(signature))) {
    return(list(valid = FALSE, reason = "signature mismatch"))
  }
  if (as.integer(state$requested_g) != as.integer(g)) {
    return(list(valid = FALSE, reason = "G mismatch"))
  }
  if (as.integer(state$start_id) != as.integer(start_id)) {
    return(list(valid = FALSE, reason = "start_id mismatch"))
  }
  if (
    length(state$p0) != g ||
    length(state$beta) != g ||
    length(state$Sigma) != g ||
    ncol(as.matrix(state$postd)) != g
  ) {
    return(list(valid = FALSE, reason = "group dimension mismatch"))
  }
  if (!is.finite(state$ll_old) || !is.finite(state$iteration)) {
    return(list(valid = FALSE, reason = "invalid numeric state"))
  }
  list(valid = TRUE, reason = "OK")
}

read_checkpoint_slot <- function(path, signature, g, start_id) {
  if (!file.exists(path)) {
    return(list(exists = FALSE, valid = FALSE, state = NULL, reason = "missing", path = path))
  }
  state <- tryCatch(readRDS(path), error = function(e) e)
  if (inherits(state, "error")) {
    return(
      list(
        exists = TRUE,
        valid = FALSE,
        state = NULL,
        reason = paste0("read error: ", conditionMessage(state)),
        path = path
      )
    )
  }
  verdict <- validate_checkpoint_state(state, signature, g, start_id)
  list(
    exists = TRUE,
    valid = isTRUE(verdict$valid),
    state = if (isTRUE(verdict$valid)) state else NULL,
    reason = verdict$reason,
    path = path
  )
}

load_latest_checkpoint <- function(paths, signature, g, start_id) {
  a <- read_checkpoint_slot(paths$checkpoint_a, signature, g, start_id)
  b <- read_checkpoint_slot(paths$checkpoint_b, signature, g, start_id)
  valid <- Filter(function(x) isTRUE(x$valid), list(a, b))
  if (length(valid) == 0L) {
    if (!a$exists && !b$exists) return(NULL)
    stop(
      paste0(
        "Checkpoint file(s) exist but neither is valid for ",
        task_suffix(g, start_id), ".\n",
        "Slot A: ", a$reason, "\nSlot B: ", b$reason,
        "\nV295 will not silently restart this task."
      )
    )
  }
  if (length(valid) == 1L) return(valid[[1L]]$state)
  iterations <- vapply(valid, function(x) as.integer(x$state$iteration), integer(1))
  if (length(unique(iterations)) > 1L) return(valid[[which.max(iterations)]]$state)
  mtimes <- vapply(valid, function(x) as.numeric(file.info(x$path)$mtime), numeric(1))
  valid[[which.max(mtimes)]]$state
}

save_checkpoint_dual <- function(state, paths) {
  target <- if (as.integer(state$iteration) %% 2L == 0L) {
    paths$checkpoint_a
  } else {
    paths$checkpoint_b
  }
  safe_save_rds(state, target, compress = FALSE)
  test <- readRDS(target)
  if (
    as.integer(test$iteration) != as.integer(state$iteration) ||
    !identical(as.character(test$signature), as.character(state$signature))
  ) {
    stop("Checkpoint verification failed after finalized write.")
  }
  invisible(target)
}


# ============================================================
# 8. Static model-data preparation
# ============================================================

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


# ============================================================
# 9. Class-parallel EM kernels executed by inner workers
# ============================================================

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


# ============================================================
# 10. Per-task memory gate and resource measurement
# ============================================================

read_observed_peak_mb <- function(output_dir = role_dir) {
  files <- list.files(
    output_dir,
    pattern = paste0("^", VERSION_TAG, "_task_summary_G[0-9]+_S[0-9]+\\.csv$"),
    full.names = TRUE
  )
  if (length(files) == 0L) return(NA_real_)
  vals <- unlist(
    lapply(
      files,
      function(path) {
        x <- tryCatch(fread(path, showProgress = FALSE), error = function(e) NULL)
        if (is.null(x) || !"peak_private_mb" %in% names(x)) return(NA_real_)
        suppressWarnings(as.numeric(x$peak_private_mb[1L]))
      }
    )
  )
  vals <- vals[is.finite(vals) & vals > 0]
  if (length(vals) == 0L) NA_real_ else max(vals)
}

wait_for_memory_gate <- function(g, start_id, stage_log) {
  repeat {
    mem <- query_memory()
    required_available_mb <- max(
      MEMORY_MIN_AVAILABLE_GB * 1024,
      MEMORY_MIN_AVAILABLE_FRACTION * mem$TotalPhysicalMB
    )
    observed_peak <- read_observed_peak_mb(role_dir)
    if (is.finite(observed_peak)) {
      required_available_mb <- max(
        required_available_mb,
        observed_peak * MEMORY_PEAK_SAFETY_FACTOR
      )
    }
    commit_fraction <- if (mem$CommitLimitMB > 0) {
      mem$CommitUsedMB / mem$CommitLimitMB
    } else {
      NA_real_
    }
    ok <- mem$AvailablePhysicalMB >= required_available_mb &&
      (is.na(commit_fraction) || commit_fraction <= MEMORY_MAX_COMMIT_FRACTION)
    if (isTRUE(ok)) return(invisible(mem))

    append_log(
      stage_log,
      paste0(
        "memory gate waiting for ", task_suffix(g, start_id),
        "; available_GB=", round(mem$AvailablePhysicalMB / 1024, 2),
        "; required_GB=", round(required_available_mb / 1024, 2),
        "; commit_fraction=", round(commit_fraction, 4)
      )
    )
    Sys.sleep(MEMORY_POLL_SECONDS)
  }
}

measure_task_resources <- function(controller_pid, inner_pids) {
  process_mem <- query_process_memory(c(controller_pid, inner_pids))
  system_mem <- query_memory()
  data.table(
    recorded_at = timestamp_now(),
    controller_pid = as.integer(controller_pid),
    inner_processes_found = as.integer(process_mem$ProcessCount),
    working_set_mb = as.numeric(process_mem$WorkingSetMB),
    private_memory_mb = as.numeric(process_mem$PrivateMemoryMB),
    available_physical_mb = as.numeric(system_mem$AvailablePhysicalMB),
    commit_used_mb = as.numeric(system_mem$CommitUsedMB),
    commit_limit_mb = as.numeric(system_mem$CommitLimitMB)
  )
}


# ============================================================
# 11. Parallel checkpoint-enabled one-start GBMTM engine
# ============================================================

run_parallel_one_start <- function(
  g,
  start_id,
  seed_used,
  paths,
  signature,
  inner_cluster,
  class_chunks,
  inner_pids,
  static,
  pause_path
) {
  resume_state <- load_latest_checkpoint(paths, signature, g, start_id)
  resumed <- !is.null(resume_state)
  resumed_from_iteration <- if (resumed) as.integer(resume_state$iteration) else 0L

  parallel::clusterEvalQ(
    inner_cluster,
    {
      if (exists("V295_DATAI_CACHE", envir = .GlobalEnv, inherits = FALSE)) {
        rm("V295_DATAI_CACHE", envir = .GlobalEnv)
      }
      NULL
    }
  )

  peak_working_set_mb <- 0
  peak_private_mb <- 0
  resource_history <- data.table()

  record_resources <- function(iteration, phase) {
    r <- measure_task_resources(Sys.getpid(), inner_pids)
    r[, iteration := as.integer(iteration)]
    r[, phase := as.character(phase)]
    resource_history <<- rbind(resource_history, r, fill = TRUE)
    peak_working_set_mb <<- max(peak_working_set_mb, r$working_set_mb, na.rm = TRUE)
    peak_private_mb <<- max(peak_private_mb, r$private_memory_mb, na.rm = TRUE)
    safe_fwrite(resource_history, paths$resource_log)
  }

  postp <- matrix(nrow = static$n, ncol = g)
  rownames(postp) <- static$cou

  if (!resumed) {
    set.seed(seed_used)
    Z <- sample(rep(seq_len(g), length = static$n))
    names(Z) <- static$cou
    p0 <- prop.table(table(factor(Z, levels = seq_len(g))))

    init_chunks <- parallel::clusterMap(
      inner_cluster,
      fun = inner_init_chunk,
      class_ids = class_chunks,
      MoreArgs = list(Z = Z, p0 = p0),
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE,
      RECYCLE = FALSE
    )
    init <- flatten_class_results(init_chunks, g, static$n)
    beta <- init$beta
    Sigma <- init$Sigma
    postd <- init$postd
    rownames(postd) <- static$cou

    ll_old <- sum(log(apply(exp(postd), 1L, sum)))
    count <- 0L
    count2 <- 0L
    em_history <- data.table(
      iteration = 0L,
      logLik = as.numeric(ll_old),
      delta_logLik = NA_real_,
      iteration_seconds = NA_real_,
      stop_reason = "INITIALIZATION",
      checkpoint_saved = 1L,
      recorded_at = timestamp_now()
    )

    record_resources(0L, "initialization")

    state <- make_checkpoint_state(
      signature, g, start_id, seed_used, count, count2, ll_old,
      p0, beta, Sigma, postd, Z, em_history,
      peak_working_set_mb, peak_private_mb
    )
    save_checkpoint_dual(state, paths)
    safe_fwrite(em_history, paths$em_history)
    append_log(paths$stage_log, paste0("fresh initialization; ll=", format(ll_old, digits = 16)))
  } else {
    count <- as.integer(resume_state$iteration)
    count2 <- as.integer(resume_state$count2)
    ll_old <- as.numeric(resume_state$ll_old)
    p0 <- resume_state$p0
    beta <- resume_state$beta
    Sigma <- resume_state$Sigma
    postd <- as.matrix(resume_state$postd)
    Z <- resume_state$Z
    names(Z) <- static$cou
    rownames(postd) <- static$cou
    em_history <- as.data.table(resume_state$em_history)
    peak_working_set_mb <- as.numeric(resume_state$peak_working_set_mb)
    peak_private_mb <- as.numeric(resume_state$peak_private_mb)
    append_log(
      paths$stage_log,
      paste0("resumed from iteration ", count, "; ll=", format(ll_old, digits = 16))
    )
  }

  if (file.exists(pause_path)) {
    return(
      list(
        paused = TRUE,
        fit = NULL,
        resumed = resumed,
        resumed_from_iteration = resumed_from_iteration,
        final_iteration = count,
        peak_working_set_mb = peak_working_set_mb,
        peak_private_mb = peak_private_mb,
        stop_reason = "PAUSE_AFTER_CHECKPOINT"
      )
    )
  }

  fine <- 0L
  stop_reason <- ""
  p0_old <- p0
  beta_old <- beta
  Sigma_old <- Sigma

  while (fine == 0L) {
    iteration_start <- Sys.time()
    count <- count + 1L

    for (i in seq_len(nrow(postd))) {
      postp[i, ] <- exp(postd[i, ]) / sum(exp(postd[i, ]))
    }

    Z_old <- Z
    Z <- apply(postp, 1L, function(x) which.max(x))
    p0 <- colMeans(postp)

    mstep_chunks <- parallel::clusterMap(
      inner_cluster,
      fun = inner_mstep_chunk,
      class_ids = class_chunks,
      MoreArgs = list(postp = postp, beta = beta, Sigma = Sigma),
      SIMPLIFY = FALSE,
      USE.NAMES = FALSE,
      RECYCLE = FALSE
    )

    mstep <- flatten_class_results(mstep_chunks, g, static$n)
    beta <- mstep$beta
    Sigma <- mstep$Sigma
    postd <- mstep$postd
    rownames(postd) <- static$cou

    ll_new <- sum(log(apply(exp(postd), 1L, sum)))
    raw_delta_ll <- ll_new - ll_old
    reasons <- character(0)

    if (ll_new >= ll_old) {
      if (ll_new - ll_old < REFINE_TOL) {
        fine <- 1L
        reasons <- c(reasons, "DELTA_LL_BELOW_TOL")
      }
    } else {
      p0 <- p0_old
      beta <- beta_old
      Sigma <- Sigma_old
      ll_new <- ll_old
      count <- count - 1L
      fine <- 1L
      reasons <- c(reasons, "LL_DECREASE_ROLLBACK")
    }

    if (count >= NONBINDING_MAXIT) {
      fine <- 1L
      reasons <- c(reasons, "MAXIT")
    }

    if (sum(Z_old - Z) == 0) {
      count2 <- count2 + 1L
      if (ll_new - ll_old < REFINE_TOL * sqrt(count2)) {
        fine <- 1L
        reasons <- c(reasons, "STABLE_ASSIGNMENT_RULE")
      }
    } else {
      count2 <- 0L
    }

    iteration_seconds <- as.numeric(
      difftime(Sys.time(), iteration_start, units = "secs")
    )

    stop_reason_iter <- if (length(reasons) == 0L) "" else paste(unique(reasons), collapse = "+")

    em_history <- rbind(
      em_history,
      data.table(
        iteration = as.integer(count),
        logLik = as.numeric(ll_new),
        delta_logLik = as.numeric(ll_new - ll_old),
        raw_delta_logLik = as.numeric(raw_delta_ll),
        iteration_seconds = as.numeric(iteration_seconds),
        stop_reason = stop_reason_iter,
        checkpoint_saved = 0L,
        recorded_at = timestamp_now()
      ),
      fill = TRUE
    )

    record_resources(count, "EM_iteration")

    if (fine == 0L) {
      ll_old <- ll_new
      p0_old <- p0
      beta_old <- beta
      Sigma_old <- Sigma

      if (count %% CHECKPOINT_EVERY_ACCEPTED_ITER == 0L) {
        em_history[.N, checkpoint_saved := 1L]
        state <- make_checkpoint_state(
          signature, g, start_id, seed_used, count, count2, ll_old,
          p0, beta, Sigma, postd, Z, em_history,
          peak_working_set_mb, peak_private_mb
        )
        saved <- save_checkpoint_dual(state, paths)
        append_log(
          paths$stage_log,
          paste0(
            "accepted iteration ", count,
            "; ll=", format(ll_old, digits = 16),
            "; checkpoint=", basename(saved)
          )
        )
      }
    }

    safe_fwrite(em_history, paths$em_history)

    if (fine == 0L && file.exists(pause_path)) {
      return(
        list(
          paused = TRUE,
          fit = NULL,
          resumed = resumed,
          resumed_from_iteration = resumed_from_iteration,
          final_iteration = count,
          peak_working_set_mb = peak_working_set_mb,
          peak_private_mb = peak_private_mb,
          stop_reason = "PAUSE_AFTER_CHECKPOINT"
        )
      )
    }

    if (fine == 1L) stop_reason <- stop_reason_iter
  }

  for (i in seq_len(nrow(postd))) {
    postp[i, ] <- exp(postd[i, ]) / sum(exp(postd[i, ]))
  }
  Z <- apply(postp, 1L, function(x) which.max(x))
  names(Z) <- static$cou

  assL <- vector("list", length = g)
  for (i in seq_len(g)) {
    assL[[i]] <- sort(names(which(Z == i)))
  }

  cache_chunks <- parallel::clusterEvalQ(
    inner_cluster,
    {
      if (exists("V295_DATAI_CACHE", envir = .GlobalEnv, inherits = FALSE)) {
        get("V295_DATAI_CACHE", envir = .GlobalEnv, inherits = FALSE)
      } else {
        list()
      }
    }
  )

  dataI_list <- vector("list", g)
  for (chunk in cache_chunks) {
    if (length(chunk) > 0L) {
      for (nm in names(chunk)) dataI_list[[as.integer(nm)]] <- chunk[[nm]]
    }
  }

  if (static$ismiss && any(vapply(dataI_list, is.null, logical(1)))) {
    stop("Final class-specific imputed data cache is incomplete.")
  }

  if (static$ismiss) {
    auxdat <- dataI_list[[1L]][, static$x_names, drop = FALSE] * p0[1L]
    if (g > 1L) {
      for (i in 2:g) {
        auxdat <- auxdat + dataI_list[[i]][, static$x_names, drop = FALSE] * p0[i]
      }
    }
    newdat <- data.frame(
      static$data[, c(static$unit_name, static$time_name)],
      auxdat
    )
  } else {
    newdat <- static$data[, c(static$unit_name, static$time_name, static$x_names)]
  }

  pruneFun <- getFromNamespace("pruneFun", "gbmt")
  nparCalc <- getFromNamespace("nparCalc", "gbmt")
  icCalc <- getFromNamespace("icCalc", "gbmt")
  appaCalc <- getFromNamespace("appaCalc", "gbmt")

  mpruned <- pruneFun(beta = beta, posterior = postp, data = newdat)
  reg <- mpruned$reg
  beta <- mpruned$beta
  Sigma <- mpruned$S
  ll_final <- mpruned$logLik

  npar <- nparCalc(beta)
  ic <- icCalc(ll_final, npar = npar, ss = nrow(newdat))

  pred <- vector("list", length = g)
  names(pred) <- seq_len(g)
  for (i in seq_len(g)) {
    pred[[i]] <- data.frame(static$Xmat %*% beta[[i]])
    rownames(pred[[i]]) <- static$tnam
    colnames(pred[[i]]) <- static$x_names
  }

  names(p0) <- names(assL) <- colnames(postp) <-
    names(beta) <- names(Sigma) <- names(reg) <- seq_len(g)

  fit <- list(
    call = list(
      x.names = static$x_names,
      unit = static$unit_name,
      time = static$time_name,
      ng = g,
      d = MODEL_DEGREE,
      scaling = 0,
      pruning = TRUE,
      nstart = 1L,
      tol = REFINE_TOL,
      maxit = NONBINDING_MAXIT,
      delete.empty = FALSE
    ),
    prior = p0,
    beta = beta,
    Sigma = Sigma,
    posterior = round(postp, 4),
    Xmat = static$Xmat,
    fitted = pred,
    reg = reg,
    assign = Z,
    assign.list = assL,
    logLik = ll_final,
    npar = npar,
    ic = ic,
    appa = appaCalc(postp),
    data.orig = static$data[, c(static$unit_name, static$time_name, static$x_names)],
    data.norm = newdat,
    data.imputed = newdat,
    em = c(
      logLik = ll_final,
      n.iter = count,
      converged = ifelse(count >= NONBINDING_MAXIT, 0, 1)
    )
  )
  class(fit) <- "gbmt"

  record_resources(count, "final_pruning_and_fit_object")

  list(
    paused = FALSE,
    fit = fit,
    resumed = resumed,
    resumed_from_iteration = resumed_from_iteration,
    final_iteration = count,
    peak_working_set_mb = peak_working_set_mb,
    peak_private_mb = peak_private_mb,
    stop_reason = stop_reason,
    em_logLik_before_pruning = ll_new,
    final_logLik_after_pruning = ll_final
  )
}


# ============================================================
# 12. Result extraction
# ============================================================

safe_ic_value <- function(fit, target) {
  x <- fit$ic
  if (is.null(x) || is.null(names(x))) return(NA_real_)
  idx <- which(tolower(names(x)) == tolower(target))
  if (length(idx) == 0L) return(NA_real_)
  as.numeric(x[idx[1L]])
}

extract_class_quality <- function(fit, g, start_id) {
  post <- as.matrix(fit$posterior)
  assign <- as.integer(fit$assign)
  max_post <- apply(post, 1L, max, na.rm = TRUE)

  rows <- rbindlist(
    lapply(
      seq_len(ncol(post)),
      function(k) {
        idx <- which(assign == k)
        n_k <- length(idx)
        appa_recalculated <- if (n_k > 0L) {
          mean(post[idx, k], na.rm = TRUE)
        } else {
          NA_real_
        }
        appa_model <- if (length(fit$appa) >= k) {
          as.numeric(fit$appa[k])
        } else {
          appa_recalculated
        }
        prior_k <- as.numeric(fit$prior[k])
        occ <- if (
          is.finite(appa_model) && is.finite(prior_k) &&
          appa_model > 0 && appa_model < 1 &&
          prior_k > 0 && prior_k < 1
        ) {
          (appa_model / (1 - appa_model)) /
            (prior_k / (1 - prior_k))
        } else {
          NA_real_
        }
        data.table(
          version = VERSION_TAG,
          requested_g = g,
          start_id = start_id,
          class = k,
          n = n_k,
          pct = 100 * n_k / length(assign),
          prior = prior_k,
          APPA = appa_model,
          APPA_recalculated_from_rounded_posterior = appa_recalculated,
          OCC = occ,
          mean_max_posterior = if (n_k > 0L) mean(max_post[idx], na.rm = TRUE) else NA_real_,
          median_max_posterior = if (n_k > 0L) median(max_post[idx], na.rm = TRUE) else NA_real_,
          pct_max_posterior_ge_070 = if (n_k > 0L) 100 * mean(max_post[idx] >= 0.70) else NA_real_,
          pct_max_posterior_ge_090 = if (n_k > 0L) 100 * mean(max_post[idx] >= 0.90) else NA_real_
        )
      }
    )
  )
  rows
}

extract_patient_class <- function(fit, static, g, start_id) {
  post <- as.matrix(fit$posterior)
  assign <- as.integer(fit$assign)
  out <- data.table(
    version = VERSION_TAG,
    requested_g = g,
    start_id = start_id,
    unit_id = rownames(post),
    assigned_class = assign,
    max_posterior = apply(post, 1L, max, na.rm = TRUE),
    assigned_posterior = post[cbind(seq_len(nrow(post)), assign)]
  )
  for (k in seq_len(ncol(post))) {
    out[, (paste0("p_class", k)) := post[, k]]
  }
  out
}

extract_trajectory <- function(fit, g, start_id) {
  rows <- list()
  idx <- 0L
  for (k in seq_len(g)) {
    m <- as.matrix(fit$fitted[[k]])
    for (v in seq_along(Z_MARKERS)) {
      variable <- Z_MARKERS[v]
      scaling_row <- SCALING_TABLE[SCALING_TABLE$variable == variable, ]
      for (t in seq_len(nrow(m))) {
        idx <- idx + 1L
        fitted_z <- as.numeric(m[t, variable])
        rows[[idx]] <- data.table(
          version = VERSION_TAG,
          requested_g = g,
          start_id = start_id,
          class = k,
          variable = variable,
          indicator = scaling_row$indicator,
          time_model = as.numeric(rownames(m)[t]),
          time_hour_center = 2 * t - 1,
          fitted_z = fitted_z,
          fitted_native = scaling_row$mean_native + scaling_row$sd_native * fitted_z
        )
      }
    }
  }
  rbindlist(rows)
}

extract_beta <- function(fit, g, start_id) {
  rows <- list()
  idx <- 0L
  for (k in seq_len(g)) {
    b <- as.matrix(fit$beta[[k]])
    for (i in seq_len(nrow(b))) {
      for (j in seq_len(ncol(b))) {
        idx <- idx + 1L
        rows[[idx]] <- data.table(
          version = VERSION_TAG,
          requested_g = g,
          start_id = start_id,
          class = k,
          term = rownames(b)[i],
          variable = colnames(b)[j],
          coefficient = as.numeric(b[i, j])
        )
      }
    }
  }
  rbindlist(rows)
}

extract_sigma <- function(fit, g, start_id) {
  rows <- list()
  idx <- 0L
  for (k in seq_len(g)) {
    s <- as.matrix(fit$Sigma[[k]])
    for (i in seq_len(nrow(s))) {
      for (j in seq_len(ncol(s))) {
        idx <- idx + 1L
        rows[[idx]] <- data.table(
          version = VERSION_TAG,
          requested_g = g,
          start_id = start_id,
          class = k,
          variable_row = rownames(s)[i],
          variable_col = colnames(s)[j],
          covariance = as.numeric(s[i, j])
        )
      }
    }
  }
  rbindlist(rows)
}


# ============================================================
# 12B. Within-G label alignment and start-consistency analysis
# ============================================================

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

read_start_outputs <- function(g, output_dir, element) {
  rows <- vector("list", 5L)
  for (start_id in 1:5) {
    path <- task_paths(g, start_id, output_dir)[[element]]
    if (!file.exists(path)) {
      stop(
        paste0(
          "Cannot analyze G", g, ": missing ", element,
          " for start ", start_id, "."
        )
      )
    }
    rows[[start_id]] <- fread(path, showProgress = FALSE)
  }
  rbindlist(rows, fill = TRUE)
}

analyze_completed_g <- function(g, output_dir = role_dir) {
  all_done <- all(
    vapply(
      1:5,
      function(start_id) task_is_complete(g, start_id, output_dir),
      logical(1)
    )
  )
  if (!all_done) return(invisible(FALSE))

  summaries <- read_start_outputs(g, output_dir, "summary")
  summaries[, start_id := as.integer(start_id)]
  setorder(summaries, start_id)

  benchmark <- as.data.table(STAGE2_BENCHMARK)[requested_g == g]
  if (nrow(benchmark) != 1L) stop("Stage-2 benchmark lookup failed.")
  summaries[
    ,
    `:=`(
      stage2_tol = benchmark$stage2_tol[1L],
      stage2_em_iterations = benchmark$stage2_em_iterations[1L],
      stage2_logLik = benchmark$stage2_logLik[1L],
      stage2_BIC = benchmark$stage2_BIC[1L],
      refined_logLik_minus_stage2 = final_logLik - benchmark$stage2_logLik[1L],
      refined_BIC_minus_stage2 = BIC - benchmark$stage2_BIC[1L]
    )
  ]

  eligible <- summaries[
    status == "CONVERGED_VALID_FULL_G" & is.finite(final_logLik)
  ]
  reference_rule <- "maximum logLik among CONVERGED_VALID_FULL_G starts"
  if (nrow(eligible) == 0L) {
    eligible <- summaries[structural_full_g == 1L & is.finite(final_logLik)]
    reference_rule <- "fallback: maximum logLik among structural full-G starts"
  }
  if (nrow(eligible) == 0L) {
    eligible <- summaries[is.finite(final_logLik)]
    reference_rule <- "fallback: maximum finite logLik"
  }
  if (nrow(eligible) == 0L) stop(paste0("No finite G", g, " fit is available."))

  reference_start <- as.integer(
    eligible[order(-final_logLik, start_id), start_id][1L]
  )

  trajectories <- read_start_outputs(g, output_dir, "trajectory")
  class_quality <- read_start_outputs(g, output_dir, "class_quality")
  betas <- read_start_outputs(g, output_dir, "beta")
  sigmas <- read_start_outputs(g, output_dir, "sigma")

  expected_trajectory_rows <- g * length(Z_MARKERS) * EXPECTED_TIME_POINTS
  trajectory_counts <- trajectories[, .N, by = start_id]
  if (
    nrow(trajectory_counts) != 5L ||
    any(trajectory_counts$N != expected_trajectory_rows)
  ) {
    stop(paste0("G", g, " trajectory files have unexpected dimensions."))
  }

  reference_trajectory <- trajectories[start_id == reference_start]
  label_maps <- vector("list", 5L)

  for (sid in 1:5) {
    current <- trajectories[start_id == sid]
    cost <- matrix(NA_real_, nrow = g, ncol = g)

    for (reference_class in seq_len(g)) {
      reference_values <- reference_trajectory[
        class == reference_class,
        .(variable, time_model, reference_z = fitted_z)
      ]
      setorder(reference_values, variable, time_model)

      for (current_class in seq_len(g)) {
        current_values <- current[
          class == current_class,
          .(variable, time_model, current_z = fitted_z)
        ]
        matched <- merge(
          reference_values,
          current_values,
          by = c("variable", "time_model"),
          all = FALSE,
          sort = TRUE
        )
        if (nrow(matched) != length(Z_MARKERS) * EXPECTED_TIME_POINTS) {
          stop("Trajectory keys were incomplete during label alignment.")
        }
        cost[reference_class, current_class] <- sqrt(
          mean((matched$reference_z - matched$current_z)^2)
        )
      }
    }

    reference_to_original <- if (sid == reference_start) {
      seq_len(g)
    } else {
      hungarian_min_cost(cost)
    }

    label_maps[[sid]] <- data.table(
      version = VERSION_TAG,
      requested_g = g,
      reference_start_id = reference_start,
      start_id = sid,
      original_class = as.integer(reference_to_original),
      aligned_class = seq_len(g),
      trajectory_rmse_to_reference = cost[
        cbind(seq_len(g), reference_to_original)
      ],
      reference_selection_rule = reference_rule
    )
  }

  label_map <- rbindlist(label_maps)

  align_table <- function(x) {
    x <- copy(x)
    x[, original_class := as.integer(class)]
    x[, class := NULL]
    x <- merge(
      x,
      label_map[
        ,
        .(start_id, original_class, aligned_class)
      ],
      by = c("start_id", "original_class"),
      all.x = TRUE,
      sort = FALSE
    )
    if (anyNA(x$aligned_class)) stop("A class could not be label-aligned.")
    setnames(x, "aligned_class", "class")
    setcolorder(
      x,
      c(
        intersect(c("version", "requested_g", "start_id"), names(x)),
        "original_class", "class",
        setdiff(
          names(x),
          c("version", "requested_g", "start_id", "original_class", "class")
        )
      )
    )
    x
  }

  trajectories_aligned <- align_table(trajectories)
  class_aligned <- align_table(class_quality)
  beta_aligned <- align_table(betas)
  sigma_aligned <- align_table(sigmas)

  reference_patient <- fread(
    task_paths(g, reference_start, output_dir)$patient_class,
    showProgress = FALSE
  )
  reference_patient[, unit_id := as.character(unit_id)]
  setorder(reference_patient, unit_id)
  reference_assign <- as.integer(reference_patient$assigned_class)
  posterior_names <- paste0("p_class", seq_len(g))
  if (!all(posterior_names %in% names(reference_patient))) {
    stop("Reference posterior columns are incomplete.")
  }
  reference_posterior <- as.matrix(reference_patient[, ..posterior_names])

  consistency_rows <- vector("list", 5L)
  variable_rows <- list()
  variable_index <- 0L

  reference_class_quality <- class_aligned[start_id == reference_start]
  setorder(reference_class_quality, class)

  for (sid in 1:5) {
    map <- label_map[start_id == sid]
    original_to_aligned <- integer(g)
    original_to_aligned[map$original_class] <- map$aligned_class

    patient <- fread(
      task_paths(g, sid, output_dir)$patient_class,
      showProgress = FALSE
    )
    patient[, unit_id := as.character(unit_id)]
    if (!all(posterior_names %in% names(patient))) {
      stop(paste0("Posterior columns are incomplete for G", g, " start ", sid, "."))
    }
    order_to_reference <- match(reference_patient$unit_id, patient$unit_id)
    if (anyNA(order_to_reference)) stop("Patient IDs differ across starts.")
    patient <- patient[order_to_reference]

    aligned_assign <- original_to_aligned[as.integer(patient$assigned_class)]
    current_posterior <- as.matrix(patient[, ..posterior_names])
    aligned_posterior <- matrix(NA_real_, nrow(current_posterior), g)
    for (original_class in seq_len(g)) {
      aligned_posterior[, original_to_aligned[original_class]] <-
        current_posterior[, original_class]
    }

    current_class_quality <- class_aligned[start_id == sid]
    setorder(current_class_quality, class)
    class_compare <- merge(
      reference_class_quality[
        ,
        .(
          class,
          reference_pct = pct,
          reference_prior = prior,
          reference_APPA = APPA
        )
      ],
      current_class_quality[
        ,
        .(
          class,
          current_pct = pct,
          current_prior = prior,
          current_APPA = APPA
        )
      ],
      by = "class",
      all = TRUE,
      sort = TRUE
    )

    current_trajectory <- trajectories_aligned[start_id == sid]
    reference_aligned_trajectory <- trajectories_aligned[
      start_id == reference_start
    ]
    trajectory_compare <- merge(
      reference_aligned_trajectory[
        ,
        .(class, variable, time_model, reference_z = fitted_z)
      ],
      current_trajectory[
        ,
        .(class, variable, time_model, current_z = fitted_z)
      ],
      by = c("class", "variable", "time_model"),
      all = FALSE,
      sort = TRUE
    )
    if (nrow(trajectory_compare) != expected_trajectory_rows) {
      stop("Aligned trajectory dimensions are incomplete.")
    }

    trajectory_by_variable <- trajectory_compare[
      ,
      .(
        trajectory_RMSE = sqrt(mean((reference_z - current_z)^2)),
        trajectory_MAE = mean(abs(reference_z - current_z)),
        trajectory_max_abs = max(abs(reference_z - current_z))
      ),
      by = variable
    ]
    trajectory_by_variable[
      ,
      `:=`(
        version = VERSION_TAG,
        requested_g = g,
        reference_start_id = reference_start,
        start_id = sid
      )
    ]
    setcolorder(
      trajectory_by_variable,
      c(
        "version", "requested_g", "reference_start_id", "start_id",
        "variable", "trajectory_RMSE", "trajectory_MAE",
        "trajectory_max_abs"
      )
    )
    variable_index <- variable_index + 1L
    variable_rows[[variable_index]] <- trajectory_by_variable

    reference_beta <- beta_aligned[start_id == reference_start]
    current_beta <- beta_aligned[start_id == sid]
    beta_compare <- merge(
      reference_beta[
        ,
        .(class, term, variable, reference_value = coefficient)
      ],
      current_beta[
        ,
        .(class, term, variable, current_value = coefficient)
      ],
      by = c("class", "term", "variable"),
      all = FALSE,
      sort = TRUE
    )

    reference_sigma <- sigma_aligned[start_id == reference_start]
    current_sigma <- sigma_aligned[start_id == sid]
    sigma_compare <- merge(
      reference_sigma[
        ,
        .(
          class, variable_row, variable_col,
          reference_value = covariance
        )
      ],
      current_sigma[
        ,
        .(
          class, variable_row, variable_col,
          current_value = covariance
        )
      ],
      by = c("class", "variable_row", "variable_col"),
      all = FALSE,
      sort = TRUE
    )

    trajectory_difference <-
      trajectory_compare$reference_z - trajectory_compare$current_z
    beta_difference <- beta_compare$reference_value - beta_compare$current_value
    sigma_difference <- sigma_compare$reference_value - sigma_compare$current_value

    consistency_rows[[sid]] <- data.table(
      version = VERSION_TAG,
      requested_g = g,
      reference_start_id = reference_start,
      start_id = sid,
      is_reference = as.integer(sid == reference_start),
      assignment_agreement = mean(aligned_assign == reference_assign),
      adjusted_rand_index = adjusted_rand_index(reference_assign, aligned_assign),
      posterior_MAE = mean(abs(reference_posterior - aligned_posterior)),
      posterior_RMSE = sqrt(mean((reference_posterior - aligned_posterior)^2)),
      trajectory_RMSE = sqrt(mean(trajectory_difference^2)),
      trajectory_MAE = mean(abs(trajectory_difference)),
      trajectory_max_abs = max(abs(trajectory_difference)),
      beta_RMSE = sqrt(mean(beta_difference^2)),
      beta_max_abs = max(abs(beta_difference)),
      Sigma_RMSE = sqrt(mean(sigma_difference^2)),
      Sigma_max_abs = max(abs(sigma_difference)),
      max_class_pct_point_difference = max(
        abs(class_compare$reference_pct - class_compare$current_pct),
        na.rm = TRUE
      ),
      max_prior_difference = max(
        abs(class_compare$reference_prior - class_compare$current_prior),
        na.rm = TRUE
      ),
      max_APPA_difference = max(
        abs(class_compare$reference_APPA - class_compare$current_APPA),
        na.rm = TRUE
      )
    )
  }

  consistency <- rbindlist(consistency_rows)
  variable_consistency <- rbindlist(variable_rows)

  best_logLik <- max(summaries$final_logLik, na.rm = TRUE)
  best_BIC <- min(summaries$BIC, na.rm = TRUE)
  multistart_summary <- merge(
    summaries,
    consistency,
    by = c("version", "requested_g", "start_id"),
    all.x = TRUE,
    sort = TRUE
  )
  multistart_summary[
    ,
    `:=`(
      logLik_rank = frank(-final_logLik, ties.method = "min"),
      delta_logLik_from_best = best_logLik - final_logLik,
      delta_BIC_from_best = BIC - best_BIC
    )
  ]
  setorder(multistart_summary, start_id)

  prefix <- file.path(output_dir, paste0(VERSION_TAG, "_"))
  safe_fwrite(
    multistart_summary,
    paste0(prefix, "multistart_summary_G", g, ".csv")
  )
  safe_fwrite(
    label_map,
    paste0(prefix, "label_map_G", g, ".csv")
  )
  safe_fwrite(
    consistency,
    paste0(prefix, "start_consistency_G", g, ".csv")
  )
  safe_fwrite(
    variable_consistency,
    paste0(prefix, "trajectory_consistency_G", g, ".csv")
  )
  safe_fwrite(
    trajectories_aligned,
    paste0(prefix, "trajectory_aligned_G", g, ".csv")
  )
  safe_fwrite(
    class_aligned,
    paste0(prefix, "class_aligned_G", g, ".csv")
  )
  safe_fwrite(
    beta_aligned,
    paste0(prefix, "beta_aligned_G", g, ".csv")
  )
  safe_fwrite(
    sigma_aligned,
    paste0(prefix, "sigma_aligned_G", g, ".csv")
  )

  writeLines(
    c(
      paste0("Version: ", VERSION_TAG),
      paste0("G: ", g),
      paste0("Reference start: ", reference_start),
      paste0("Reference rule: ", reference_rule),
      paste0("Completed: ", timestamp_now()),
      "Alignment cost: RMSE across all five fitted z-scale trajectories",
      "No outcome variable was read or used"
    ),
    file.path(output_dir, paste0(VERSION_TAG, "_alignment_report_G", g, ".txt"))
  )

  invisible(TRUE)
}

build_cross_g_comparison <- function(root = output_root) {
  paths <- list.files(
    root,
    pattern = paste0("^", VERSION_TAG, "_multistart_summary_G(3|4|5)\\.csv$"),
    recursive = TRUE,
    full.names = TRUE
  )
  if (length(paths) == 0L) return(invisible(FALSE))

  tables <- lapply(
    paths,
    function(path) {
      x <- fread(path, showProgress = FALSE)
      x[, source_file := relative_output_path(path, root)]
      x
    }
  )
  all_starts <- rbindlist(tables, fill = TRUE)
  setorder(all_starts, requested_g, start_id)
  safe_fwrite(
    all_starts,
    file.path(root, paste0(VERSION_TAG, "_multistart_summary_all.csv"))
  )

  comparison <- all_starts[
    ,
    {
      reference_rows <- which(is_reference == 1L)
      reference_start <- if (length(reference_rows) > 0L) {
        start_id[reference_rows[1L]]
      } else {
        NA_integer_
      }
      nonreference <- is_reference == 0L
      list(
        starts_completed = .N,
        converged_valid_full_g = sum(status == "CONVERGED_VALID_FULL_G"),
        reference_start_id = reference_start,
        stage2_logLik = stage2_logLik[1L],
        best_refined_logLik_minus_stage2 =
          max(refined_logLik_minus_stage2, na.rm = TRUE),
        best_logLik = max(final_logLik, na.rm = TRUE),
        median_logLik = median(final_logLik, na.rm = TRUE),
        worst_delta_logLik = max(delta_logLik_from_best, na.rm = TRUE),
        best_BIC = min(BIC, na.rm = TRUE),
        best_refined_BIC_minus_stage2 =
          min(refined_BIC_minus_stage2, na.rm = TRUE),
        median_BIC = median(BIC, na.rm = TRUE),
        min_class_pct_across_starts = min(min_class_pct_nonempty, na.rm = TRUE),
        min_APPA_across_starts = min(min_APPA_nonempty, na.rm = TRUE),
        min_assignment_agreement = if (any(nonreference)) {
          min(assignment_agreement[nonreference], na.rm = TRUE)
        } else {
          NA_real_
        },
        min_adjusted_rand_index = if (any(nonreference)) {
          min(adjusted_rand_index[nonreference], na.rm = TRUE)
        } else {
          NA_real_
        },
        max_trajectory_RMSE = if (any(nonreference)) {
          max(trajectory_RMSE[nonreference], na.rm = TRUE)
        } else {
          NA_real_
        },
        max_posterior_MAE = if (any(nonreference)) {
          max(posterior_MAE[nonreference], na.rm = TRUE)
        } else {
          NA_real_
        }
      )
    },
    by = requested_g
  ]
  setorder(comparison, requested_g)
  safe_fwrite(
    comparison,
    file.path(root, paste0(VERSION_TAG, "_G_comparison_all.csv"))
  )
  invisible(TRUE)
}


# ============================================================
# 13. One complete G x start_id task
# ============================================================

run_one_task <- function(task_row) {
  g <- as.integer(task_row$requested_g)
  start_id <- as.integer(task_row$start_id)
  paths <- task_paths(g, start_id, role_dir)

  if (task_is_complete(g, start_id, role_dir)) {
    x <- fread(paths$summary, showProgress = FALSE)
    x[, resumed_skip := 1L]
    return(x)
  }

  if (file.exists(pause_request_path)) {
    return(
      data.table(
        version = VERSION_TAG,
        machine_role = MACHINE_ROLE,
        requested_g = g,
        start_id = start_id,
        status = "NOT_STARTED_PAUSE_REQUESTED"
      )
    )
  }

  # Offset the second persistent lane so both fits do not allocate their
  # largest working objects in the same instant. Subsequent admission is
  # decided from live Windows memory plus observed prior-task peak memory.
  if (as.integer(LANE_ID) > 1L) {
    Sys.sleep((as.integer(LANE_ID) - 1L) * START_STAGGER_SECONDS)
  }

  wait_for_memory_gate(g, start_id, paths$stage_log)

  if (file.exists(paths$done)) file.remove(paths$done)
  if (file.exists(paths$error)) file.remove(paths$error)

  start_time <- Sys.time()
  seed_used <- seed_for_task(g, start_id)
  signature <- make_checkpoint_signature(g, start_id, seed_used)

  writeLines(
    c(
      paste0("Version: ", VERSION_TAG),
      paste0("Machine role: ", MACHINE_ROLE),
      paste0("Machine name: ", Sys.info()[["nodename"]]),
      paste0("CPU: ", cpu_name),
      paste0("Task: ", task_suffix(g, start_id)),
      paste0("Seed: ", seed_used),
      paste0("Controller PID: ", Sys.getpid()),
      paste0("Inner workers: ", length(INNER_PIDS)),
      paste0("Inner logical CPUs: ", paste(INNER_SLOTS, collapse = ",")),
      paste0("Start: ", timestamp_now()),
      "Status: RUNNING",
      "GPU used: FALSE",
      "Checkpoint: every accepted non-final EM iteration"
    ),
    paths$running
  )

  append_log(paths$stage_log, "task entered")

  fit_result <- tryCatch(
    run_parallel_one_start(
      g = g,
      start_id = start_id,
      seed_used = seed_used,
      paths = paths,
      signature = signature,
      inner_cluster = INNER_CL,
      class_chunks = split(
        seq_len(g),
        rep(seq_along(INNER_SLOTS), length.out = g)
      ),
      inner_pids = INNER_PIDS,
      static = V295_STATIC,
      pause_path = pause_request_path
    ),
    error = function(e) e
  )

  if (inherits(fit_result, "error")) {
    if (file.exists(paths$running)) file.remove(paths$running)
    error_text <- conditionMessage(fit_result)
    writeLines(
      c(
        paste0("Task: ", task_suffix(g, start_id)),
        "Status: ERROR",
        paste0("Error: ", error_text),
        paste0("Time: ", timestamp_now())
      ),
      paths$error
    )
    append_log(paths$stage_log, paste0("ERROR: ", error_text))
    return(
      data.table(
        version = VERSION_TAG,
        machine_role = MACHINE_ROLE,
        machine_name = Sys.info()[["nodename"]],
        requested_g = g,
        start_id = start_id,
        seed = seed_used,
        status = "ERROR",
        error_message = error_text
      )
    )
  }

  if (isTRUE(fit_result$paused)) {
    if (file.exists(paths$running)) file.remove(paths$running)
    append_log(paths$stage_log, "paused after safe checkpoint")
    return(
      data.table(
        version = VERSION_TAG,
        machine_role = MACHINE_ROLE,
        machine_name = Sys.info()[["nodename"]],
        requested_g = g,
        start_id = start_id,
        seed = seed_used,
        status = "PAUSED_AT_CHECKPOINT",
        em_iterations = fit_result$final_iteration,
        peak_working_set_mb = fit_result$peak_working_set_mb,
        peak_private_mb = fit_result$peak_private_mb
      )
    )
  }

  fit <- fit_result$fit

  safe_save_rds(fit, paths$full_fit, compress = TRUE)

  compact_fit <- fit
  compact_fit$reg <- NULL
  compact_fit$data.orig <- NULL
  compact_fit$data.norm <- NULL
  compact_fit$data.imputed <- NULL
  compact_fit$v284_metadata <- list(
    version = VERSION_TAG,
    machine_role = MACHINE_ROLE,
    machine_name = Sys.info()[["nodename"]],
    requested_g = g,
    start_id = start_id,
    seed = seed_used,
    input_md5 = input_md5,
    stop_reason = fit_result$stop_reason
  )
  safe_save_rds(compact_fit, paths$compact_fit, compress = TRUE)

  class_dt <- extract_class_quality(fit, g, start_id)
  patient_dt <- extract_patient_class(fit, V295_STATIC, g, start_id)
  trajectory_dt <- extract_trajectory(fit, g, start_id)
  beta_dt <- extract_beta(fit, g, start_id)
  sigma_dt <- extract_sigma(fit, g, start_id)

  safe_fwrite(class_dt, paths$class_quality)
  safe_fwrite(patient_dt, paths$patient_class)
  safe_fwrite(trajectory_dt, paths$trajectory)
  safe_fwrite(beta_dt, paths$beta)
  safe_fwrite(sigma_dt, paths$sigma)

  stored_g <- as.integer(fit$call$ng)
  posterior_g <- ncol(as.matrix(fit$posterior))
  occupied_g <- sum(class_dt$n > 0L)
  empty_classes <- class_dt[n == 0L, class]
  structure_exact <- stored_g == g && posterior_g == g
  all_occupied <- occupied_g == g
  full_g <- structure_exact && all_occupied
  package_converged <- as.numeric(fit$em["converged"]) == 1

  status <- if (package_converged && full_g) {
    "CONVERGED_VALID_FULL_G"
  } else if (package_converged && structure_exact && !all_occupied) {
    "CONVERGED_EMPTY_CLASS"
  } else if (package_converged && !structure_exact) {
    "CONVERGED_STRUCTURE_MISMATCH"
  } else if (!package_converged && full_g) {
    "NONCONVERGED_FULL_G"
  } else if (!package_converged && structure_exact && !all_occupied) {
    "NONCONVERGED_EMPTY_CLASS"
  } else {
    "NONCONVERGED_STRUCTURE_MISMATCH"
  }

  end_time <- Sys.time()
  max_post <- apply(as.matrix(fit$posterior), 1L, max, na.rm = TRUE)
  nonempty_class <- class_dt[n > 0L]

  summary_dt <- data.table(
    version = VERSION_TAG,
    engine = ENGINE_TAG,
    machine_role = MACHINE_ROLE,
    machine_name = Sys.info()[["nodename"]],
    cpu_name = cpu_name,
    requested_g = g,
    start_id = start_id,
    seed = seed_used,
    status = status,
    package_converged = as.integer(package_converged),
    stored_g = stored_g,
    posterior_g = posterior_g,
    occupied_g = occupied_g,
    empty_classes = if (length(empty_classes) == 0L) "" else paste(empty_classes, collapse = ","),
    structural_full_g = as.integer(full_g),
    d = MODEL_DEGREE,
    scaling = 0L,
    pruning = 1L,
    delete_empty = 0L,
    nstart_within_task = 1L,
    tol = REFINE_TOL,
    technical_maxit = NONBINDING_MAXIT,
    em_iterations = fit_result$final_iteration,
    stop_reason = fit_result$stop_reason,
    em_logLik_before_pruning = fit_result$em_logLik_before_pruning,
    final_logLik = as.numeric(fit$logLik),
    AIC = safe_ic_value(fit, "AIC"),
    BIC = safe_ic_value(fit, "BIC"),
    CAIC = safe_ic_value(fit, "CAIC"),
    ssBIC = safe_ic_value(fit, "ssBIC"),
    HQIC = safe_ic_value(fit, "HQIC"),
    npar = as.integer(fit$npar),
    min_class_pct_nonempty = min(nonempty_class$pct),
    min_APPA_nonempty = min(nonempty_class$APPA),
    mean_max_posterior = mean(max_post),
    median_max_posterior = median(max_post),
    pct_max_posterior_ge_070 = 100 * mean(max_post >= 0.70),
    pct_max_posterior_ge_090 = 100 * mean(max_post >= 0.90),
    preferred_class_size_ok = as.integer(min(nonempty_class$pct) >= PREFERRED_MIN_CLASS_PCT),
    preferred_APPA_ok = as.integer(min(nonempty_class$APPA) >= PREFERRED_MIN_APPA),
    resumed = as.integer(fit_result$resumed),
    resumed_from_iteration = as.integer(fit_result$resumed_from_iteration),
    checkpoint_every_accepted_iteration = CHECKPOINT_EVERY_ACCEPTED_ITER,
    controller_pid = Sys.getpid(),
    inner_worker_count = length(INNER_PIDS),
    inner_logical_cpus = paste(INNER_SLOTS, collapse = ","),
    peak_working_set_mb = fit_result$peak_working_set_mb,
    peak_private_mb = fit_result$peak_private_mb,
    gpu_used = 0L,
    input_md5 = input_md5,
    n_patients = EXPECTED_PATIENTS,
    n_rows = EXPECTED_LONG_ROWS,
    start_time = format(start_time, "%Y-%m-%d %H:%M:%S"),
    end_time = format(end_time, "%Y-%m-%d %H:%M:%S"),
    elapsed_seconds = as.numeric(difftime(end_time, start_time, units = "secs")),
    full_model_file = relative_output_path(paths$full_fit),
    compact_model_file = relative_output_path(paths$compact_fit),
    resumed_skip = 0L,
    error_message = NA_character_
  )

  safe_fwrite(summary_dt, paths$summary)

  resume_dt <- summary_dt[
    ,
    .(
      version,
      machine_role,
      machine_name,
      requested_g,
      start_id,
      seed,
      resumed,
      resumed_from_iteration,
      em_iterations,
      checkpoint_every_accepted_iteration,
      input_md5,
      completed_at = timestamp_now()
    )
  ]
  safe_fwrite(resume_dt, paths$resume_audit)

  writeLines(
    c(
      paste0("Version: ", VERSION_TAG),
      paste0("Task: ", task_suffix(g, start_id)),
      paste0("Status: ", status),
      paste0("Seed: ", seed_used),
      paste0("Iterations: ", fit_result$final_iteration),
      paste0("Final logLik: ", format(as.numeric(fit$logLik), digits = 16)),
      paste0("Stop reason: ", fit_result$stop_reason),
      paste0("End: ", timestamp_now())
    ),
    paths$done
  )

  for (cp in c(paths$checkpoint_a, paths$checkpoint_b)) {
    if (file.exists(cp)) file.remove(cp)
  }
  if (file.exists(paths$running)) file.remove(paths$running)
  append_log(paths$stage_log, paste0("task finished: ", status))
  summary_dt
}


# ============================================================
# 14. Initialize each persistent outer lane and inner CPU pool
# ============================================================

initialize_outer_lane <- function(
  lane_id,
  inner_slots,
  controller_slot,
  input_path,
  helper_path,
  powershell_path
) {
  Sys.setenv(
    OMP_NUM_THREADS = "1",
    OPENBLAS_NUM_THREADS = "1",
    MKL_NUM_THREADS = "1",
    BLAS_NUM_THREADS = "1",
    VECLIB_MAXIMUM_THREADS = "1",
    NUMEXPR_NUM_THREADS = "1"
  )
  suppressPackageStartupMessages({
    library(data.table)
    library(gbmt)
    library(parallel)
  })
  data.table::setDTthreads(1L)

  assign("LANE_ID", as.integer(lane_id), envir = .GlobalEnv)
  assign("INNER_SLOTS", as.integer(inner_slots), envir = .GlobalEnv)

  apply_affinity(Sys.getpid(), as.integer(controller_slot))

  static <- prepare_static_data(input_path)
  assign("V295_STATIC", static, envir = .GlobalEnv)

  inner_console <- file.path(
    role_dir,
    paste0(VERSION_TAG, "_inner_console_lane", lane_id, ".log")
  )
  inner_cl <- parallel::makePSOCKcluster(length(inner_slots), outfile = inner_console)
  inner_pids <- unlist(
    parallel::clusterEvalQ(inner_cl, Sys.getpid()),
    use.names = FALSE
  )

  if (length(inner_pids) != length(inner_slots)) {
    parallel::stopCluster(inner_cl)
    stop("Inner PID count differs from allocated CPU-slot count.")
  }

  for (i in seq_along(inner_pids)) {
    apply_affinity(inner_pids[i], inner_slots[i])
  }

  parallel::clusterEvalQ(
    inner_cl,
    {
      Sys.setenv(
        OMP_NUM_THREADS = "1",
        OPENBLAS_NUM_THREADS = "1",
        MKL_NUM_THREADS = "1",
        BLAS_NUM_THREADS = "1",
        VECLIB_MAXIMUM_THREADS = "1",
        NUMEXPR_NUM_THREADS = "1"
      )
      suppressPackageStartupMessages({
        library(gbmt)
        library(data.table)
      })
      data.table::setDTthreads(1L)
      NULL
    }
  )

  assign("V295_STATIC", static, envir = .GlobalEnv)
  parallel::clusterExport(
    inner_cl,
    varlist = c("V295_STATIC", "inner_init_chunk", "inner_mstep_chunk"),
    envir = .GlobalEnv
  )

  assign("INNER_CL", inner_cl, envir = .GlobalEnv)
  assign("INNER_PIDS", as.integer(inner_pids), envir = .GlobalEnv)

  list(
    lane_id = lane_id,
    controller_pid = Sys.getpid(),
    controller_slot = controller_slot,
    inner_pids = inner_pids,
    inner_slots = inner_slots
  )
}

shutdown_outer_lane <- function() {
  if (exists("INNER_CL", envir = .GlobalEnv, inherits = FALSE)) {
    try(parallel::stopCluster(get("INNER_CL", envir = .GlobalEnv)), silent = TRUE)
    rm("INNER_CL", envir = .GlobalEnv)
  }
  TRUE
}

run_task_pool <- function(task_subset, lane_slots, controller_slot_values, phase_tag) {
  task_subset <- as.data.table(task_subset)
  if (nrow(task_subset) == 0L) return(list())

  lane_count <- length(lane_slots)
  if (
    lane_count < 1L ||
    lane_count != length(controller_slot_values) ||
    lane_count > MAX_SUPPORTED_SIMULTANEOUS_FITS
  ) {
    stop("Invalid outer-lane configuration.")
  }
  if (any(lengths(lane_slots) < 1L)) stop("An outer lane has no inner CPU slot.")
  all_assigned_slots <- c(unlist(lane_slots), controller_slot_values)
  if (anyDuplicated(all_assigned_slots)) stop("A logical CPU was assigned more than once.")

  outer_console <- file.path(
    role_dir,
    paste0(VERSION_TAG, "_outer_console_", phase_tag, ".log")
  )
  outer_cl <- parallel::makePSOCKcluster(lane_count, outfile = outer_console)
  stopped <- FALSE

  on.exit(
    {
      if (!isTRUE(stopped)) {
        try(parallel::clusterEvalQ(outer_cl, shutdown_outer_lane()), silent = TRUE)
        try(parallel::stopCluster(outer_cl), silent = TRUE)
      }
    },
    add = TRUE
  )

  outer_exports <- c(
    "VERSION_TAG", "ENGINE_TAG", "MODEL_DEGREE", "REFINE_TOL",
    "NONBINDING_MAXIT", "Z_MARKERS", "REQUIRED_COLS",
    "EXPECTED_PATIENTS", "EXPECTED_LONG_ROWS", "STAGE2_MASTER_SEED",
    "START_SEED_STEP", "CHECKPOINT_STATE_VERSION",
    "CHECKPOINT_EVERY_ACCEPTED_ITER", "START_STAGGER_SECONDS",
    "MEMORY_MIN_AVAILABLE_GB",
    "MEMORY_MIN_AVAILABLE_FRACTION", "MEMORY_MAX_COMMIT_FRACTION",
    "MEMORY_PEAK_SAFETY_FACTOR", "MEMORY_POLL_SECONDS",
    "PROCESS_PRIORITY_CLASS", "PREFERRED_MIN_CLASS_PCT",
    "PREFERRED_MIN_APPA", "SCALING_TABLE", "MACHINE_ROLE", "cpu_name",
    "input_rds", "input_md5", "output_root", "role_dir",
    "pause_request_path", "windows_helper", "powershell_exe",
    "timestamp_now", "append_log", "safe_save_rds", "safe_fwrite",
    "relative_output_path", "seed_for_task", "task_suffix", "task_paths",
    "task_is_complete", "make_checkpoint_signature",
    "make_checkpoint_state", "validate_checkpoint_state",
    "read_checkpoint_slot", "load_latest_checkpoint", "save_checkpoint_dual",
    "prepare_static_data", "inner_init_chunk", "inner_mstep_chunk",
    "flatten_class_results", "read_ps_csv", "run_ps_helper", "query_memory",
    "query_process_memory", "apply_affinity", "read_observed_peak_mb",
    "wait_for_memory_gate", "measure_task_resources",
    "run_parallel_one_start", "safe_ic_value", "extract_class_quality",
    "extract_patient_class", "extract_trajectory", "extract_beta",
    "extract_sigma", "run_one_task", "initialize_outer_lane",
    "shutdown_outer_lane"
  )

  parallel::clusterExport(outer_cl, outer_exports, envir = .GlobalEnv)

  lane_meta <- parallel::clusterMap(
    outer_cl,
    fun = initialize_outer_lane,
    lane_id = seq_len(lane_count),
    inner_slots = lane_slots,
    controller_slot = as.integer(controller_slot_values),
    MoreArgs = list(
      input_path = input_rds,
      helper_path = windows_helper,
      powershell_path = powershell_exe
    ),
    SIMPLIFY = FALSE,
    USE.NAMES = FALSE,
    RECYCLE = FALSE
  )

  lane_dt <- rbindlist(
    lapply(
      lane_meta,
      function(x) {
        data.table(
          phase = phase_tag,
          lane_id = x$lane_id,
          controller_pid = x$controller_pid,
          controller_slot = x$controller_slot,
          inner_worker_count = length(x$inner_pids),
          inner_pids = paste(x$inner_pids, collapse = ","),
          inner_slots = paste(x$inner_slots, collapse = ",")
        )
      }
    )
  )
  safe_fwrite(
    lane_dt,
    file.path(role_dir, paste0(VERSION_TAG, "_cpu_lanes_", phase_tag, ".csv"))
  )

  task_rows <- lapply(seq_len(nrow(task_subset)), function(i) task_subset[i])
  results <- parallel::parLapplyLB(outer_cl, task_rows, run_one_task)

  parallel::clusterEvalQ(outer_cl, shutdown_outer_lane())
  parallel::stopCluster(outer_cl)
  stopped <- TRUE
  results
}

choose_post_calibration_lanes <- function(observed_peak_private_mb) {
  if (!is.finite(observed_peak_private_mb) || observed_peak_private_mb <= 0) {
    stop("A finite positive calibration peak_private_mb is required.")
  }

  gc(full = TRUE)
  Sys.sleep(MEMORY_POST_CALIBRATION_SETTLE_SECONDS)
  mem <- query_memory()
  system_reserve_mb <- max(
    MEMORY_MIN_AVAILABLE_GB * 1024,
    MEMORY_MIN_AVAILABLE_FRACTION * mem$TotalPhysicalMB
  )
  per_fit_reserve_mb <- observed_peak_private_mb * MEMORY_PEAK_SAFETY_FACTOR
  physical_capacity <- floor(
    max(0, mem$TotalPhysicalMB - system_reserve_mb) / per_fit_reserve_mb
  )
  available_capacity <- floor(
    max(0, mem$AvailablePhysicalMB - system_reserve_mb) / per_fit_reserve_mb
  )
  commit_headroom_mb <- max(
    0,
    MEMORY_MAX_COMMIT_FRACTION * mem$CommitLimitMB - mem$CommitUsedMB
  )
  commit_capacity <- floor(commit_headroom_mb / per_fit_reserve_mb)
  safe_capacity <- min(
    MAX_SUPPORTED_SIMULTANEOUS_FITS,
    physical_capacity,
    available_capacity,
    commit_capacity
  )
  chosen <- 2L

  lane_rows <- split(
    seq_len(nrow(heavy_slots)),
    rep(seq_len(2L), length.out = nrow(heavy_slots))
  )
  chosen_lane_slots <- lapply(
    lane_rows,
    function(idx) as.integer(heavy_slots$LogicalProcessorIndex[idx])
  )
  chosen_controller_slots <- as.integer(controller_slots$LogicalProcessorIndex)

  audit <- data.table(
    measured_at = timestamp_now(),
    observed_peak_private_mb = observed_peak_private_mb,
    safety_factor = MEMORY_PEAK_SAFETY_FACTOR,
    per_fit_reserve_mb = per_fit_reserve_mb,
    total_physical_mb = mem$TotalPhysicalMB,
    available_physical_mb = mem$AvailablePhysicalMB,
    system_reserve_mb = system_reserve_mb,
    commit_used_mb = mem$CommitUsedMB,
    commit_limit_mb = mem$CommitLimitMB,
    max_commit_fraction = MEMORY_MAX_COMMIT_FRACTION,
    physical_capacity = physical_capacity,
    available_capacity = available_capacity,
    commit_capacity = commit_capacity,
    safe_capacity_by_original_rule = safe_capacity,
    forced_two_lanes = as.integer(FORCE_TWO_LANES),
    chosen_simultaneous_fits = chosen,
    lane_1_inner_slots = paste(chosen_lane_slots[[1L]], collapse = ","),
    lane_2_inner_slots = if (chosen == 2L) {
      paste(chosen_lane_slots[[2L]], collapse = ",")
    } else {
      ""
    },
    controller_slots = paste(chosen_controller_slots, collapse = ",")
  )

  list(
    lane_count = chosen,
    lane_slots = chosen_lane_slots,
    controller_slots = chosen_controller_slots,
    audit = audit
  )
}


# ============================================================
# 15. Save preflight configuration
# ============================================================

task_table[, seed := mapply(seed_for_task, requested_g, start_id)]
task_table[, task := mapply(task_suffix, requested_g, start_id)]
task_table[, machine_role := MACHINE_ROLE]
task_table[, priority_order := .I]

safe_fwrite(task_table, file.path(role_dir, "V295_task_manifest.csv"))
safe_fwrite(
  as.data.table(STAGE2_BENCHMARK),
  file.path(output_root, "V295_stage2_benchmark.csv")
)

preflight_config <- data.table(
  parameter = c(
    "VERSION_TAG", "ENGINE_TAG", "EXPECTED_MACHINE_ROLE",
    "MACHINE_ROLE", "cpu_name",
    "MODEL_DEGREE", "REFINE_TOL", "NONBINDING_MAXIT",
    "SCALING", "PRUNING", "DELETE_EMPTY", "NSTART_PER_TASK",
    "MAX_SUPPORTED_SIMULTANEOUS_FITS", "FORCE_TWO_LANES", "CALIBRATION_G",
    "CALIBRATION_START_ID", "STEP4_G_VALUES", "STEP4_PRIORITY_G",
    "MAX_CLASS_WORKERS_PER_FIT", "HEAVY_CPU_SLOTS",
    "CHECKPOINT_EVERY_ACCEPTED_ITER", "INPUT_FILE", "INPUT_MD5",
    "EXPECTED_INPUT_MD5", "OUTPUT_ROOT",
    "GPU_USED", "R_VERSION", "GBMT_VERSION", "DATATABLE_VERSION"
  ),
  value = c(
    VERSION_TAG, ENGINE_TAG, EXPECTED_MACHINE_ROLE,
    MACHINE_ROLE, cpu_name,
    MODEL_DEGREE, REFINE_TOL, NONBINDING_MAXIT,
    0, TRUE, FALSE, 1,
    MAX_SUPPORTED_SIMULTANEOUS_FITS, FORCE_TWO_LANES, CALIBRATION_G,
    CALIBRATION_START_ID,
    paste(STEP4_G_VALUES, collapse = ","),
    paste(STEP4_PRIORITY_G, collapse = ","),
    MAX_CLASS_WORKERS_PER_FIT,
    nrow(heavy_slots),
    CHECKPOINT_EVERY_ACCEPTED_ITER,
    basename(input_rds), input_md5,
    EXPECTED_INPUT_MD5,
    normalizePath(output_root, winslash = "/", mustWork = TRUE),
    FALSE, R.version.string, as.character(packageVersion("gbmt")),
    as.character(packageVersion("data.table"))
  )
)

safe_fwrite(preflight_config, file.path(role_dir, "V295_preflight_config.csv"))


# ============================================================
# 16. Calibrate memory, then run all pending starts
# ============================================================

initial_pending <- task_table[
  !vapply(
    seq_len(.N),
    function(i) task_is_complete(requested_g[i], start_id[i], role_dir),
    logical(1)
  )
]

cat("============================================================\n")
cat("V295 STEP 4A MULTIPLE-START REFINEMENT\n")
cat("============================================================\n")
cat("Machine role: ", MACHINE_ROLE, "\n", sep = "")
cat("CPU: ", cpu_name, "\n", sep = "")
cat("Selected G values: ", paste(STEP4_G_VALUES, collapse = ", "), "\n", sep = "")
cat("Tasks assigned: ", paste(task_table$task, collapse = ", "), "\n", sep = "")
cat(
  "Pending: ",
  if (nrow(initial_pending) == 0L) "none" else paste(initial_pending$task, collapse = ", "),
  "\n",
  sep = ""
)
cat("GPU used by gbmt: no\n")
cat("Output: ", normalizePath(role_dir, winslash = "/", mustWork = TRUE), "\n\n", sep = "")

master_start <- Sys.time()
results <- list()
calibration_complete <- task_is_complete(
  CALIBRATION_G,
  CALIBRATION_START_ID,
  role_dir
)

if (!calibration_complete) {
  stop(
    paste0(
      "The completed V295 G6-S01 memory-reference task was not found.\n",
      "This G3-G5 continuation script will not schedule or refit G6.\n",
      "Required existing files are checked through task_is_complete(6, 1):\n",
      normalizePath(role_dir, winslash = "/", mustWork = TRUE)
    )
  )
}

allocation <- NULL
if (calibration_complete) {
  observed_peak <- read_observed_peak_mb(role_dir)
  allocation <- choose_post_calibration_lanes(observed_peak)
  safe_fwrite(
    allocation$audit,
    file.path(role_dir, "V295_memory_allocation.csv")
  )

  remaining_pending <- task_table[
    !vapply(
      seq_len(.N),
      function(i) task_is_complete(requested_g[i], start_id[i], role_dir),
      logical(1)
    )
  ]

  if (nrow(remaining_pending) > 0L && !file.exists(pause_request_path)) {
    cat(
      "Running pending G3-G5 starts with ", allocation$lane_count,
      " simultaneous fits.\n",
      sep = ""
    )
    refinement_results <- run_task_pool(
      task_subset = remaining_pending,
      lane_slots = allocation$lane_slots,
      controller_slot_values = allocation$controller_slots,
      phase_tag = paste0("refinement_", allocation$lane_count, "lane")
    )
    results <- c(results, refinement_results)
  }
} else {
  stop("Unexpected unreachable state: G6-S01 memory reference is incomplete.")
}

execution_config <- copy(preflight_config)
execution_config <- rbind(
  execution_config,
  data.table(
    parameter = c(
      "CALIBRATION_COMPLETE",
      "CHOSEN_SIMULTANEOUS_FITS",
      "CHOSEN_LANE_1_SLOTS",
      "CHOSEN_LANE_2_SLOTS",
      "CHOSEN_CONTROLLER_SLOTS"
    ),
    value = c(
      calibration_complete,
      if (is.null(allocation)) NA_character_ else as.character(allocation$lane_count),
      if (is.null(allocation)) NA_character_ else paste(allocation$lane_slots[[1L]], collapse = ","),
      if (is.null(allocation) || allocation$lane_count < 2L) "" else paste(allocation$lane_slots[[2L]], collapse = ","),
      if (is.null(allocation)) NA_character_ else paste(allocation$controller_slots, collapse = ",")
    )
  )
)
safe_fwrite(execution_config, file.path(role_dir, "V295_execution_config.csv"))


# ============================================================
# 17. Analyze completed starts and aggregate both machines
# ============================================================

analysis_dirs <- file.path(output_root, c("laptop", "desktop"))
for (analysis_dir in analysis_dirs[dir.exists(analysis_dirs)]) {
  for (g_value in 3:5) {
    analyze_completed_g(as.integer(g_value), analysis_dir)
  }
}

summary_files <- list.files(
  output_root,
  pattern = paste0("^", VERSION_TAG, "_task_summary_G(3|4|5)_S(01|02|03|04|05)\\.csv$"),
  recursive = TRUE,
  full.names = TRUE
)

summary_list <- lapply(
  summary_files,
  function(path) tryCatch(fread(path, showProgress = FALSE), error = function(e) NULL)
)
summary_list <- Filter(Negate(is.null), summary_list)

if (length(summary_list) > 0L) {
  combined <- rbindlist(summary_list, fill = TRUE)
  setorder(combined, requested_g, start_id)
  safe_fwrite(combined, file.path(output_root, "V295_task_summary_all.csv"))
} else {
  combined <- data.table()
}

build_cross_g_comparison(output_root)

master_end <- Sys.time()

writeLines(
  c(
    paste0("Version: ", VERSION_TAG),
    paste0("Machine role: ", MACHINE_ROLE),
    paste0("Machine: ", Sys.info()[["nodename"]]),
    paste0("CPU: ", cpu_name),
    paste0("Start: ", format(master_start, "%Y-%m-%d %H:%M:%S")),
    paste0("End: ", format(master_end, "%Y-%m-%d %H:%M:%S")),
    paste0(
      "Elapsed hours: ",
      round(as.numeric(difftime(master_end, master_start, units = "hours")), 3)
    ),
    paste0("Selected G values: ", paste(STEP4_G_VALUES, collapse = ",")),
    paste0("Assigned tasks: ", paste(task_table$task, collapse = ",")),
    paste0(
      "Completed assigned tasks: ",
      sum(
        vapply(
          seq_len(nrow(task_table)),
          function(i) task_is_complete(task_table$requested_g[i], task_table$start_id[i], role_dir),
          logical(1)
        )
      ),
      "/", nrow(task_table)
    ),
    "GPU used: FALSE (gbmt 0.1.4 has no GPU implementation)",
    paste0("Output: ", normalizePath(role_dir, winslash = "/", mustWork = TRUE))
  ),
  file.path(role_dir, "V295_machine_report.txt")
)

cat("\n============================================================\n")
cat("V295 SESSION FINISHED\n")
cat("============================================================\n")
cat("Machine role: ", MACHINE_ROLE, "\n", sep = "")
cat("Output: ", normalizePath(role_dir, winslash = "/", mustWork = TRUE), "\n", sep = "")
cat("To request a safe pause during a future session, create:\n")
cat(normalizePath(pause_request_path, winslash = "/", mustWork = FALSE), "\n")
cat("Delete that pause file before restarting.\n")
cat("The combined summary is V295_multistart_summary_all.csv.\n")
cat("============================================================\n")
