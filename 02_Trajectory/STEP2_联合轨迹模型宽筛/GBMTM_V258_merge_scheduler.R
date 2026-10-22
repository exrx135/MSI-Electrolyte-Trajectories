# ============================================================
# V258 - GBMTM G5-G20 DYNAMIC-SCHEDULER CHECKPOINT/RESUME SCREEN
# ============================================================
#
# Study model:
#   MSI + Na + K + Cl + HCO3 jointly enter ONE GBMTM model.
#
# Purpose:
#   Continue the frozen V251 broad screen with power-safe EM checkpoints
#   plus hardware-aware dynamic CPU scheduling and memory-aware task launch.
#   Final output filenames and schemas remain V251-compatible so old and
#   newly completed G models can be merged directly.
#
# Scientific settings inherited unchanged from V251:
#   - Requested G values: 5:20 (16 models)
#   - worker count is derived from pending tasks and detected logical CPUs
#   - 5-second global stagger between actual task starts
#   - physical cores are used before SMT sibling threads
#   - on hybrid CPUs: P-core > E-core > LP-E core > P-core SMT
#   - safe-boundary core promotion after completed+checkpointed EM iterations
#   - d = 2
#   - scaling = 0
#   - pruning = TRUE
#   - nstart = 1 (one randomized initialization)
#   - tol = 10 (deliberately relaxed FIRST-STAGE screening tolerance)
#   - delete.empty = FALSE
#   - NO automatic random-seed retry
#   - NO wall-clock timeout
#   - NO small/G-specific maxit budget
#   - maxit uses .Machine$integer.max only because gbmt() requires
#     a finite integer argument. It is intended to be non-binding.
#
# IMPORTANT:
#   1) A requested GX will NOT be silently converted into X-1 by this
#      script. delete.empty=FALSE prevents the gbmt empty-class deletion
#      step from intentionally removing an empty component.
#   2) If GX finishes with X stored components but one or more modal
#      groups have zero patients, the fit is SAVED and labeled
#      CONVERGED_EMPTY_CLASS. It is not accepted as a valid GX for IC
#      ranking, but it is also NOT relabeled as X-1 and NOT rerun.
#   3) If the stored model structure itself does not equal requested X
#      despite delete.empty=FALSE, the fit is saved and labeled
#      CONVERGED_STRUCTURE_MISMATCH for diagnosis.
#   4) tol=10 is for broad screening only. Final shortlisted models must
#      later be refitted with strict tolerance (e.g. 1e-4) and multiple
#      independent starts.
#
# Windows path rule:
#   Never hard-code a Windows username.
#   Sys.getenv("USERPROFILE") is used as the current Windows user profile root.
#   Required relative hierarchy is preserved exactly:
#   Desktop/Research/MSI/提取轨迹文件/GBMTM_stage1_output
#   Existing and new Stage-2 results are merged in exactly:
#   Desktop/Research/MSI/提取轨迹文件/GBMTM_stage1_output/GBMTM_V251_G5_G20_screen
#
# Output naming rule:
#   Version_Content_GX
#   e.g.
#     V251_task_summary_G5.csv
#     V251_class_quality_G5.csv
#     V251_EM_runs_G5.csv
#     V251_EM_progress_G5.log
#
# ============================================================


options(stringsAsFactors = FALSE)


# ============================================================
# 0. Fixed configuration
# ============================================================

VERSION_TAG <- "V251"
ENGINE_TAG <- "V258"

GROUPS_TO_SCREEN <- 5:20

# No arbitrary fixed worker count is imposed. The active PSOCK pool is
# min(number of pending G tasks, number of detected logical processors).
# Physical cores are always assigned before SMT sibling threads.
MODEL_DEGREE <- 2L

SCREEN_TOL <- 10

NSTART <- 1L

# Frozen project rule: actual task launches are separated globally by 5 s.
START_STAGGER_SECONDS <- 5L

# Memory-aware launch guard. This does not kill or suspend running fits;
# it only delays the next not-yet-loaded model when Windows RAM/commit
# headroom is critically low. The thresholds are hardware-relative.
MEMORY_MIN_AVAILABLE_GB <- 4
MEMORY_MIN_AVAILABLE_FRACTION <- 0.125
MEMORY_MAX_COMMIT_FRACTION <- 0.88
MEMORY_POLL_SECONDS <- 10L

# Worker priority is elevated but never set to Realtime.
PROCESS_PRIORITY_CLASS <- "AboveNormal"

# gbmt() requires finite integer maxit.
# This is a technical sentinel rather than a practical stopping budget.
NONBINDING_MAXIT <- .Machine$integer.max

MASTER_SEED <- 20260904L

EXPECTED_PATIENTS <- 9804L
EXPECTED_LONG_ROWS <- 117648L

PREFERRED_MIN_CLASS_PCT <- 5
PREFERRED_MIN_APPA <- 0.70

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


# ============================================================
# 1. Required packages
# ============================================================

required_packages <- c(
  "data.table",
  "gbmt"
)

missing_packages <- required_packages[
  !vapply(
    required_packages,
    requireNamespace,
    logical(1),
    quietly = TRUE
  )
]

if (length(missing_packages) > 0L) {
  stop(
    paste0(
      "Missing required R package(s): ",
      paste(missing_packages, collapse = ", "),
      "\nInstall them before running V251."
    )
  )
}

suppressPackageStartupMessages({
  library(data.table)
  library(gbmt)
  library(parallel)
})


# ============================================================
# 2. Prevent nested threading inside each PSOCK worker
# ============================================================

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1",
  BLAS_NUM_THREADS = "1",
  VECLIB_MAXIMUM_THREADS = "1",
  NUMEXPR_NUM_THREADS = "1"
)

data.table::setDTthreads(1L)


# ============================================================
# 3. Verify gbmt() interface before running
# ============================================================

gbmt_formals <- names(
  formals(
    gbmt::gbmt
  )
)

required_gbmt_args <- c(
  "nstart",
  "tol",
  "maxit",
  "delete.empty"
)

missing_gbmt_args <- setdiff(
  required_gbmt_args,
  gbmt_formals
)

if (length(missing_gbmt_args) > 0L) {
  stop(
    paste0(
      "Installed gbmt::gbmt() does not expose required argument(s): ",
      paste(missing_gbmt_args, collapse = ", "),
      "\nInstalled gbmt version: ",
      as.character(packageVersion("gbmt")),
      "\nStop here. Do not run with a different interface."
    )
  )
}


if (!identical(
  as.character(packageVersion("gbmt")),
  "0.1.4"
)) {
  stop(
    paste0(
      "V252 checkpoint engine was source-matched to gbmt 0.1.4.\n",
      "Installed gbmt version: ",
      as.character(packageVersion("gbmt")),
      "\nStop here rather than silently changing the model engine."
    )
  )
}

required_gbmt_internals <- c(
  "scaleFun",
  "scaleFun_inv",
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
      !inherits(
        try(
          getFromNamespace(nm, "gbmt"),
          silent = TRUE
        ),
        "try-error"
      )
    },
    logical(1)
  )
]

if (length(missing_gbmt_internals) > 0L) {
  stop(
    paste0(
      "Required gbmt 0.1.4 internal function(s) are unavailable: ",
      paste(missing_gbmt_internals, collapse = ", "),
      "\nStop here; checkpoint execution will not guess replacement logic."
    )
  )
}


# ============================================================
# 4. Resolve exact Windows profile root using USERPROFILE
# ============================================================

user_home <- Sys.getenv(
  "USERPROFILE",
  unset = ""
)

if (!nzchar(user_home)) {
  stop(
    paste0(
      "Windows USERPROFILE environment variable is empty.\n",
      "The script will not guess OneDrive, another drive, ",
      "another username, Documents, or another folder hierarchy."
    )
  )
}

if (!dir.exists(user_home)) {
  stop(
    paste0(
      "Resolved USERPROFILE directory does not exist:\n",
      normalizePath(
        user_home,
        winslash = "/",
        mustWork = FALSE
      ),
      "\nThe script will not guess OneDrive, another drive, ",
      "another username, Documents, or another folder hierarchy."
    )
  )
}

base_dir <- file.path(
  user_home,
  "Desktop",
  "Research",
  "MSI",
  "提取轨迹文件",
  "GBMTM_stage1_output"
)

input_rds <- file.path(
  base_dir,
  "MIMIC_GBMTM_long_scaled.rds"
)

output_dir <- file.path(
  base_dir,
  "GBMTM_V251_G5_G20_screen"
)

if (!dir.exists(base_dir)) {
  stop(
    paste0(
      "Required base directory does not exist.\n",
      "Checked exactly:\n",
      normalizePath(
        base_dir,
        winslash = "/",
        mustWork = FALSE
      ),
      "\nThe script will not guess OneDrive, another drive, ",
      "another username, or another folder hierarchy."
    )
  )
}

if (!file.exists(input_rds)) {
  stop(
    paste0(
      "Required input file does not exist.\n",
      "Checked exactly:\n",
      normalizePath(
        input_rds,
        winslash = "/",
        mustWork = FALSE
      )
    )
  )
}

dir.create(
  output_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(output_dir)) {
  stop(
    paste0(
      "Could not create output directory:\n",
      normalizePath(
        output_dir,
        winslash = "/",
        mustWork = FALSE
      )
    )
  )
}


# ============================================================
# 4A. Windows hardware-aware dynamic scheduler
# ============================================================
#
# Frozen scheduling policy implemented here:
#   - no arbitrary fixed task count;
#   - actual model starts are globally staggered by 5 seconds;
#   - when tasks <= physical cores, one task per physical core;
#   - SMT sibling threads are used only after all physical cores are occupied;
#   - on supported Intel hybrid CPUs: P > E > LP-E > P-SMT;
#   - when a better core becomes free, the lower-tier task with the lowest
#     completed EM iteration is reserved for promotion; ties are broken by
#     larger recent mean delta-LL, then longer mean iteration time;
#   - migration happens only after a completed accepted EM iteration has been
#     checkpointed successfully;
#   - CPU scheduling changes execution speed only. They do not alter model
#     parameters, initialization, stopping rules, seeds, or V251 outputs.
#
# Windows GetSystemCpuSetInformation is used to discover actual CPU topology.
# ProcessorAffinity is then applied to individual PSOCK worker PIDs.
# ============================================================

if (.Platform$OS.type != "windows") {
  stop(
    paste0(
      "V258 dynamic scheduler is designed for the project's Windows machines.\n",
      "Detected OS type: ",
      .Platform$OS.type,
      "\nStop here rather than silently disabling the frozen scheduling rules."
    )
  )
}

powershell_exe <- Sys.which("powershell")

if (!nzchar(powershell_exe)) {
  powershell_exe <- Sys.which("pwsh")
}

if (!nzchar(powershell_exe)) {
  stop(
    "PowerShell/pwsh was not found. V258 cannot query topology or set affinity."
  )
}

scheduler_dir <- file.path(
  output_dir,
  "V258_scheduler"
)

dir.create(
  scheduler_dir,
  recursive = TRUE,
  showWarnings = FALSE
)

if (!dir.exists(scheduler_dir)) {
  stop(
    paste0(
      "Could not create scheduler directory:\n",
      normalizePath(
        scheduler_dir,
        winslash = "/",
        mustWork = FALSE
      )
    )
  )
}

scheduler_helper_ps1 <- file.path(
  scheduler_dir,
  "V258_windows_cpu_helper.ps1"
)

scheduler_state_path <- file.path(
  scheduler_dir,
  "V258_scheduler_state.rds"
)

scheduler_lock_dir <- file.path(
  scheduler_dir,
  "V258_scheduler_state.lock"
)

scheduler_launch_lock_dir <- file.path(
  scheduler_dir,
  "V258_launch.lock"
)

scheduler_last_launch_path <- file.path(
  scheduler_dir,
  "V258_last_launch.rds"
)

scheduler_event_log <- file.path(
  output_dir,
  "V258_scheduler_events.csv"
)

scheduler_memory_log <- file.path(
  output_dir,
  "V258_memory_gate.log"
)

scheduler_topology_path <- file.path(
  output_dir,
  "V258_cpu_topology.csv"
)

scheduler_slot_path <- file.path(
  output_dir,
  "V258_cpu_slots.csv"
)

# A small PowerShell helper is generated at runtime so the R script remains
# a single self-contained user-delivered file and requires no Rtools/Rcpp.
ps_helper_lines <- c(
  'param(',
  '  [Parameter(Mandatory=$true)][ValidateSet("topology","affinity","memory","process")][string]$Mode,',
  '  [int]$TargetPid = 0,',
  '  [int]$LogicalIndex = -1,',
  '  [string]$PriorityClass = "AboveNormal"',
  ')',
  '',
  '$ErrorActionPreference = "Stop"',
  '',
  'if ($Mode -eq "process") {',
  '  if ($TargetPid -le 0) { throw "TargetPid must be positive." }',
  '  $p = Get-Process -Id $TargetPid -ErrorAction SilentlyContinue',
  '  if ($null -eq $p) {',
  '    [PSCustomObject]@{',
  '      PID = $TargetPid',
  '      ProcessName = ""',
  '      ExecutablePath = ""',
  '      StartTimeLocal = ""',
  '      CommandLine = ""',
  '      Status = "NOT_FOUND"',
  '    } | ConvertTo-Csv -NoTypeInformation',
  '  } else {',
  '    $exePath = ""',
  '    $startLocal = ""',
  '    $commandLine = ""',
  '    try {',
  '      if ($null -ne $p.Path) { $exePath = [string]$p.Path }',
  '    } catch {}',
  '    try {',
  '      $startLocal = $p.StartTime.ToString("yyyy-MM-dd HH:mm:ss.fff")',
  '    } catch {}',
  '    try {',
  '      $cim = Get-CimInstance Win32_Process -Filter ("ProcessId = " + $TargetPid) -ErrorAction SilentlyContinue',
  '      if ($null -ne $cim) {',
  '        if ($cim.ExecutablePath) { $exePath = [string]$cim.ExecutablePath }',
  '        if ($cim.CommandLine) { $commandLine = [string]$cim.CommandLine }',
  '      }',
  '    } catch {}',
  '    [PSCustomObject]@{',
  '      PID = $TargetPid',
  '      ProcessName = $p.ProcessName',
  '      ExecutablePath = $exePath',
  '      StartTimeLocal = $startLocal',
  '      CommandLine = $commandLine',
  '      Status = "FOUND"',
  '    } | ConvertTo-Csv -NoTypeInformation',
  '  }',
  '  exit 0',
  '}',
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
  'if ($Mode -eq "affinity") {',
  '  if ($TargetPid -le 0) { throw "TargetPid must be positive." }',
  '  if ($LogicalIndex -lt 0 -or $LogicalIndex -gt 62) { throw "LogicalIndex must be between 0 and 62." }',
  '  $p = Get-Process -Id $TargetPid',
  '  $mask = ([Int64]1 -shl $LogicalIndex)',
  '  $p.ProcessorAffinity = [IntPtr]$mask',
  '  if ($PriorityClass -and $PriorityClass -ne "Normal") {',
  '    $p.PriorityClass = $PriorityClass',
  '  }',
  '  [PSCustomObject]@{',
  '    PID = $TargetPid',
  '    LogicalIndex = $LogicalIndex',
  '    AffinityMask = $mask',
  '    PriorityClass = $p.PriorityClass.ToString()',
  '    Status = "OK"',
  '  } | ConvertTo-Csv -NoTypeInformation',
  '  exit 0',
  '}',
  '',
  '$src = @"',
  'using System;',
  'using System.Runtime.InteropServices;',
  'public static class V258CpuSetNative {',
  '  [DllImport("kernel32.dll", SetLastError=true)]',
  '  public static extern bool GetSystemCpuSetInformation(',
  '    IntPtr Information,',
  '    uint BufferLength,',
  '    out uint ReturnedLength,',
  '    IntPtr Process,',
  '    uint Flags);',
  '}',
  '"@',
  '',
  'Add-Type -TypeDefinition $src',
  '$needed = [uint32]0',
  '[void][V258CpuSetNative]::GetSystemCpuSetInformation([IntPtr]::Zero, 0, [ref]$needed, [IntPtr]::Zero, 0)',
  'if ($needed -le 0) { throw "GetSystemCpuSetInformation returned zero required bytes." }',
  '$buf = [Runtime.InteropServices.Marshal]::AllocHGlobal([int]$needed)',
  'try {',
  '  $returned = [uint32]0',
  '  $ok = [V258CpuSetNative]::GetSystemCpuSetInformation($buf, $needed, [ref]$returned, [IntPtr]::Zero, 0)',
  '  if (-not $ok) { throw "GetSystemCpuSetInformation failed." }',
  '  $cpuName = (Get-CimInstance Win32_Processor | Select-Object -First 1 -ExpandProperty Name)',
  '  $offset = 0',
  '  $rows = @()',
  '  while ($offset -lt $returned) {',
  '    $ptr = [IntPtr]::Add($buf, $offset)',
  '    $size = [Runtime.InteropServices.Marshal]::ReadInt32($ptr, 0)',
  '    $type = [Runtime.InteropServices.Marshal]::ReadInt32($ptr, 4)',
  '    if ($size -le 0) { throw "Invalid CPU set structure size." }',
  '    if ($type -eq 0) {',
  '      $id = [Runtime.InteropServices.Marshal]::ReadInt32($ptr, 8)',
  '      $group = [Runtime.InteropServices.Marshal]::ReadInt16($ptr, 12)',
  '      $logical = [Runtime.InteropServices.Marshal]::ReadByte($ptr, 14)',
  '      $core = [Runtime.InteropServices.Marshal]::ReadByte($ptr, 15)',
  '      $llc = [Runtime.InteropServices.Marshal]::ReadByte($ptr, 16)',
  '      $numa = [Runtime.InteropServices.Marshal]::ReadByte($ptr, 17)',
  '      $eff = [Runtime.InteropServices.Marshal]::ReadByte($ptr, 18)',
  '      $flags = [Runtime.InteropServices.Marshal]::ReadByte($ptr, 19)',
  '      $sched = [Runtime.InteropServices.Marshal]::ReadByte($ptr, 20)',
  '      $rows += [PSCustomObject]@{',
  '        ProcessorName = $cpuName',
  '        CpuSetId = [uint32]$id',
  '        Group = [int]$group',
  '        LogicalProcessorIndex = [int]$logical',
  '        CoreIndex = [int]$core',
  '        LastLevelCacheIndex = [int]$llc',
  '        NumaNodeIndex = [int]$numa',
  '        EfficiencyClass = [int]$eff',
  '        SchedulingClass = [int]$sched',
  '        AllFlags = [int]$flags',
  '      }',
  '    }',
  '    $offset += $size',
  '  }',
  '  $rows | Sort-Object Group,LogicalProcessorIndex | ConvertTo-Csv -NoTypeInformation',
  '} finally {',
  '  [Runtime.InteropServices.Marshal]::FreeHGlobal($buf)',
  '}'
)

writeLines(
  ps_helper_lines,
  scheduler_helper_ps1,
  useBytes = TRUE
)

run_ps_helper <- function(
  mode,
  target_pid = 0L,
  logical_index = -1L,
  priority_class = PROCESS_PRIORITY_CLASS,
  helper_path = scheduler_helper_ps1
) {

  args <- c(
    "-NoProfile",
    "-ExecutionPolicy",
    "Bypass",
    "-File",
    shQuote(helper_path),
    "-Mode",
    mode
  )

  if (identical(mode, "affinity")) {
    args <- c(
      args,
      "-TargetPid",
      as.character(as.integer(target_pid)),
      "-LogicalIndex",
      as.character(as.integer(logical_index)),
      "-PriorityClass",
      priority_class
    )
  } else if (identical(mode, "process")) {
    args <- c(
      args,
      "-TargetPid",
      as.character(as.integer(target_pid))
    )
  }

  out <- tryCatch(
    system2(
      powershell_exe,
      args = args,
      stdout = TRUE,
      stderr = TRUE
    ),
    error = function(e) e
  )

  if (inherits(out, "error")) {
    stop(
      paste0(
        "PowerShell helper failed in mode '",
        mode,
        "': ",
        conditionMessage(out)
      )
    )
  }

  status <- attr(out, "status")

  if (!is.null(status) && as.integer(status) != 0L) {
    stop(
      paste0(
        "PowerShell helper returned status ",
        status,
        " in mode '",
        mode,
        "'.\nOutput:\n",
        paste(out, collapse = "\n")
      )
    )
  }

  out
}

read_ps_csv <- function(lines) {

  if (length(lines) == 0L) {
    stop("PowerShell helper returned no CSV output.")
  }

  txt <- paste(lines, collapse = "\n")

  out <- tryCatch(
    data.table::fread(
      text = txt,
      showProgress = FALSE
    ),
    error = function(e) e
  )

  if (inherits(out, "error")) {
    stop(
      paste0(
        "Could not parse PowerShell CSV output: ",
        conditionMessage(out),
        "\nRaw output:\n",
        txt
      )
    )
  }

  out
}

query_cpu_topology <- function() {

  top <- read_ps_csv(
    run_ps_helper("topology")
  )

  required <- c(
    "ProcessorName",
    "CpuSetId",
    "Group",
    "LogicalProcessorIndex",
    "CoreIndex",
    "LastLevelCacheIndex",
    "NumaNodeIndex",
    "EfficiencyClass",
    "SchedulingClass",
    "AllFlags"
  )

  missing <- setdiff(required, names(top))

  if (length(missing) > 0L) {
    stop(
      paste0(
        "CPU topology helper is missing field(s): ",
        paste(missing, collapse = ", ")
      )
    )
  }

  for (v in setdiff(required, "ProcessorName")) {
    top[[v]] <- as.integer(top[[v]])
  }

  if (any(top$Group != 0L)) {
    stop(
      paste0(
        "V258 currently requires all CPUs to be in Windows processor group 0.\n",
        "Detected group(s): ",
        paste(sort(unique(top$Group)), collapse = ", "),
        "\nThe project machines are expected to fit within one processor group."
      )
    )
  }

  if (max(top$LogicalProcessorIndex) > 62L) {
    stop(
      "V258 ProcessorAffinity helper currently supports logical CPU indices 0-62."
    )
  }

  top
}

classify_cpu_slots <- function(top) {

  dt <- data.table::as.data.table(top)
  dt[, core_key := paste(Group, CoreIndex, sep = "_")]

  cores <- dt[
    ,
    .(
      ProcessorName = ProcessorName[1L],
      Group = Group[1L],
      CoreIndex = CoreIndex[1L],
      n_logical = .N,
      logical_indices = list(sort(LogicalProcessorIndex)),
      EfficiencyClass = max(EfficiencyClass),
      SchedulingClass = max(SchedulingClass),
      LastLevelCacheIndex = LastLevelCacheIndex[1L]
    ),
    by = core_key
  ]

  n_physical <- nrow(cores)
  n_logical <- nrow(dt)
  hybrid <- any(cores$n_logical > 1L) && any(cores$n_logical == 1L)

  cores[, core_type := NA_character_]
  classification_note <- ""

  if (!hybrid) {

    cores[, core_type := "P"]
    classification_note <- paste0(
      "Homogeneous topology detected: all physical cores share one priority tier; ",
      "primary hardware threads are treated as physical-core slots and siblings as SMT."
    )

  } else {

    cores[n_logical > 1L, core_type := "P"]

    singleton_idx <- which(cores$n_logical == 1L)
    single <- cores[singleton_idx]

    eff_values <- sort(unique(single$EfficiencyClass), decreasing = TRUE)
    sched_values <- sort(unique(single$SchedulingClass), decreasing = TRUE)

    cpu_name <- paste(unique(dt$ProcessorName), collapse = " ")

    if (length(eff_values) >= 2L) {

      best_eff <- max(eff_values)
      cores[singleton_idx[single$EfficiencyClass == best_eff], core_type := "E"]
      cores[singleton_idx[single$EfficiencyClass < best_eff], core_type := "LP-E"]
      classification_note <- "Hybrid singleton cores split by Windows EfficiencyClass."

    } else if (
      grepl("155H", cpu_name, ignore.case = TRUE) &&
      nrow(cores[n_logical > 1L]) == 6L &&
      nrow(single) == 10L &&
      length(sched_values) >= 2L
    ) {

      sched_counts <- single[, .N, by = SchedulingClass][order(SchedulingClass)]
      lp_candidates <- sched_counts[N == 2L]

      if (nrow(lp_candidates) >= 1L) {
        lp_sched <- min(lp_candidates$SchedulingClass)
        cores[
          singleton_idx[single$SchedulingClass == lp_sched],
          core_type := "LP-E"
        ]
        cores[
          singleton_idx[single$SchedulingClass != lp_sched],
          core_type := "E"
        ]
        classification_note <- paste0(
          "Core Ultra 7 155H: E and LP-E share singleton topology; the unique lowest ",
          "2-core Windows SchedulingClass group was classified as LP-E, remaining ",
          "singleton cores as E."
        )
      } else {
        stop(
          paste0(
            "Core Ultra 7 155H topology was detected, but SchedulingClass did not ",
            "contain a 2-core LP-E candidate. Stop rather than guess."
          )
        )
      }

    } else if (length(sched_values) >= 2L) {

      best_sched <- max(sched_values)
      cores[singleton_idx[single$SchedulingClass == best_sched], core_type := "E"]
      cores[singleton_idx[single$SchedulingClass < best_sched], core_type := "LP-E"]
      classification_note <- "Hybrid singleton cores split by Windows SchedulingClass."

    } else {

      # Deterministic fallback for the known project laptop: Core Ultra 7 155H
      # has 6 SMT P cores, 8 E cores and 2 LP-E cores. We still use actual
      # topology metadata rather than guessed logical CPU numbers. If a unique
      # 2-core LLC island exists among the 10 singleton cores, that island is
      # identified as LP-E. Otherwise the script stops rather than guessing.
      llc_sizes <- single[, .N, by = LastLevelCacheIndex]
      two_core_llc <- llc_sizes[N == 2L, LastLevelCacheIndex]

      if (
        grepl("155H", cpu_name, ignore.case = TRUE) &&
        nrow(cores[n_logical > 1L]) == 6L &&
        nrow(single) == 10L &&
        length(two_core_llc) == 1L
      ) {

        lp_llc <- two_core_llc[1L]
        cores[
          singleton_idx[single$LastLevelCacheIndex == lp_llc],
          core_type := "LP-E"
        ]
        cores[
          singleton_idx[single$LastLevelCacheIndex != lp_llc],
          core_type := "E"
        ]
        classification_note <- paste0(
          "Core Ultra 7 155H fallback: actual CPU-set topology found one 2-core ",
          "singleton LLC island, classified as LP-E; remaining singleton cores as E."
        )

      } else {

        fwrite(
          dt,
          scheduler_topology_path
        )

        stop(
          paste0(
            "Windows CPU-set metadata identified P versus singleton cores, but could not ",
            "reliably distinguish E from LP-E.\n",
            "Topology was written to:\n",
            normalizePath(
              scheduler_topology_path,
              winslash = "/",
              mustWork = FALSE
            ),
            "\nStop here rather than guess core types."
          )
        )

      }

    }

  }

  if (anyNA(cores$core_type)) {
    stop("CPU core classification produced NA core_type values.")
  }

  slot_list <- list()
  k <- 0L

  for (i in seq_len(nrow(cores))) {

    logicals <- unlist(cores$logical_indices[[i]])
    logicals <- sort(as.integer(logicals))

    primary <- logicals[1L]
    type <- cores$core_type[i]

    rank <- switch(
      type,
      "P" = 1L,
      "E" = 2L,
      "LP-E" = 3L,
      stop("Unknown core_type during slot construction.")
    )

    k <- k + 1L
    slot_list[[k]] <- data.table(
      slot_id = paste0(
        "G",
        cores$Group[i],
        "_C",
        cores$CoreIndex[i],
        "_L",
        primary
      ),
      Group = as.integer(cores$Group[i]),
      CoreIndex = as.integer(cores$CoreIndex[i]),
      LogicalProcessorIndex = as.integer(primary),
      tier_name = type,
      tier_rank = as.integer(rank),
      is_smt = 0L
    )

    if (length(logicals) > 1L) {

      for (sib in logicals[-1L]) {
        k <- k + 1L
        slot_list[[k]] <- data.table(
          slot_id = paste0(
            "G",
            cores$Group[i],
            "_C",
            cores$CoreIndex[i],
            "_L",
            sib,
            "_SMT"
          ),
          Group = as.integer(cores$Group[i]),
          CoreIndex = as.integer(cores$CoreIndex[i]),
          LogicalProcessorIndex = as.integer(sib),
          tier_name = "SMT",
          tier_rank = 4L,
          is_smt = 1L
        )
      }

    }

  }

  slots <- data.table::rbindlist(slot_list)
  data.table::setorder(slots, tier_rank, CoreIndex, LogicalProcessorIndex)

  list(
    topology = dt,
    cores = cores,
    slots = slots,
    n_physical = as.integer(n_physical),
    n_logical = as.integer(n_logical),
    hybrid = as.integer(hybrid),
    classification_note = classification_note
  )
}


check_existing_v251_workers <- function(
  output_dir,
  helper_path
) {

  marker_files <- list.files(
    output_dir,
    pattern = "^V251_RUNNING_G[0-9]+\\.txt$",
    full.names = TRUE
  )

  if (length(marker_files) == 0L) {
    return(
      data.table()
    )
  }

  # Build the exact set of R executables from the currently installed R tree.
  # No R process name is guessed or hard-coded.
  r_executable_files <- list.files(
    R.home(),
    pattern = "\\.exe$",
    recursive = TRUE,
    full.names = TRUE,
    ignore.case = TRUE
  )

  if (length(r_executable_files) == 0L) {
    stop(
      paste0(
        "Could not enumerate any R executables under the exact current R.home():\\n",
        normalizePath(
          R.home(),
          winslash = "/",
          mustWork = TRUE
        ),
        "\\nStop here rather than guess process names."
      )
    )
  }

  r_executable_paths <- unique(
    tolower(
      normalizePath(
        r_executable_files,
        winslash = "/",
        mustWork = TRUE
      )
    )
  )

  r_executable_names <- unique(
    tolower(
      tools::file_path_sans_ext(
        basename(
          r_executable_files
        )
      )
    )
  )

  rows <- lapply(
    marker_files,
    function(path) {

      lines <- readLines(
        path,
        warn = FALSE,
        encoding = "UTF-8"
      )

      pid_lines <- grep(
        "^PID: [0-9]+$",
        lines,
        value = TRUE
      )

      if (length(pid_lines) != 1L) {
        stop(
          paste0(
            "Existing RUNNING marker cannot be parsed exactly:\\n",
            normalizePath(
              path,
              winslash = "/",
              mustWork = TRUE
            ),
            "\\nExpected exactly one line of the form 'PID: <integer>'. ",
            "Stop here rather than guess whether another GBMTM worker is active."
          )
        )
      }

      start_lines <- grep(
        "^Start: [0-9]{4}-[0-9]{2}-[0-9]{2} [0-9]{2}:[0-9]{2}:[0-9]{2}$",
        lines,
        value = TRUE
      )

      if (length(start_lines) != 1L) {
        stop(
          paste0(
            "Existing RUNNING marker has no exact task start timestamp:\\n",
            normalizePath(
              path,
              winslash = "/",
              mustWork = TRUE
            ),
            "\\nExpected exactly one line 'Start: YYYY-MM-DD HH:MM:SS'. ",
            "Stop here rather than guess process identity."
          )
        )
      }

      marker_start_text <- sub(
        "^Start: ",
        "",
        start_lines
      )

      marker_start <- as.POSIXct(
        marker_start_text,
        format = "%Y-%m-%d %H:%M:%S"
      )

      if (is.na(marker_start)) {
        stop(
          paste0(
            "Could not parse marker Start timestamp exactly: ",
            marker_start_text,
            "\\nMarker: ",
            normalizePath(
              path,
              winslash = "/",
              mustWork = TRUE
            )
          )
        )
      }

      pid <- as.integer(
        sub(
          "^PID: ",
          "",
          pid_lines
        )
      )

      g_match <- regexec(
        "^V251_RUNNING_G([0-9]+)\\.txt$",
        basename(path)
      )

      g_parts <- regmatches(
        basename(path),
        g_match
      )[[1L]]

      if (length(g_parts) != 2L) {
        stop(
          paste0(
            "Could not parse G from existing RUNNING marker exactly:\\n",
            basename(path)
          )
        )
      }

      g <- as.integer(
        g_parts[2L]
      )

      proc <- read_ps_csv(
        run_ps_helper(
          mode = "process",
          target_pid = pid,
          helper_path = helper_path
        )
      )

      required_proc_fields <- c(
        "PID",
        "ProcessName",
        "ExecutablePath",
        "StartTimeLocal",
        "CommandLine",
        "Status"
      )

      if (
        nrow(proc) != 1L ||
        !all(
          required_proc_fields %in%
            names(proc)
        )
      ) {
        stop(
          paste0(
            "Process check returned an unexpected structure for PID ",
            pid,
            ". Stop here rather than guess."
          )
        )
      }

      process_status <- as.character(
        proc$Status[1L]
      )

      process_name <- as.character(
        proc$ProcessName[1L]
      )

      executable_path <- as.character(
        proc$ExecutablePath[1L]
      )

      process_start_text <- as.character(
        proc$StartTimeLocal[1L]
      )

      command_line <- as.character(
        proc$CommandLine[1L]
      )

      if (is.na(process_name)) process_name <- ""
      if (is.na(executable_path)) executable_path <- ""
      if (is.na(process_start_text)) process_start_text <- ""
      if (is.na(command_line)) command_line <- ""

      identity_status <- ""
      identity_reason <- ""
      process_start <- as.POSIXct(NA)
      exact_r_executable <- FALSE

      if (identical(process_status, "NOT_FOUND")) {

        identity_status <- "STALE"
        identity_reason <- "PID_NOT_FOUND"

      } else if (!identical(process_status, "FOUND")) {

        stop(
          paste0(
            "Unexpected process Status for PID ",
            pid,
            ": ",
            process_status,
            ". Stop here rather than guess."
          )
        )

      } else {

        process_name_lower <- tolower(
          process_name
        )

        process_name_is_installed_r_executable <- (
          nzchar(process_name_lower) &&
            process_name_lower %in%
              r_executable_names
        )

        if (nzchar(executable_path)) {

          executable_path_normalized <- tolower(
            normalizePath(
              executable_path,
              winslash = "/",
              mustWork = FALSE
            )
          )

          exact_r_executable <- executable_path_normalized %in%
            r_executable_paths

          if (
            !exact_r_executable &&
            process_name_is_installed_r_executable
          ) {
            stop(
              paste0(
                "PID ",
                pid,
                " is currently an R-named executable, but its exact executable path is outside current R.home().\\n",
                "Process: ",
                process_name,
                "\\nExecutable: ",
                executable_path,
                "\\nMarker: ",
                normalizePath(
                  path,
                  winslash = "/",
                  mustWork = TRUE
                ),
                "\\nStop here rather than assume it is stale or active."
              )
            )
          }

        } else {

          # When Windows does not expose ExecutablePath, only accept the process
          # as an R process if its exact ProcessName is present among executables
          # enumerated from the current R installation.
          exact_r_executable <- process_name_is_installed_r_executable

        }

        if (!exact_r_executable) {

          identity_status <- "STALE"
          identity_reason <- "PID_REUSED_BY_NON_R_PROCESS"

        } else {

          if (!nzchar(process_start_text)) {
            stop(
              paste0(
                "PID ",
                pid,
                " is an exact current-R executable, but Windows did not expose its process start time.\\n",
                "Marker: ",
                normalizePath(
                  path,
                  winslash = "/",
                  mustWork = TRUE
                ),
                "\\nStop here rather than guess whether this PID was reused."
              )
            )
          }

          process_start <- as.POSIXct(
            process_start_text,
            format = "%Y-%m-%d %H:%M:%OS"
          )

          if (is.na(process_start)) {
            stop(
              paste0(
                "Could not parse current process start time exactly for PID ",
                pid,
                ": ",
                process_start_text
              )
            )
          }

          # A genuine worker process already exists before it writes its own
          # RUNNING marker. Therefore a process whose creation time is later
          # than the marker's task Start timestamp is a reused PID, not the
          # original worker. Two seconds only absorbs timestamp precision.
          if (
            as.numeric(
              difftime(
                process_start,
                marker_start,
                units = "secs"
              )
            ) > 2
          ) {

            identity_status <- "STALE"
            identity_reason <- "PID_REUSED_AFTER_MARKER_START"

          } else {

            marker_worker_start_lines <- grep(
              "^Worker process start local: .+$",
              lines,
              value = TRUE
            )

            marker_worker_exe_lines <- grep(
              "^Worker executable: .+$",
              lines,
              value = TRUE
            )

            if (
              length(marker_worker_start_lines) > 1L ||
              length(marker_worker_exe_lines) > 1L
            ) {
              stop(
                paste0(
                  "RUNNING marker contains duplicate V258 process-identity lines:\\n",
                  normalizePath(
                    path,
                    winslash = "/",
                    mustWork = TRUE
                  )
                )
              )
            }

            if (length(marker_worker_start_lines) == 1L) {

              recorded_start_text <- sub(
                "^Worker process start local: ",
                "",
                marker_worker_start_lines
              )

              recorded_start <- as.POSIXct(
                recorded_start_text,
                format = "%Y-%m-%d %H:%M:%OS"
              )

              if (is.na(recorded_start)) {
                stop(
                  paste0(
                    "Could not parse recorded worker process start time in marker:\\n",
                    normalizePath(
                      path,
                      winslash = "/",
                      mustWork = TRUE
                    )
                  )
                )
              }

              if (
                abs(
                  as.numeric(
                    difftime(
                      process_start,
                      recorded_start,
                      units = "secs"
                    )
                  )
                ) > 1
              ) {
                identity_status <- "STALE"
                identity_reason <- "PID_PROCESS_START_DOES_NOT_MATCH_MARKER"
              }

            }

            if (
              identity_status == "" &&
              length(marker_worker_exe_lines) == 1L &&
              nzchar(executable_path)
            ) {

              recorded_executable <- sub(
                "^Worker executable: ",
                "",
                marker_worker_exe_lines
              )

              current_executable_normalized <- tolower(
                normalizePath(
                  executable_path,
                  winslash = "/",
                  mustWork = FALSE
                )
              )

              recorded_executable_normalized <- tolower(
                normalizePath(
                  recorded_executable,
                  winslash = "/",
                  mustWork = FALSE
                )
              )

              if (!identical(
                current_executable_normalized,
                recorded_executable_normalized
              )) {
                identity_status <- "STALE"
                identity_reason <- "PID_EXECUTABLE_DOES_NOT_MATCH_MARKER"
              }

            }

            if (identity_status == "") {
              identity_status <- "ACTIVE"
              identity_reason <- "EXACT_R_PROCESS_IDENTITY_CONFIRMED"
            }

          }

        }

      }

      data.table(
        requested_g = as.integer(g),
        pid = as.integer(pid),
        process_name = process_name,
        executable_path = executable_path,
        process_start_local = process_start_text,
        process_status = process_status,
        identity_status = identity_status,
        identity_reason = identity_reason,
        marker_start = marker_start_text,
        marker_file = normalizePath(
          path,
          winslash = "/",
          mustWork = TRUE
        )
      )
    }
  )

  out <- rbindlist(
    rows,
    fill = TRUE
  )

  cleanup_audit_path <- file.path(
    output_dir,
    "V258_stale_marker_cleanup.csv"
  )

  stale <- out[
    identity_status == "STALE"
  ]

  if (nrow(stale) > 0L) {

    stale_audit <- copy(stale)
    stale_audit[, removed_at := format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )]

    fwrite(
      stale_audit,
      cleanup_audit_path,
      append = file.exists(cleanup_audit_path),
      col.names = !file.exists(cleanup_audit_path)
    )

    for (f in stale$marker_file) {
      if (file.exists(f)) {
        file.remove(f)
      }
    }

  }

  active <- out[
    identity_status == "ACTIVE"
  ]

  if (nrow(active) > 0L) {
    stop(
      paste0(
        "Another V251-compatible GBMTM worker process is still active.\\n",
        "Do NOT run V258 concurrently with V251/V254/V258 workers that write to the same folder.\\n",
        paste0(
          "G",
          active$requested_g,
          " | PID ",
          active$pid,
          " | process ",
          active$process_name,
          " | executable ",
          active$executable_path,
          " | process start ",
          active$process_start_local,
          " | marker start ",
          active$marker_start,
          " | marker ",
          active$marker_file,
          collapse = "\\n"
        ),
        "\\nStop the existing run first, then start V258."
      )
    )
  }

  invisible(out)
}

existing_worker_check <- check_existing_v251_workers(
  output_dir = output_dir,
  helper_path = scheduler_helper_ps1
)

cpu_info <- classify_cpu_slots(
  query_cpu_topology()
)

cpu_topology <- cpu_info$topology
cpu_slots <- cpu_info$slots
physical_cores <- cpu_info$n_physical
logical_cores <- cpu_info$n_logical
MAX_ACTIVE_TASKS <- logical_cores

fwrite(
  cpu_topology,
  scheduler_topology_path
)

fwrite(
  cpu_slots,
  scheduler_slot_path
)

# Scheduler state is session-local. Model checkpoints are persistent.
# On each new master session we intentionally rebuild CPU leases because
# worker PIDs change after reboot/restart.
initial_slots <- data.table::copy(cpu_slots)
initial_slots[, `:=`(
  assigned_g = NA_integer_,
  assigned_pid = NA_integer_,
  reserved_g = NA_integer_,
  reserved_pid = NA_integer_
)]

initial_tasks <- data.table(
  g = integer(),
  pid = integer(),
  status = character(),
  slot_id = character(),
  reserved_slot_id = character(),
  iteration = integer(),
  recent_delta_mean = numeric(),
  avg_iter_seconds = numeric(),
  updated_at = character()
)

scheduler_state <- list(
  engine = ENGINE_TAG,
  created_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
  classification_note = cpu_info$classification_note,
  slots = initial_slots,
  tasks = initial_tasks,
  delta_history = list(),
  time_history = list()
)

saveRDS(
  scheduler_state,
  scheduler_state_path,
  compress = FALSE
)

if (dir.exists(scheduler_lock_dir)) {
  unlink(scheduler_lock_dir, recursive = TRUE, force = TRUE)
}

if (dir.exists(scheduler_launch_lock_dir)) {
  unlink(scheduler_launch_lock_dir, recursive = TRUE, force = TRUE)
}

if (file.exists(scheduler_last_launch_path)) {
  file.remove(scheduler_last_launch_path)
}

if (file.exists(scheduler_event_log)) {
  file.remove(scheduler_event_log)
}

scheduler_acquire_lock <- function(
  lock_dir,
  timeout_seconds = 120
) {

  start <- Sys.time()

  repeat {

    ok <- dir.create(
      lock_dir,
      showWarnings = FALSE,
      recursive = FALSE
    )

    if (isTRUE(ok)) {
      return(invisible(TRUE))
    }

    elapsed <- as.numeric(
      difftime(
        Sys.time(),
        start,
        units = "secs"
      )
    )

    if (elapsed > timeout_seconds) {
      stop(
        paste0(
          "Timed out waiting for scheduler lock:\n",
          lock_dir
        )
      )
    }

    Sys.sleep(runif(1L, 0.05, 0.15))

  }
}

scheduler_release_lock <- function(lock_dir) {
  if (dir.exists(lock_dir)) {
    unlink(lock_dir, recursive = TRUE, force = TRUE)
  }
  invisible(TRUE)
}

scheduler_save_state_locked <- function(
  state,
  state_path
) {

  tmp <- paste0(
    state_path,
    ".tmp_",
    Sys.getpid()
  )

  if (file.exists(tmp)) {
    file.remove(tmp)
  }

  saveRDS(
    state,
    tmp,
    compress = FALSE
  )

  test <- tryCatch(
    readRDS(tmp),
    error = function(e) e
  )

  if (inherits(test, "error")) {
    stop(
      paste0(
        "Scheduler-state verification failed: ",
        conditionMessage(test)
      )
    )
  }

  rm(test)

  if (file.exists(state_path)) {
    file.remove(state_path)
  }

  if (!file.rename(tmp, state_path)) {
    if (!file.copy(tmp, state_path, overwrite = TRUE)) {
      stop("Could not finalize scheduler state file.")
    }
    file.remove(tmp)
  }

  invisible(TRUE)
}

scheduler_log_event_locked <- function(
  event_log,
  event,
  g = NA_integer_,
  pid = NA_integer_,
  from_slot = NA_character_,
  to_slot = NA_character_,
  iteration = NA_integer_,
  detail = ""
) {

  row <- data.table(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    event = as.character(event),
    g = as.integer(g),
    pid = as.integer(pid),
    from_slot = as.character(from_slot),
    to_slot = as.character(to_slot),
    iteration = as.integer(iteration),
    detail = as.character(detail)
  )

  fwrite(
    row,
    event_log,
    append = file.exists(event_log),
    col.names = !file.exists(event_log)
  )

  invisible(TRUE)
}

scheduler_apply_affinity <- function(
  pid,
  logical_index,
  helper_path,
  priority_class
) {

  out <- run_ps_helper(
    mode = "affinity",
    target_pid = as.integer(pid),
    logical_index = as.integer(logical_index),
    priority_class = priority_class,
    helper_path = helper_path
  )

  parsed <- read_ps_csv(out)

  if (
    nrow(parsed) < 1L ||
    !"Status" %in% names(parsed) ||
    !identical(as.character(parsed$Status[1L]), "OK")
  ) {
    stop(
      paste0(
        "Affinity assignment did not return OK for PID ",
        pid,
        " to logical CPU ",
        logical_index,
        "."
      )
    )
  }

  invisible(TRUE)
}

query_memory_status <- function(helper_path) {

  mem <- read_ps_csv(
    run_ps_helper(
      mode = "memory",
      helper_path = helper_path
    )
  )

  required <- c(
    "TotalPhysicalMB",
    "AvailablePhysicalMB",
    "CommitLimitMB",
    "CommitUsedMB"
  )

  missing <- setdiff(required, names(mem))

  if (length(missing) > 0L) {
    stop(
      paste0(
        "Memory helper missing field(s): ",
        paste(missing, collapse = ", ")
      )
    )
  }

  for (v in required) {
    mem[[v]] <- as.numeric(mem[[v]])
  }

  mem[1L]
}

scheduler_wait_for_memory <- function(
  g,
  helper_path,
  memory_log,
  min_available_gb,
  min_available_fraction,
  max_commit_fraction,
  poll_seconds
) {

  repeat {

    mem <- query_memory_status(helper_path)

    required_available_mb <- max(
      as.numeric(min_available_gb) * 1024,
      as.numeric(min_available_fraction) * mem$TotalPhysicalMB
    )

    commit_fraction <- if (
      is.finite(mem$CommitLimitMB) &&
      mem$CommitLimitMB > 0
    ) {
      mem$CommitUsedMB / mem$CommitLimitMB
    } else {
      NA_real_
    }

    ok_available <- mem$AvailablePhysicalMB >= required_available_mb
    ok_commit <- is.na(commit_fraction) || commit_fraction <= max_commit_fraction

    if (isTRUE(ok_available) && isTRUE(ok_commit)) {
      return(invisible(mem))
    }

    cat(
      paste0(
        format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
        " | G",
        g,
        " launch delayed by RAM guard; available=",
        round(mem$AvailablePhysicalMB / 1024, 2),
        " GB; required>=",
        round(required_available_mb / 1024, 2),
        " GB; commit=",
        if (is.na(commit_fraction)) "NA" else paste0(round(100 * commit_fraction, 1), "%"),
        "\n"
      ),
      file = memory_log,
      append = TRUE
    )

    Sys.sleep(as.numeric(poll_seconds))

  }
}

scheduler_wait_for_launch_slot <- function(
  g,
  spacing_seconds,
  launch_lock_dir,
  last_launch_path
) {

  scheduler_acquire_lock(
    launch_lock_dir,
    timeout_seconds = 600
  )

  on.exit(
    scheduler_release_lock(launch_lock_dir),
    add = TRUE
  )

  now <- Sys.time()
  wait_seconds <- 0

  if (file.exists(last_launch_path)) {

    last <- tryCatch(
      readRDS(last_launch_path),
      error = function(e) NULL
    )

    if (!is.null(last)) {
      elapsed <- as.numeric(
        difftime(
          now,
          last,
          units = "secs"
        )
      )
      wait_seconds <- max(
        0,
        as.numeric(spacing_seconds) - elapsed
      )
    }

  }

  if (wait_seconds > 0) {
    Sys.sleep(wait_seconds)
  }

  actual <- Sys.time()

  saveRDS(
    actual,
    last_launch_path,
    compress = FALSE
  )

  invisible(
    list(
      g = as.integer(g),
      waited_seconds = as.numeric(wait_seconds),
      actual_start = actual
    )
  )
}

scheduler_reserve_free_slots_locked <- function(
  state,
  event_log
) {

  repeat {

    slots <- data.table::as.data.table(state$slots)
    tasks <- data.table::as.data.table(state$tasks)

    free <- slots[
      is.na(assigned_g) &
      is.na(reserved_g)
    ]

    if (nrow(free) == 0L || nrow(tasks) == 0L) {
      break
    }

    data.table::setorder(
      free,
      tier_rank,
      CoreIndex,
      LogicalProcessorIndex
    )

    made_reservation <- FALSE

    for (i in seq_len(nrow(free))) {

      free_slot <- free[i]

      running <- tasks[
        status == "RUNNING" &
        (is.na(reserved_slot_id) | reserved_slot_id == "")
      ]

      if (nrow(running) == 0L) {
        next
      }

      running <- merge(
        running,
        slots[
          ,
          .(
            slot_id,
            current_tier_rank = tier_rank
          )
        ],
        by = "slot_id",
        all.x = TRUE,
        sort = FALSE
      )

      candidates <- running[
        current_tier_rank > free_slot$tier_rank
      ]

      if (nrow(candidates) == 0L) {
        next
      }

      candidates[
        ,
        iteration_sort := fifelse(
          is.na(iteration),
          0L,
          iteration
        )
      ]

      candidates[
        ,
        delta_sort := fifelse(
          is.na(recent_delta_mean),
          -Inf,
          recent_delta_mean
        )
      ]

      candidates[
        ,
        time_sort := fifelse(
          is.na(avg_iter_seconds),
          -Inf,
          avg_iter_seconds
        )
      ]

      data.table::setorder(
        candidates,
        iteration_sort,
        -delta_sort,
        -time_sort,
        g
      )

      target <- candidates[1L]

      slot_idx <- which(
        state$slots$slot_id == free_slot$slot_id
      )

      task_idx <- which(
        state$tasks$g == target$g &
        state$tasks$pid == target$pid
      )

      if (length(slot_idx) != 1L || length(task_idx) != 1L) {
        stop("Scheduler reservation index mismatch.")
      }

      state$slots$reserved_g[slot_idx] <- as.integer(target$g)
      state$slots$reserved_pid[slot_idx] <- as.integer(target$pid)
      state$tasks$reserved_slot_id[task_idx] <- as.character(free_slot$slot_id)

      scheduler_log_event_locked(
        event_log = event_log,
        event = "RESERVE_PROMOTION",
        g = target$g,
        pid = target$pid,
        from_slot = target$slot_id,
        to_slot = free_slot$slot_id,
        iteration = target$iteration,
        detail = paste0(
          "free tier=",
          free_slot$tier_name,
          "; candidate chosen by low iteration -> high recent deltaLL -> long iteration time"
        )
      )

      made_reservation <- TRUE
      break

    }

    if (!made_reservation) {
      break
    }

  }

  state
}

scheduler_register_task <- function(
  g,
  pid,
  state_path,
  lock_dir,
  event_log,
  helper_path,
  priority_class
) {

  scheduler_acquire_lock(lock_dir)

  assignment <- NULL

  tryCatch(
    {
      state <- readRDS(state_path)
      slots <- data.table::as.data.table(state$slots)

      free <- slots[
        is.na(assigned_g) &
        is.na(reserved_g)
      ]

      if (nrow(free) == 0L) {
        stop(
          paste0(
            "No CPU slot is free for G",
            g,
            ". Active worker count exceeded detected logical CPUs."
          )
        )
      }

      data.table::setorder(
        free,
        tier_rank,
        CoreIndex,
        LogicalProcessorIndex
      )

      chosen <- free[1L]
      slot_idx <- which(state$slots$slot_id == chosen$slot_id)

      state$slots$assigned_g[slot_idx] <- as.integer(g)
      state$slots$assigned_pid[slot_idx] <- as.integer(pid)

      new_task <- data.table(
        g = as.integer(g),
        pid = as.integer(pid),
        status = "RUNNING",
        slot_id = as.character(chosen$slot_id),
        reserved_slot_id = NA_character_,
        iteration = 0L,
        recent_delta_mean = NA_real_,
        avg_iter_seconds = NA_real_,
        updated_at = format(Sys.time(), "%Y-%m-%d %H:%M:%S")
      )

      g_value <- as.integer(g)
      pid_value <- as.integer(pid)

      state$tasks <- data.table::rbindlist(
        list(
          data.table::as.data.table(state$tasks)[
            !(g == g_value & pid == pid_value)
          ],
          new_task
        ),
        fill = TRUE,
        use.names = TRUE
      )

      scheduler_log_event_locked(
        event_log = event_log,
        event = "ASSIGN_INITIAL",
        g = g,
        pid = pid,
        to_slot = chosen$slot_id,
        iteration = 0L,
        detail = paste0(
          "tier=",
          chosen$tier_name,
          "; logicalCPU=",
          chosen$LogicalProcessorIndex
        )
      )

      scheduler_save_state_locked(state, state_path)

      assignment <- list(
        slot_id = as.character(chosen$slot_id),
        logical_index = as.integer(chosen$LogicalProcessorIndex),
        tier_name = as.character(chosen$tier_name),
        tier_rank = as.integer(chosen$tier_rank)
      )
    },
    finally = {
      scheduler_release_lock(lock_dir)
    }
  )

  affinity_result <- tryCatch(
    {
      scheduler_apply_affinity(
        pid = pid,
        logical_index = assignment$logical_index,
        helper_path = helper_path,
        priority_class = priority_class
      )
      TRUE
    },
    error = function(e) e
  )

  if (inherits(affinity_result, "error")) {
    try(
      scheduler_release_task(
        g = g,
        pid = pid,
        state_path = state_path,
        lock_dir = lock_dir,
        event_log = event_log
      ),
      silent = TRUE
    )
    stop(
      paste0(
        "Initial CPU affinity assignment failed for G",
        g,
        ": ",
        conditionMessage(affinity_result)
      )
    )
  }

  assignment
}

scheduler_iteration_boundary <- function(
  g,
  pid,
  iteration,
  delta_ll,
  iter_seconds,
  state_path,
  lock_dir,
  event_log,
  helper_path,
  priority_class
) {

  scheduler_acquire_lock(lock_dir)

  new_assignment <- NULL

  tryCatch(
    {
      state <- readRDS(state_path)
      task_idx <- which(
        state$tasks$g == as.integer(g) &
        state$tasks$pid == as.integer(pid) &
        state$tasks$status == "RUNNING"
      )

      if (length(task_idx) != 1L) {
        stop(
          paste0(
            "Scheduler cannot find exactly one RUNNING task row for G",
            g,
            " PID ",
            pid,
            "."
          )
        )
      }

      key <- as.character(g)

      dh <- state$delta_history[[key]]
      th <- state$time_history[[key]]

      if (is.null(dh)) dh <- numeric()
      if (is.null(th)) th <- numeric()

      if (is.finite(delta_ll)) {
        dh <- c(dh, as.numeric(delta_ll))
        if (length(dh) > 5L) dh <- tail(dh, 5L)
      }

      if (is.finite(iter_seconds) && iter_seconds >= 0) {
        th <- c(th, as.numeric(iter_seconds))
        if (length(th) > 5L) th <- tail(th, 5L)
      }

      state$delta_history[[key]] <- dh
      state$time_history[[key]] <- th

      state$tasks$iteration[task_idx] <- as.integer(iteration)
      state$tasks$recent_delta_mean[task_idx] <- if (length(dh) > 0L) mean(dh) else NA_real_
      state$tasks$avg_iter_seconds[task_idx] <- if (length(th) > 0L) mean(th) else NA_real_
      state$tasks$updated_at[task_idx] <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

      # First reserve any currently free better slots. A reservation may be
      # made for this task itself because it is known to be at a safe boundary.
      state <- scheduler_reserve_free_slots_locked(
        state,
        event_log
      )

      task_idx <- which(
        state$tasks$g == as.integer(g) &
        state$tasks$pid == as.integer(pid) &
        state$tasks$status == "RUNNING"
      )

      reserved <- state$tasks$reserved_slot_id[task_idx]

      if (!is.na(reserved) && nzchar(reserved)) {

        old_slot <- as.character(state$tasks$slot_id[task_idx])
        old_idx <- which(state$slots$slot_id == old_slot)
        new_idx <- which(state$slots$slot_id == reserved)

        if (length(old_idx) != 1L || length(new_idx) != 1L) {
          stop("Scheduler promotion slot lookup failed.")
        }

        if (
          !is.na(state$slots$assigned_g[new_idx]) ||
          !identical(as.integer(state$slots$reserved_g[new_idx]), as.integer(g)) ||
          !identical(as.integer(state$slots$reserved_pid[new_idx]), as.integer(pid))
        ) {
          stop("Reserved promotion slot is not exclusively reserved for this task.")
        }

        state$slots$assigned_g[old_idx] <- NA_integer_
        state$slots$assigned_pid[old_idx] <- NA_integer_

        state$slots$assigned_g[new_idx] <- as.integer(g)
        state$slots$assigned_pid[new_idx] <- as.integer(pid)
        state$slots$reserved_g[new_idx] <- NA_integer_
        state$slots$reserved_pid[new_idx] <- NA_integer_

        state$tasks$slot_id[task_idx] <- reserved
        state$tasks$reserved_slot_id[task_idx] <- NA_character_

        new_assignment <- list(
          slot_id = reserved,
          logical_index = as.integer(state$slots$LogicalProcessorIndex[new_idx]),
          tier_name = as.character(state$slots$tier_name[new_idx]),
          tier_rank = as.integer(state$slots$tier_rank[new_idx])
        )

        scheduler_log_event_locked(
          event_log = event_log,
          event = "PROMOTE_AT_SAFE_BOUNDARY",
          g = g,
          pid = pid,
          from_slot = old_slot,
          to_slot = reserved,
          iteration = iteration,
          detail = paste0(
            "checkpointed iteration boundary; new tier=",
            new_assignment$tier_name,
            "; logicalCPU=",
            new_assignment$logical_index
          )
        )

        # The old slot is now free. Reserve it for the next eligible slower-tier
        # task if doing so improves that task's tier.
        state <- scheduler_reserve_free_slots_locked(
          state,
          event_log
        )
      }

      scheduler_save_state_locked(state, state_path)
    },
    finally = {
      scheduler_release_lock(lock_dir)
    }
  )

  if (!is.null(new_assignment)) {
    scheduler_apply_affinity(
      pid = pid,
      logical_index = new_assignment$logical_index,
      helper_path = helper_path,
      priority_class = priority_class
    )
  }

  invisible(new_assignment)
}

scheduler_release_task <- function(
  g,
  pid,
  state_path,
  lock_dir,
  event_log
) {

  if (!file.exists(state_path)) {
    return(invisible(FALSE))
  }

  scheduler_acquire_lock(lock_dir)

  tryCatch(
    {
      state <- readRDS(state_path)
      task_idx <- which(
        state$tasks$g == as.integer(g) &
        state$tasks$pid == as.integer(pid) &
        state$tasks$status == "RUNNING"
      )

      if (length(task_idx) == 1L) {

        old_slot <- as.character(state$tasks$slot_id[task_idx])
        slot_idx <- which(state$slots$slot_id == old_slot)

        if (length(slot_idx) == 1L) {
          state$slots$assigned_g[slot_idx] <- NA_integer_
          state$slots$assigned_pid[slot_idx] <- NA_integer_
        }

        # Clear any reservation made for this task but not yet consumed.
        reserved_idx <- which(
          state$slots$reserved_g == as.integer(g) &
          state$slots$reserved_pid == as.integer(pid)
        )

        if (length(reserved_idx) > 0L) {
          state$slots$reserved_g[reserved_idx] <- NA_integer_
          state$slots$reserved_pid[reserved_idx] <- NA_integer_
        }

        state$tasks$status[task_idx] <- "FINISHED"
        state$tasks$reserved_slot_id[task_idx] <- NA_character_
        state$tasks$updated_at[task_idx] <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")

        scheduler_log_event_locked(
          event_log = event_log,
          event = "RELEASE_TASK",
          g = g,
          pid = pid,
          from_slot = old_slot,
          iteration = state$tasks$iteration[task_idx],
          detail = "worker task returned; slot released"
        )

        state <- scheduler_reserve_free_slots_locked(
          state,
          event_log
        )

        scheduler_save_state_locked(state, state_path)
      }
    },
    finally = {
      scheduler_release_lock(lock_dir)
    }
  )

  invisible(TRUE)
}

cat(
  "Detected CPU: ",
  paste(unique(cpu_topology$ProcessorName), collapse = " | "),
  "\n",
  sep = ""
)
cat(
  "Physical / logical cores: ",
  physical_cores,
  " / ",
  logical_cores,
  "\n",
  sep = ""
)
cat(
  "CPU slot order: ",
  paste(
    paste0(
      cpu_slots$tier_name,
      "@L",
      cpu_slots$LogicalProcessorIndex
    ),
    collapse = ", "
  ),
  "\n",
  sep = ""
)
cat(
  "Topology classification: ",
  cpu_info$classification_note,
  "\n\n",
  sep = ""
)


# ============================================================
# 5. Console header
# ============================================================

cat("\n")
cat("============================================================\n")
cat("V258 GBMTM G5-G20 DYNAMIC-SCHEDULER CHECKPOINT/RESUME SCREEN\n")
cat("============================================================\n")
cat(
  "Input : ",
  normalizePath(
    input_rds,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)
cat(
  "Output: ",
  normalizePath(
    output_dir,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n",
  sep = ""
)
cat("Requested G            : G5-G20\n")
cat(paste0("Max active tasks        : ", MAX_ACTIVE_TASKS, " (detected logical CPUs; no arbitrary cap)\n"))
cat(paste0("Global task stagger     : ", START_STAGGER_SECONDS, " sec between actual starts\n"))
cat("Screening tolerance    : 10\n")
cat("nstart                  : 1\n")
cat("delete.empty            : FALSE\n")
cat("Wall-clock timeout      : NONE\n")
cat("Practical maxit budget : NONE\n")
cat(
  "Technical maxit value  : ",
  NONBINDING_MAXIT,
  "\n",
  sep = ""
)
cat("Automatic retry        : NONE\n")
cat("EM checkpoint/resume   : ENABLED; every accepted iteration\n")
cat("CPU scheduler          : P > E > LP-E > SMT; safe-boundary promotion\n")
cat("RAM launch guard       : ENABLED; does not kill running tasks\n")
cat("Compatible final files : V251_* in the existing GBMTM_V251_G5_G20_screen folder\n")
cat("============================================================\n\n")


# ============================================================
# 6. Validate full input
# ============================================================

cat("Reading and validating MIMIC long-format input...\n")

full_dt <- as.data.table(
  readRDS(
    input_rds
  )
)

missing_cols <- setdiff(
  REQUIRED_COLS,
  names(full_dt)
)

if (length(missing_cols) > 0L) {
  stop(
    paste0(
      "Missing required column(s): ",
      paste(missing_cols, collapse = ", ")
    )
  )
}

full_dt[
  ,
  unit_id := as.character(unit_id)
]

full_dt[
  ,
  time_model := as.numeric(time_model)
]

for (v in Z_MARKERS) {
  full_dt[
    ,
    (v) := as.numeric(get(v))
  ]
}

n_patients <- uniqueN(
  full_dt$unit_id
)

n_rows <- nrow(
  full_dt
)

if (n_patients != EXPECTED_PATIENTS) {
  stop(
    paste0(
      "Unexpected patient count: ",
      n_patients,
      "; expected ",
      EXPECTED_PATIENTS,
      "."
    )
  )
}

if (n_rows != EXPECTED_LONG_ROWS) {
  stop(
    paste0(
      "Unexpected row count: ",
      n_rows,
      "; expected ",
      EXPECTED_LONG_ROWS,
      "."
    )
  )
}

if (anyNA(full_dt$unit_id)) {
  stop("unit_id contains missing values.")
}

if (anyNA(full_dt$time_model)) {
  stop("time_model contains missing values.")
}

duplicate_time_rows <- full_dt[
  ,
  .N,
  by = .(
    unit_id,
    time_model
  )
][
  N > 1L
]

if (nrow(duplicate_time_rows) > 0L) {
  stop(
    "Duplicate unit_id + time_model rows detected."
  )
}

model_input <- as.data.frame(
  full_dt[
    ,
    ..REQUIRED_COLS
  ]
)

minimal_input_rds <- file.path(
  output_dir,
  "V251_model_input.rds"
)

saveRDS(
  model_input,
  minimal_input_rds,
  compress = FALSE
)

if (!file.exists(minimal_input_rds)) {
  stop(
    "Failed to create V251_model_input.rds."
  )
}

cat(
  "Input validation passed: ",
  n_patients,
  " patients; ",
  n_rows,
  " rows.\n",
  sep = ""
)

cat(
  "Minimal worker input size: ",
  round(
    file.info(minimal_input_rds)$size / 1024^2,
    2
  ),
  " MB\n\n",
  sep = ""
)

input_md5 <- unname(
  tools::md5sum(
    input_rds
  )
)

cat(
  "Source input MD5        : ",
  input_md5,
  "\n\n",
  sep = ""
)

rm(
  full_dt,
  model_input,
  duplicate_time_rows
)

gc(full = TRUE)


# ============================================================
# 7. Save configuration
# ============================================================

config_dt <- data.table(
  parameter = c(
    "VERSION_TAG",
    "ENGINE_TAG",
    "input_rds",
    "output_dir",
    "GROUPS_TO_SCREEN",
    "MAX_ACTIVE_TASKS",
    "MODEL_DEGREE",
    "SCREEN_TOL",
    "NSTART",
    "START_STAGGER_SECONDS",
    "NONBINDING_MAXIT",
    "DELETE_EMPTY",
    "PRUNING",
    "MASTER_SEED",
    "EXPECTED_PATIENTS",
    "EXPECTED_LONG_ROWS",
    "PREFERRED_MIN_CLASS_PCT",
    "PREFERRED_MIN_APPA",
    "MEMORY_MIN_AVAILABLE_GB",
    "MEMORY_MIN_AVAILABLE_FRACTION",
    "MEMORY_MAX_COMMIT_FRACTION",
    "MEMORY_POLL_SECONDS",
    "PROCESS_PRIORITY_CLASS",
    "logical_cores",
    "physical_cores",
    "cpu_classification_note",
    "R_version",
    "gbmt_version",
    "data.table_version"
  ),
  value = c(
    VERSION_TAG,
    ENGINE_TAG,
    input_rds,
    output_dir,
    paste(GROUPS_TO_SCREEN, collapse = ","),
    as.character(MAX_ACTIVE_TASKS),
    as.character(MODEL_DEGREE),
    as.character(SCREEN_TOL),
    as.character(NSTART),
    as.character(START_STAGGER_SECONDS),
    as.character(NONBINDING_MAXIT),
    "FALSE",
    "TRUE",
    as.character(MASTER_SEED),
    as.character(EXPECTED_PATIENTS),
    as.character(EXPECTED_LONG_ROWS),
    as.character(PREFERRED_MIN_CLASS_PCT),
    as.character(PREFERRED_MIN_APPA),
    as.character(MEMORY_MIN_AVAILABLE_GB),
    as.character(MEMORY_MIN_AVAILABLE_FRACTION),
    as.character(MEMORY_MAX_COMMIT_FRACTION),
    as.character(MEMORY_POLL_SECONDS),
    PROCESS_PRIORITY_CLASS,
    as.character(logical_cores),
    as.character(physical_cores),
    cpu_info$classification_note,
    R.version.string,
    as.character(packageVersion("gbmt")),
    as.character(packageVersion("data.table"))
  )
)

# Preserve the historical V251_config.csv exactly if it already exists.
# V258 execution metadata is written to a separate engine-level file.
fwrite(
  config_dt,
  file.path(
    output_dir,
    "V258_execution_config.csv"
  )
)

resume_config_dt <- data.table(
  parameter = c(
    "ENGINE_TAG",
    "COMPATIBLE_OUTPUT_VERSION",
    "CHECKPOINT_ENABLED",
    "CHECKPOINT_EVERY_ACCEPTED_ITER",
    "CHECKPOINT_SLOTS",
    "DYNAMIC_CPU_SCHEDULER",
    "PHYSICAL_BEFORE_SMT",
    "HYBRID_PRIORITY",
    "SAFE_BOUNDARY_MIGRATION",
    "GLOBAL_START_STAGGER_SECONDS",
    "MAX_ACTIVE_TASKS",
    "RAM_LAUNCH_GUARD",
    "input_md5",
    "gbmt_version"
  ),
  value = c(
    ENGINE_TAG,
    VERSION_TAG,
    "TRUE",
    "1",
    "2",
    "TRUE",
    "TRUE",
    "P>E>LP-E>SMT",
    "TRUE",
    as.character(START_STAGGER_SECONDS),
    as.character(MAX_ACTIVE_TASKS),
    "TRUE",
    as.character(input_md5),
    as.character(packageVersion("gbmt"))
  )
)

fwrite(
  resume_config_dt,
  file.path(
    output_dir,
    "V258_resume_config.csv"
  )
)

# ============================================================
# 8. Helper functions
# ============================================================

append_stage <- function(
  path,
  text
) {

  cat(
    paste0(
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      ),
      " | ",
      text,
      "\n"
    ),
    file = path,
    append = TRUE
  )

}


safe_ic_value <- function(
  fit,
  target
) {

  x <- fit$ic

  if (
    is.null(x) ||
    is.null(names(x))
  ) {
    return(NA_real_)
  }

  idx <- which(
    tolower(names(x)) ==
      tolower(target)
  )

  if (length(idx) == 0L) {
    return(NA_real_)
  }

  as.numeric(
    x[idx[1L]]
  )

}


extract_em_table <- function(
  fit
) {

  x <- fit$em

  if (is.null(x)) {
    return(
      data.table()
    )
  }

  if (is.null(dim(x))) {
    x <- matrix(
      x,
      nrow = 1L
    )
  }

  as.data.table(
    as.data.frame(
      x,
      stringsAsFactors = FALSE
    )
  )

}


extract_em_status <- function(
  fit
) {

  em_dt <- extract_em_table(
    fit
  )

  if (nrow(em_dt) == 0L) {
    return(
      list(
        converged = FALSE,
        iterations = NA_real_
      )
    )
  }

  nm <- names(em_dt)

  iter_col <- grep(
    "iter",
    nm,
    ignore.case = TRUE,
    value = TRUE
  )

  conv_col <- grep(
    "converg",
    nm,
    ignore.case = TRUE,
    value = TRUE
  )

  # gbmt 0.1.4 with nstart=1 may be converted to unnamed
  # V1/V2/V3 columns by extract_em_table(). In that case the
  # source-defined order is logLik, n.iter, converged.
  if (length(iter_col) == 0L && ncol(em_dt) >= 2L) {
    iter_col <- names(em_dt)[2L]
  }

  if (length(conv_col) == 0L && ncol(em_dt) >= 3L) {
    conv_col <- names(em_dt)[3L]
  }

  iterations <- NA_real_

  if (length(iter_col) > 0L) {

    vals <- suppressWarnings(
      as.numeric(
        em_dt[[iter_col[1L]]]
      )
    )

    vals <- vals[
      is.finite(vals)
    ]

    if (length(vals) > 0L) {
      iterations <- max(vals)
    }

  }

  converged <- FALSE

  if (length(conv_col) > 0L) {

    vals <- suppressWarnings(
      as.numeric(
        em_dt[[conv_col[1L]]]
      )
    )

    vals <- vals[
      is.finite(vals)
    ]

    if (length(vals) > 0L) {
      converged <- all(vals == 1)
    }

  }

  list(
    converged = isTRUE(converged),
    iterations = as.numeric(iterations)
  )

}


extract_class_quality <- function(
  fit,
  requested_g
) {

  post <- as.matrix(
    fit$posterior
  )

  assign_vec <- as.integer(
    fit$assign
  )

  if (nrow(post) != length(assign_vec)) {
    stop(
      "Posterior row count differs from assignment length."
    )
  }

  posterior_ncol <- ncol(post)

  stored_ng <- suppressWarnings(
    as.integer(
      fit$call$ng
    )
  )

  if (
    length(stored_ng) == 0L ||
    is.na(stored_ng)
  ) {
    stored_ng <- posterior_ncol
  }

  if (posterior_ncol < 1L) {
    stop(
      "Returned posterior matrix has zero columns."
    )
  }

  modal_post <- max.col(
    post,
    ties.method = "first"
  )

  mismatch_n <- sum(
    assign_vec != modal_post,
    na.rm = TRUE
  )

  class_dt <- rbindlist(
    lapply(
      seq_len(posterior_ncol),
      function(k) {

        idx <- which(
          assign_vec == k
        )

        n_k <- length(idx)

        appa <- if (n_k > 0L) {
          mean(
            post[
              idx,
              k
            ],
            na.rm = TRUE
          )
        } else {
          NA_real_
        }

        soft_mass <- sum(
          post[
            ,
            k
          ],
          na.rm = TRUE
        )

        data.table(
          requested_g = as.integer(requested_g),
          stored_ng = as.integer(stored_ng),
          class = as.integer(k),
          n = as.integer(n_k),
          pct = as.numeric(
            100 * n_k / length(assign_vec)
          ),
          APPA = as.numeric(appa),
          posterior_soft_mass = as.numeric(soft_mass),
          posterior_soft_mass_pct = as.numeric(
            100 * soft_mass / nrow(post)
          ),
          max_posterior_for_class = as.numeric(
            max(
              post[
                ,
                k
              ],
              na.rm = TRUE
            )
          )
        )

      }
    )
  )

  occupied_groups <- sum(
    class_dt$n > 0L
  )

  empty_classes <- class_dt[
    n == 0L,
    class
  ]

  nonempty_dt <- class_dt[
    n > 0L
  ]

  min_class_pct <- if (nrow(nonempty_dt) > 0L) {
    min(
      nonempty_dt$pct,
      na.rm = TRUE
    )
  } else {
    NA_real_
  }

  min_appa <- if (
    nrow(nonempty_dt) > 0L &&
    any(is.finite(nonempty_dt$APPA))
  ) {
    min(
      nonempty_dt$APPA[
        is.finite(nonempty_dt$APPA)
      ]
    )
  } else {
    NA_real_
  }

  max_post <- apply(
    post,
    1L,
    max,
    na.rm = TRUE
  )

  structure_exact <- (
    stored_ng == requested_g
  ) &&
    (
      posterior_ncol == requested_g
    )

  all_occupied <- (
    occupied_groups == requested_g
  )

  full_g <- structure_exact &&
    all_occupied

  summary_dt <- data.table(
    requested_g = as.integer(requested_g),
    stored_ng = as.integer(stored_ng),
    posterior_ncol = as.integer(posterior_ncol),
    occupied_groups = as.integer(occupied_groups),
    structure_exact = as.integer(structure_exact),
    all_classes_occupied = as.integer(all_occupied),
    full_g = as.integer(full_g),
    empty_classes = if (
      length(empty_classes) > 0L
    ) {
      paste(
        empty_classes,
        collapse = ","
      )
    } else {
      ""
    },
    assignment_posterior_mismatch_n = as.integer(mismatch_n),
    min_class_pct_nonempty = as.numeric(min_class_pct),
    min_APPA_nonempty = as.numeric(min_appa),
    mean_max_posterior = as.numeric(
      mean(
        max_post,
        na.rm = TRUE
      )
    ),
    median_max_posterior = as.numeric(
      median(
        max_post,
        na.rm = TRUE
      )
    ),
    pct_max_posterior_ge_0_70 = as.numeric(
      100 *
        mean(
          max_post >= 0.70,
          na.rm = TRUE
        )
    ),
    pct_max_posterior_ge_0_90 = as.numeric(
      100 *
        mean(
          max_post >= 0.90,
          na.rm = TRUE
        )
    )
  )

  list(
    class_table = class_dt,
    summary = summary_dt
  )

}



# ============================================================
# 8A. V252 checkpoint/resume engine
# ============================================================
#
# Scientific model is intentionally unchanged from V251:
#   - same gbmt 0.1.4 EM update equations
#   - same random initialization for nstart = 1
#   - same d = 2, scaling = 0, pruning = TRUE
#   - same tol = 10 and delete.empty = FALSE
#
# The only change is execution resilience. The accepted EM state is
# written after initialization and after every accepted non-final EM
# iteration. Two alternating checkpoint files are used so that a
# power loss during one checkpoint write still leaves the other slot.
#
# Checkpoints store only the minimal state needed to continue the EM:
#   p0, beta, Sigma, postd, Z, ll_old, iteration, count2.
# Large lm/regression objects are deliberately NOT checkpointed.
# With pruning=TRUE they are rebuilt by gbmt's own pruneFun() during
# finalization, so this keeps checkpoint files small.
# ============================================================

CHECKPOINT_STATE_VERSION <- 1L
CHECKPOINT_EVERY_ACCEPTED_ITER <- 1L

make_checkpoint_signature <- function(
  g,
  seed_used,
  input_md5,
  z_markers,
  model_degree,
  screen_tol,
  nonbinding_maxit
) {

  paste(
    "GBMTM_V252",
    paste0("G=", as.integer(g)),
    paste0("seed=", as.integer(seed_used)),
    paste0("input_md5=", as.character(input_md5)),
    paste0("x=", paste(z_markers, collapse = ",")),
    paste0("d=", as.integer(model_degree)),
    paste0("tol=", format(as.numeric(screen_tol), scientific = FALSE, trim = TRUE)),
    paste0("maxit=", as.integer(nonbinding_maxit)),
    "scaling=0",
    "pruning=TRUE",
    "delete.empty=FALSE",
    "nstart=1",
    sep = "|"
  )

}

append_progress_line <- function(
  path,
  text
) {

  cat(
    paste0(
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      ),
      " | ",
      text,
      "\n"
    ),
    file = path,
    append = TRUE
  )

}

validate_checkpoint_state <- function(
  state,
  expected_signature,
  expected_g
) {

  required_names <- c(
    "state_version",
    "signature",
    "requested_g",
    "iteration",
    "count2",
    "ll_old",
    "p0",
    "beta",
    "Sigma",
    "postd",
    "Z",
    "saved_at"
  )

  if (!is.list(state)) {
    return(
      list(
        valid = FALSE,
        reason = "checkpoint is not a list"
      )
    )
  }

  missing_names <- setdiff(
    required_names,
    names(state)
  )

  if (length(missing_names) > 0L) {
    return(
      list(
        valid = FALSE,
        reason = paste0(
          "missing checkpoint field(s): ",
          paste(missing_names, collapse = ", ")
        )
      )
    )
  }

  if (!identical(
    as.integer(state$state_version),
    as.integer(CHECKPOINT_STATE_VERSION)
  )) {
    return(
      list(
        valid = FALSE,
        reason = paste0(
          "checkpoint state version mismatch: ",
          state$state_version
        )
      )
    )
  }

  if (!identical(
    as.character(state$signature),
    as.character(expected_signature)
  )) {
    return(
      list(
        valid = FALSE,
        reason = "checkpoint signature does not match current data/settings"
      )
    )
  }

  if (!identical(
    as.integer(state$requested_g),
    as.integer(expected_g)
  )) {
    return(
      list(
        valid = FALSE,
        reason = "checkpoint G does not match requested G"
      )
    )
  }

  if (
    length(state$p0) != expected_g ||
    length(state$beta) != expected_g ||
    length(state$Sigma) != expected_g ||
    ncol(as.matrix(state$postd)) != expected_g
  ) {
    return(
      list(
        valid = FALSE,
        reason = "checkpoint group dimension does not match requested G"
      )
    )
  }

  if (
    !is.numeric(state$iteration) ||
    length(state$iteration) != 1L ||
    !is.finite(state$iteration) ||
    state$iteration < 0
  ) {
    return(
      list(
        valid = FALSE,
        reason = "checkpoint iteration is invalid"
      )
    )
  }

  if (
    !is.numeric(state$ll_old) ||
    length(state$ll_old) != 1L ||
    !is.finite(state$ll_old)
  ) {
    return(
      list(
        valid = FALSE,
        reason = "checkpoint ll_old is invalid"
      )
    )
  }

  list(
    valid = TRUE,
    reason = "OK"
  )

}

read_checkpoint_slot <- function(
  path,
  expected_signature,
  expected_g
) {

  if (!file.exists(path)) {
    return(
      list(
        exists = FALSE,
        valid = FALSE,
        path = path,
        state = NULL,
        reason = "missing"
      )
    )
  }

  state <- tryCatch(
    readRDS(path),
    error = function(e) e
  )

  if (inherits(state, "error")) {
    return(
      list(
        exists = TRUE,
        valid = FALSE,
        path = path,
        state = NULL,
        reason = paste0(
          "readRDS failed: ",
          conditionMessage(state)
        )
      )
    )
  }

  verdict <- validate_checkpoint_state(
    state = state,
    expected_signature = expected_signature,
    expected_g = expected_g
  )

  list(
    exists = TRUE,
    valid = isTRUE(verdict$valid),
    path = path,
    state = if (isTRUE(verdict$valid)) state else NULL,
    reason = verdict$reason
  )

}

load_latest_checkpoint <- function(
  checkpoint_a,
  checkpoint_b,
  expected_signature,
  expected_g,
  stage_path
) {

  a <- read_checkpoint_slot(
    path = checkpoint_a,
    expected_signature = expected_signature,
    expected_g = expected_g
  )

  b <- read_checkpoint_slot(
    path = checkpoint_b,
    expected_signature = expected_signature,
    expected_g = expected_g
  )

  any_exists <- isTRUE(a$exists) || isTRUE(b$exists)

  valid_slots <- list()

  if (isTRUE(a$valid)) {
    valid_slots[[length(valid_slots) + 1L]] <- a
  }

  if (isTRUE(b$valid)) {
    valid_slots[[length(valid_slots) + 1L]] <- b
  }

  if (length(valid_slots) == 0L) {

    if (!any_exists) {
      return(NULL)
    }

    stop(
      paste0(
        "Checkpoint file(s) exist for G",
        expected_g,
        " but no valid checkpoint can be loaded.\n",
        "Slot A: ",
        a$reason,
        "\nSlot B: ",
        b$reason,
        "\nThe script will NOT silently start from scratch."
      )
    )

  }

  if (length(valid_slots) == 1L) {

    chosen <- valid_slots[[1L]]

  } else {

    iter_a <- as.integer(valid_slots[[1L]]$state$iteration)
    iter_b <- as.integer(valid_slots[[2L]]$state$iteration)

    if (iter_a > iter_b) {
      chosen <- valid_slots[[1L]]
    } else if (iter_b > iter_a) {
      chosen <- valid_slots[[2L]]
    } else {

      mt_a <- file.info(valid_slots[[1L]]$path)$mtime
      mt_b <- file.info(valid_slots[[2L]]$path)$mtime

      if (is.na(mt_b) || (!is.na(mt_a) && mt_a >= mt_b)) {
        chosen <- valid_slots[[1L]]
      } else {
        chosen <- valid_slots[[2L]]
      }

    }

  }

  append_stage(
    stage_path,
    paste0(
      "valid checkpoint loaded: ",
      basename(chosen$path),
      "; iteration=",
      chosen$state$iteration,
      "; ll=",
      format(chosen$state$ll_old, digits = 15)
    )
  )

  chosen$state

}

save_checkpoint_dual_slot <- function(
  state,
  checkpoint_a,
  checkpoint_b
) {

  iteration <- as.integer(state$iteration)

  target <- if (
    iteration %% 2L == 0L
  ) {
    checkpoint_a
  } else {
    checkpoint_b
  }

  tmp <- paste0(
    target,
    ".tmp_",
    Sys.getpid()
  )

  if (file.exists(tmp)) {
    file.remove(tmp)
  }

  saveRDS(
    state,
    tmp,
    compress = FALSE
  )

  if (
    !file.exists(tmp) ||
    is.na(file.info(tmp)$size) ||
    file.info(tmp)$size <= 0
  ) {
    stop(
      paste0(
        "Checkpoint temporary file was not written correctly: ",
        tmp
      )
    )
  }

  test_state <- tryCatch(
    readRDS(tmp),
    error = function(e) e
  )

  if (inherits(test_state, "error")) {
    stop(
      paste0(
        "Checkpoint verification read failed: ",
        conditionMessage(test_state)
      )
    )
  }

  if (
    !identical(
      as.integer(test_state$iteration),
      iteration
    ) ||
    !identical(
      as.character(test_state$signature),
      as.character(state$signature)
    )
  ) {
    stop(
      "Checkpoint verification failed after writing temporary file."
    )
  }

  rm(test_state)

  if (file.exists(target)) {
    file.remove(target)
  }

  renamed <- file.rename(
    tmp,
    target
  )

  if (!renamed) {

    copied <- file.copy(
      tmp,
      target,
      overwrite = TRUE
    )

    if (!copied) {
      stop(
        paste0(
          "Could not finalize checkpoint file: ",
          target
        )
      )
    }

    file.remove(tmp)

  }

  invisible(target)

}

make_checkpoint_state <- function(
  signature,
  requested_g,
  iteration,
  count2,
  ll_old,
  p0,
  beta,
  Sigma,
  postd,
  Z
) {

  list(
    state_version = as.integer(CHECKPOINT_STATE_VERSION),
    signature = as.character(signature),
    requested_g = as.integer(requested_g),
    iteration = as.integer(iteration),
    count2 = as.integer(count2),
    ll_old = as.numeric(ll_old),
    p0 = p0,
    beta = beta,
    Sigma = Sigma,
    postd = postd,
    Z = Z,
    saved_at = format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  )

}

# ------------------------------------------------------------
# Checkpoint-enabled implementation of gbmt 0.1.4 gbmtFit().
# The EM mathematics and stopping rules follow the package source.
# ------------------------------------------------------------

gbmt_fit_checkpoint_v252 <- function(
  x.names,
  unit,
  time,
  ng,
  d,
  data,
  pruning,
  delete.empty,
  tol,
  maxit,
  init,
  quiet,
  id_restart,
  nch0,
  checkpoint_a,
  checkpoint_b,
  checkpoint_signature,
  checkpoint_every,
  em_progress_path,
  stage_path,
  scheduler_context = NULL
) {

  if (!isTRUE(pruning)) {
    stop(
      "V252 checkpoint engine is intentionally restricted to pruning=TRUE."
    )
  }

  if (!identical(delete.empty, FALSE)) {
    stop(
      "V252 checkpoint engine is intentionally restricted to delete.empty=FALSE."
    )
  }

  wlsFit <- getFromNamespace(
    "wlsFit",
    "gbmt"
  )

  imputFun <- getFromNamespace(
    "imputFun",
    "gbmt"
  )

  logdmN <- getFromNamespace(
    "logdmN",
    "gbmt"
  )

  pldBound <- getFromNamespace(
    "pldBound",
    "gbmt"
  )

  pruneFun <- getFromNamespace(
    "pruneFun",
    "gbmt"
  )

  nparCalc <- getFromNamespace(
    "nparCalc",
    "gbmt"
  )

  icCalc <- getFromNamespace(
    "icCalc",
    "gbmt"
  )

  appaCalc <- getFromNamespace(
    "appaCalc",
    "gbmt"
  )

  data[, unit] <- factor(
    data[, unit]
  )

  tnam <- levels(
    factor(
      data[, time]
    )
  )

  tval <- as.numeric(
    data[, time] -
      min(data[, time]) +
      1
  )

  ismiss <- sum(
    is.na(
      data[, x.names]
    )
  ) > 0

  cou <- levels(
    factor(
      data[, unit]
    )
  )

  n <- length(cou)
  np <- length(x.names)
  nt <- length(tnam)

  Xmat <- matrix(
    nrow = nt,
    ncol = d + 1L
  )

  for (i in 0:d) {
    Xmat[, i + 1L] <- sort(
      unique(tval)
    )^i
  }

  rownames(Xmat) <- tnam
  colnames(Xmat) <- paste0(
    "I(t^",
    0:d,
    ")"
  )

  postp <- matrix(
    nrow = n,
    ncol = ng
  )

  rownames(postp) <- cou

  dataI <- data

  if (ismiss) {

    for (j in seq_along(x.names)) {

      ijna <- which(
        is.na(
          data[, x.names[j]]
        )
      )

      if (length(ijna) > 0L) {
        dataI[ijna, x.names[j]] <- mean(
          data[, x.names[j]],
          na.rm = TRUE
        )
      }

    }

    dataI_list <- vector(
      "list",
      length = ng
    )

  } else {

    dataI_list <- NULL

  }

  resume_state <- load_latest_checkpoint(
    checkpoint_a = checkpoint_a,
    checkpoint_b = checkpoint_b,
    expected_signature = checkpoint_signature,
    expected_g = ng,
    stage_path = stage_path
  )

  resumed <- !is.null(resume_state)
  resumed_from_iteration <- 0L

  if (!resumed) {

    if (is.null(init)) {
      Z <- Z0 <- sample(
        rep(
          1:ng,
          length = n
        )
      )
    } else {
      Z <- Z0 <- init
    }

    names(Z) <- names(Z0) <- cou

    postd <- matrix(
      nrow = n,
      ncol = ng
    )

    rownames(postd) <- cou

    p0 <- prop.table(
      table(
        factor(
          Z,
          levels = 1:ng
        )
      )
    )

    beta <- Sigma <- reg <- vector(
      "list",
      length = ng
    )

    append_progress_line(
      em_progress_path,
      paste0(
        "G",
        ng,
        " fresh initialization started"
      )
    )

    for (i in 1:ng) {

      iwei <- rep(
        0,
        nrow(data)
      )

      iwei[
        which(
          data[, unit] %in%
            cou[
              which(
                Z == i
              )
            ]
        )
      ] <- 1

      imod <- wlsFit(
        x.names = x.names,
        data = dataI,
        t = tval,
        weights = iwei,
        d = d
      )

      reg[[i]] <- imod$fit
      beta[[i]] <- imod$beta
      Sigma[[i]] <- imod$S

      for (w in 1:n) {

        auxw <- which(
          data[, unit] == cou[w]
        )

        wdat <- as.matrix(
          dataI[
            auxw,
            x.names,
            drop = FALSE
          ]
        )

        if (np > 1L) {

          ijD <- sum(
            logdmN(
              wdat,
              Xmat %*% beta[[i]],
              Sigma[[i]]
            )
          )

        } else {

          ijD <- sum(
            dnorm(
              dataI[auxw, x.names],
              Xmat %*% beta[[i]],
              sqrt(Sigma[[i]]),
              log = TRUE
            )
          )

        }

        postd[w, i] <- pldBound(
          ijD + log(p0[i])
        )

      }

    }

    ll_old <- sum(
      log(
        apply(
          exp(postd),
          1L,
          sum
        )
      )
    )

    count <- 0L
    count2 <- 0L

    init_state <- make_checkpoint_state(
      signature = checkpoint_signature,
      requested_g = ng,
      iteration = count,
      count2 = count2,
      ll_old = ll_old,
      p0 = p0,
      beta = beta,
      Sigma = Sigma,
      postd = postd,
      Z = Z
    )

    save_checkpoint_dual_slot(
      state = init_state,
      checkpoint_a = checkpoint_a,
      checkpoint_b = checkpoint_b
    )

    append_progress_line(
      em_progress_path,
      paste0(
        "Restart 1 - EM initialization. Log likelihood: ",
        round(ll_old, 4),
        " | checkpoint saved"
      )
    )

    if (!is.null(scheduler_context) && isTRUE(scheduler_context$enabled)) {
      scheduler_iteration_boundary(
        g = scheduler_context$g,
        pid = scheduler_context$pid,
        iteration = 0L,
        delta_ll = NA_real_,
        iter_seconds = 0,
        state_path = scheduler_context$state_path,
        lock_dir = scheduler_context$lock_dir,
        event_log = scheduler_context$event_log,
        helper_path = scheduler_context$helper_path,
        priority_class = scheduler_context$priority_class
      )
    }

  } else {

    count <- as.integer(
      resume_state$iteration
    )

    count2 <- as.integer(
      resume_state$count2
    )

    ll_old <- as.numeric(
      resume_state$ll_old
    )

    p0 <- resume_state$p0
    beta <- resume_state$beta
    Sigma <- resume_state$Sigma
    postd <- as.matrix(
      resume_state$postd
    )
    Z <- resume_state$Z

    names(Z) <- cou
    rownames(postd) <- cou

    reg <- vector(
      "list",
      length = ng
    )

    resumed_from_iteration <- count

    append_progress_line(
      em_progress_path,
      paste0(
        "Restart 1 - RESUMED from checkpoint iteration ",
        count,
        ". Log likelihood: ",
        round(ll_old, 4)
      )
    )

    if (!is.null(scheduler_context) && isTRUE(scheduler_context$enabled)) {
      scheduler_iteration_boundary(
        g = scheduler_context$g,
        pid = scheduler_context$pid,
        iteration = count,
        delta_ll = NA_real_,
        iter_seconds = 0,
        state_path = scheduler_context$state_path,
        lock_dir = scheduler_context$lock_dir,
        event_log = scheduler_context$event_log,
        helper_path = scheduler_context$helper_path,
        priority_class = scheduler_context$priority_class
      )
    }

  }

  fine <- 0L

  p0_old <- p0
  beta_old <- beta
  Sigma_old <- Sigma
  reg_old <- reg

  while (fine == 0L) {

    iteration_start_time <- Sys.time()
    count <- count + 1L

    # E-step: compute posterior from accepted postd.
    for (i in 1:nrow(postd)) {
      postp[i, ] <- exp(postd[i, ]) /
        sum(
          exp(postd[i, ])
        )
    }

    Z_old <- Z

    Z <- apply(
      postp,
      1L,
      function(x) {
        which.max(x)
      }
    )

    # M-step: prior, beta, Sigma.
    p0 <- colMeans(postp)

    for (i in 1:ng) {

      iwei <- numeric(0)

      for (j in 1:length(cou)) {
        iwei[
          which(
            data[, unit] == cou[j]
          )
        ] <- postp[j, i]
      }

      if (ismiss) {

        idat <- dataI_list[[i]] <- imputFun(
          xdat = data[, x.names, drop = FALSE],
          tval = tval,
          beta = beta[[i]],
          S = Sigma[[i]],
          d = d
        )

      } else {

        idat <- data[, x.names, drop = FALSE]

      }

      imod <- wlsFit(
        x.names = x.names,
        data = idat,
        t = tval,
        weights = iwei,
        d = d
      )

      reg[[i]] <- imod$fit
      beta[[i]] <- imod$beta
      Sigma[[i]] <- imod$S

      for (w in 1:n) {

        auxw <- which(
          data[, unit] == cou[w]
        )

        wdat <- as.matrix(
          idat[
            auxw,
            x.names,
            drop = FALSE
          ]
        )

        if (np > 1L) {

          ijD <- sum(
            logdmN(
              wdat,
              Xmat %*% beta[[i]],
              Sigma[[i]]
            )
          )

        } else {

          ijD <- sum(
            dnorm(
              wdat,
              Xmat %*% beta[[i]],
              sqrt(Sigma[[i]]),
              log = TRUE
            )
          )

        }

        postd[w, i] <- pldBound(
          ijD + log(p0[i])
        )

      }

    }

    ll_new <- sum(
      log(
        apply(
          exp(postd),
          1L,
          sum
        )
      )
    )

    delta_ll <- ll_new - ll_old

    # Exact gbmt 0.1.4 stopping logic.
    if (ll_new >= ll_old) {
      if (ll_new - ll_old < tol) {
        fine <- 1L
      }
    } else {

      p0 <- p0_old
      beta <- beta_old
      Sigma <- Sigma_old
      reg <- reg_old
      ll_new <- ll_old
      count <- count - 1L
      fine <- 1L

    }

    append_progress_line(
      em_progress_path,
      paste0(
        "Restart 1 - EM iteration ",
        count,
        ". Log likelihood: ",
        round(ll_new, 4),
        "; delta: ",
        round(ll_new - ll_old, 6),
        if (fine == 1L) " | stopping condition reached" else ""
      )
    )

    if (count >= maxit) {
      fine <- 1L
    }

    if (sum(Z_old - Z) == 0) {

      count2 <- count2 + 1L

      if (
        ll_new - ll_old <
          tol * sqrt(count2)
      ) {
        fine <- 1L
      }

    } else {

      count2 <- 0L

    }

    iter_seconds_current <- as.numeric(
      difftime(
        Sys.time(),
        iteration_start_time,
        units = "secs"
      )
    )

    if (fine == 0L) {

      ll_old <- ll_new
      p0_old <- p0
      beta_old <- beta
      Sigma_old <- Sigma
      reg_old <- reg

      if (
        count %%
          as.integer(checkpoint_every) ==
          0L
      ) {

        checkpoint_state <- make_checkpoint_state(
          signature = checkpoint_signature,
          requested_g = ng,
          iteration = count,
          count2 = count2,
          ll_old = ll_old,
          p0 = p0,
          beta = beta,
          Sigma = Sigma,
          postd = postd,
          Z = Z
        )

        saved_checkpoint <- save_checkpoint_dual_slot(
          state = checkpoint_state,
          checkpoint_a = checkpoint_a,
          checkpoint_b = checkpoint_b
        )

        append_progress_line(
          em_progress_path,
          paste0(
            "Checkpoint saved after accepted iteration ",
            count,
            ": ",
            basename(saved_checkpoint)
          )
        )

        if (!is.null(scheduler_context) && isTRUE(scheduler_context$enabled)) {
          scheduler_iteration_boundary(
            g = scheduler_context$g,
            pid = scheduler_context$pid,
            iteration = count,
            delta_ll = delta_ll,
            iter_seconds = iter_seconds_current,
            state_path = scheduler_context$state_path,
            lock_dir = scheduler_context$lock_dir,
            event_log = scheduler_context$event_log,
            helper_path = scheduler_context$helper_path,
            priority_class = scheduler_context$priority_class
          )
        }

      }

    }

  }

  # Final posterior and assignments, matching gbmt 0.1.4.
  for (i in 1:nrow(postd)) {
    postp[i, ] <- exp(postd[i, ]) /
      sum(
        exp(postd[i, ])
      )
  }

  Z <- apply(
    postp,
    1L,
    function(x) {
      which.max(x)
    }
  )

  assL <- vector(
    "list",
    length = ng
  )

  for (i in 1:ng) {
    assL[[i]] <- sort(
      names(
        which(
          Z == i
        )
      )
    )
  }

  # Imputed data, matching gbmt 0.1.4.
  if (ismiss) {

    auxdat <- dataI_list[[1]][
      ,
      x.names,
      drop = FALSE
    ] * p0[1]

    if (length(dataI_list) > 1L) {

      for (i in 2:length(dataI_list)) {
        auxdat <- auxdat +
          dataI_list[[i]][
            ,
            x.names,
            drop = FALSE
          ] * p0[i]
      }

    }

    newdat <- data.frame(
      data[, c(unit, time)],
      auxdat
    )

  } else {

    newdat <- data[
      ,
      c(
        unit,
        time,
        x.names
      )
    ]

  }

  # delete.empty is frozen FALSE, but retain source-compatible branch.
  if (delete.empty) {

    auxchk <- which(
      sapply(
        assL,
        length
      ) == 0
    )

    if (length(auxchk) > 0L) {

      gOK <- setdiff(
        1:ng,
        auxchk
      )

      beta <- beta[gOK]
      Sigma <- Sigma[gOK]
      postp <- postp[, gOK, drop = FALSE]
      reg <- reg[gOK]
      p0 <- p0[gOK] /
        sum(p0[gOK])
      assL <- assL[gOK]
      postd <- postd[, gOK, drop = FALSE]

      ll_new <- sum(
        log(
          apply(
            exp(postd),
            1L,
            sum
          )
        )
      )

      warning(
        ng - length(gOK),
        " empty groups have been deleted"
      )

      ng <- length(gOK)

    }

  }

  # V251 uses pruning=TRUE. This rebuilds regression objects and is
  # why they do not need to be stored in the checkpoint.
  if (pruning) {

    mpruned <- pruneFun(
      beta = beta,
      posterior = postp,
      data = newdat
    )

    reg <- mpruned$reg
    beta <- mpruned$beta
    Sigma <- mpruned$S
    ll_new <- mpruned$logLik

  }

  npar <- nparCalc(beta)

  ic <- icCalc(
    ll_new,
    npar = npar,
    ss = nrow(newdat)
  )

  pred <- vector(
    "list",
    length = ng
  )

  names(pred) <- 1:ng

  for (i in 1:ng) {

    pred[[i]] <- data.frame(
      Xmat %*% beta[[i]]
    )

    rownames(pred[[i]]) <- tnam
    colnames(pred[[i]]) <- x.names

  }

  names(p0) <-
    names(assL) <-
    colnames(postp) <-
    names(beta) <-
    names(Sigma) <-
    names(reg) <-
    1:ng

  datOK <- data[
    ,
    c(
      unit,
      time,
      x.names
    )
  ]

  rownames(datOK) <- NULL

  res <- list(
    call = list(
      x.names = x.names,
      unit = unit,
      time = time,
      ng = ng,
      d = d
    ),
    prior = p0,
    beta = beta,
    Sigma = Sigma,
    posterior = round(
      postp,
      4
    ),
    Xmat = Xmat,
    fitted = pred,
    reg = reg,
    assign = Z,
    assign.list = assL,
    logLik = ll_new,
    npar = npar,
    ic = ic,
    appa = appaCalc(postp),
    data.orig = NULL,
    data.norm = newdat,
    data.imputed = NULL,
    em = c(
      n.iter = count,
      converged = ifelse(
        count >= maxit,
        0,
        1
      )
    ),
    nch = 0L
  )

  class(res) <- "gbmt"

  list(
    fit = res,
    resume_meta = list(
      resumed = as.integer(resumed),
      resumed_from_iteration = as.integer(
        resumed_from_iteration
      ),
      final_em_iteration = as.integer(count)
    )
  )

}

# ------------------------------------------------------------
# One-start wrapper reproducing the relevant gbmt 0.1.4 front end.
# It is intentionally specialized to the frozen V251 model settings.
# ------------------------------------------------------------

gbmt_checkpoint_one_start_v252 <- function(
  x.names,
  unit,
  time,
  ng,
  d,
  data,
  scaling,
  pruning,
  delete.empty,
  nstart,
  tol,
  maxit,
  quiet,
  checkpoint_a,
  checkpoint_b,
  checkpoint_signature,
  checkpoint_every,
  em_progress_path,
  stage_path,
  scheduler_context = NULL
) {

  if (!identical(as.integer(nstart), 1L)) {
    stop(
      "V252 checkpoint engine is restricted to nstart=1 to match V251."
    )
  }

  if (!identical(as.numeric(scaling), 0)) {
    stop(
      "V252 checkpoint engine is restricted to scaling=0 to match V251."
    )
  }

  scaleFun <- getFromNamespace(
    "scaleFun",
    "gbmt"
  )

  scaleFun_inv <- getFromNamespace(
    "scaleFun_inv",
    "gbmt"
  )

  if (!identical(class(data), "data.frame")) {
    stop(
      "Argument 'data' must be a data.frame."
    )
  }

  data[, unit] <- factor(
    data[, unit]
  )

  cou <- levels(
    data[, unit]
  )

  # Same within-unit time sort as gbmt 0.1.4.
  for (i in 1:length(cou)) {

    ind <- which(
      data[, unit] == cou[i]
    )

    idat <- data[ind, ]

    if (
      sum(
        duplicated(
          idat[, time]
        )
      ) > 0
    ) {
      stop(
        paste0(
          "Unit '",
          cou[i],
          "' has duplicated time points"
        )
      )
    }

    data[ind, ] <- idat[
      order(
        idat[, time]
      ),
    ]

  }

  # Same scaling front end. For this study scaling=0, hence identity.
  dataOK <- data

  for (i in 1:length(cou)) {

    ind <- which(
      data[, unit] == cou[i]
    )

    for (j in 1:length(x.names)) {

      ijdat <- data[
        ind,
        x.names[j]
      ]

      ijpar <- c(
        mean(
          ijdat,
          na.rm = TRUE
        ),
        sd(
          ijdat,
          na.rm = TRUE
        )
      )

      if (is.na(ijpar[1])) {
        ijpar[1] <- mean(
          data[, x.names[j]],
          na.rm = TRUE
        )
      }

      if (
        is.na(ijpar[2]) ||
        ijpar[2] == 0
      ) {
        ijpar[2] <- sd(
          data[, x.names[j]],
          na.rm = TRUE
        )
      }

      dataOK[
        ind,
        x.names[j]
      ] <- scaleFun(
        ijdat,
        scaling = scaling,
        par = ijpar
      )

    }

  }

  fit_result <- gbmt_fit_checkpoint_v252(
    x.names = x.names,
    unit = unit,
    time = time,
    ng = ng,
    d = d,
    data = dataOK,
    pruning = pruning,
    delete.empty = delete.empty,
    tol = tol,
    maxit = maxit,
    init = NULL,
    quiet = quiet,
    id_restart = 1L,
    nch0 = 0L,
    checkpoint_a = checkpoint_a,
    checkpoint_b = checkpoint_b,
    checkpoint_signature = checkpoint_signature,
    checkpoint_every = checkpoint_every,
    em_progress_path = em_progress_path,
    stage_path = stage_path,
    scheduler_context = scheduler_context
  )

  mOK <- fit_result$fit

  res <- c(
    logLik = mOK$logLik,
    mOK$em
  )

  rownames(res) <- NULL

  mOK$data.orig <- data[
    ,
    c(
      unit,
      time,
      x.names
    )
  ]

  if (
    sum(
      is.na(
        data[, x.names]
      )
    ) > 0
  ) {

    dataI <- mOK$data.norm

    for (i in 1:length(cou)) {

      ind <- which(
        data[, unit] == cou[i]
      )

      for (j in 1:length(x.names)) {

        ijdat <- data[
          ind,
          x.names[j]
        ]

        ijpar <- c(
          mean(
            ijdat,
            na.rm = TRUE
          ),
          sd(
            ijdat,
            na.rm = TRUE
          )
        )

        if (is.na(ijpar[1])) {
          ijpar[1] <- mean(
            data[, x.names[j]],
            na.rm = TRUE
          )
        }

        if (
          is.na(ijpar[2]) ||
          ijpar[2] == 0
        ) {
          ijpar[2] <- sd(
            data[, x.names[j]],
            na.rm = TRUE
          )
        }

        dataI[
          ind,
          x.names[j]
        ] <- scaleFun_inv(
          mOK$data.norm[
            ind,
            x.names[j]
          ],
          scaling = scaling,
          par = ijpar
        )

      }

    }

    mOK$data.imputed <- dataI

  } else {

    mOK$data.imputed <- data[
      ,
      c(
        unit,
        time,
        x.names
      )
    ]

  }

  mOK$em <- res
  mOK$call$scaling <- scaling
  mOK$call$pruning <- pruning
  mOK$call$nstart <- nstart
  mOK$call$tol <- tol
  mOK$call$maxit <- maxit
  mOK$nch <- NULL

  list(
    fit = mOK,
    resume_meta = fit_result$resume_meta
  )

}



task_is_complete <- function(
  g,
  output_dir
) {

  summary_path <- file.path(
    output_dir,
    paste0(
      "V251_task_summary_G",
      g,
      ".csv"
    )
  )

  fit_path <- file.path(
    output_dir,
    paste0(
      "V251_screen_fit_G",
      g,
      ".rds"
    )
  )

  if (
    !file.exists(summary_path) ||
    !file.exists(fit_path)
  ) {
    return(FALSE)
  }

  d <- tryCatch(
    fread(summary_path),
    error = function(e) NULL
  )

  if (
    is.null(d) ||
    nrow(d) < 1L ||
    !"status" %in% names(d)
  ) {
    return(FALSE)
  }

  # Any completed fit is intentionally treated as finished.
  # V251 does not silently retry empty-class solutions.
  completed_statuses <- c(
    "CONVERGED_VALID_FULL_G",
    "CONVERGED_EMPTY_CLASS",
    "CONVERGED_STRUCTURE_MISMATCH",
    "NONCONVERGED_FULL_G",
    "NONCONVERGED_EMPTY_CLASS",
    "NONCONVERGED_STRUCTURE_MISMATCH"
  )

  as.character(
    d$status[1L]
  ) %in%
    completed_statuses

}



v251_summary_col_classes <- function() {

  list(
    character = c(
      "version",
      "status",
      "empty_classes",
      "start_time",
      "end_time",
      "model_file",
      "error_message"
    ),
    integer = c(
      "requested_g",
      "package_converged",
      "stored_ng",
      "posterior_ncol",
      "occupied_groups",
      "structure_exact",
      "all_classes_occupied",
      "full_g",
      "nstart",
      "delete_empty",
      "technical_maxit",
      "automatic_retry",
      "assignment_posterior_mismatch_n",
      "preferred_class_size_ok",
      "preferred_APPA_ok",
      "IC_comparison_eligible",
      "seed",
      "n_units",
      "n_rows",
      "pid",
      "resumed_skip"
    ),
    numeric = c(
      "screening_tol",
      "em_iterations",
      "final_logLik",
      "AIC",
      "BIC",
      "CAIC",
      "ssBIC",
      "HQIC",
      "min_class_pct_nonempty",
      "min_APPA_nonempty",
      "mean_max_posterior",
      "median_max_posterior",
      "pct_max_posterior_ge_0_70",
      "pct_max_posterior_ge_0_90",
      "fit_seconds",
      "total_task_seconds"
    )
  )
}


read_v251_summary_exact <- function(path) {

  data.table::fread(
    path,
    colClasses = v251_summary_col_classes(),
    showProgress = FALSE
  )
}


# V258 regression guard for the exact bug observed in V257:
# data.table j can resolve a column named `status` instead of a function-local
# variable also named `status`.  This test is intentionally independent of
# the user's files and must pass before any legacy V251 metadata is rewritten.
run_v258_nse_regression_guard <- function() {

  old_status_value <- "NONCONVERGED_FULL_G"
  computed_status_value <- "CONVERGED_VALID_FULL_G"
  computed_iterations_value <- 103

  probe <- data.table(
    status = as.character(computed_status_value),
    package_converged = 1L,
    em_iterations = as.numeric(computed_iterations_value)
  )

  if (
    nrow(probe) != 1L ||
    !identical(
      as.character(probe$status[1L]),
      computed_status_value
    ) ||
    identical(
      as.character(probe$status[1L]),
      old_status_value
    ) ||
    !identical(
      as.integer(probe$package_converged[1L]),
      1L
    ) ||
    !isTRUE(all.equal(
      as.numeric(probe$em_iterations[1L]),
      as.numeric(computed_iterations_value),
      tolerance = 0
    ))
  ) {
    stop(
      "V258 internal NSE regression guard failed. No V251 metadata was modified."
    )
  }

  invisible(TRUE)
}

run_v258_nse_regression_guard()


refresh_existing_v251_metadata <- function(
  g,
  output_dir,
  preferred_min_class_pct,
  preferred_min_appa
) {

  summary_path <- file.path(
    output_dir,
    paste0(
      "V251_task_summary_G",
      g,
      ".csv"
    )
  )

  fit_path <- file.path(
    output_dir,
    paste0(
      "V251_screen_fit_G",
      g,
      ".rds"
    )
  )

  class_path <- file.path(
    output_dir,
    paste0(
      "V251_class_quality_G",
      g,
      ".csv"
    )
  )

  em_runs_path <- file.path(
    output_dir,
    paste0(
      "V251_EM_runs_G",
      g,
      ".csv"
    )
  )

  done_path <- file.path(
    output_dir,
    paste0(
      "V251_DONE_G",
      g,
      ".txt"
    )
  )

  audit_path <- file.path(
    output_dir,
    paste0(
      "V258_legacy_refresh_G",
      g,
      ".csv"
    )
  )

  summary_backup_path <- file.path(
    output_dir,
    paste0(
      "V258_legacy_summary_G",
      g,
      ".csv"
    )
  )

  done_backup_path <- file.path(
    output_dir,
    paste0(
      "V258_legacy_done_G",
      g,
      ".txt"
    )
  )

  if (
    !file.exists(summary_path) ||
    !file.exists(fit_path)
  ) {
    return(
      invisible(FALSE)
    )
  }

  old_summary <- read_v251_summary_exact(
    summary_path
  )

  if (nrow(old_summary) != 1L) {
    stop(
      paste0(
        "Existing V251 summary must contain exactly one row:\n",
        normalizePath(
          summary_path,
          winslash = "/",
          mustWork = TRUE
        )
      )
    )
  }

  required_summary_fields <- c(
    "version",
    "requested_g",
    "status",
    "package_converged",
    "stored_ng",
    "posterior_ncol",
    "occupied_groups",
    "structure_exact",
    "all_classes_occupied",
    "full_g",
    "empty_classes",
    "screening_tol",
    "nstart",
    "delete_empty",
    "technical_maxit",
    "automatic_retry",
    "em_iterations",
    "final_logLik",
    "AIC",
    "BIC",
    "CAIC",
    "ssBIC",
    "HQIC",
    "min_class_pct_nonempty",
    "min_APPA_nonempty",
    "mean_max_posterior",
    "median_max_posterior",
    "pct_max_posterior_ge_0_70",
    "pct_max_posterior_ge_0_90",
    "assignment_posterior_mismatch_n",
    "preferred_class_size_ok",
    "preferred_APPA_ok",
    "IC_comparison_eligible",
    "seed",
    "fit_seconds",
    "total_task_seconds",
    "n_units",
    "n_rows",
    "pid",
    "start_time",
    "end_time",
    "model_file",
    "resumed_skip",
    "error_message"
  )

  missing_summary_fields <- setdiff(
    required_summary_fields,
    names(old_summary)
  )

  if (length(missing_summary_fields) > 0L) {
    stop(
      paste0(
        "Existing V251 summary is missing exact expected field(s): ",
        paste(
          missing_summary_fields,
          collapse = ", "
        ),
        "\nFile:\n",
        normalizePath(
          summary_path,
          winslash = "/",
          mustWork = TRUE
        ),
        "\nStop here rather than guess a legacy schema."
      )
    )
  }

  if (
    as.integer(
      old_summary$requested_g[1L]
    ) !=
      as.integer(g)
  ) {
    stop(
      paste0(
        "Existing summary requested_g does not match filename G",
        g,
        "."
      )
    )
  }

  fit <- readRDS(
    fit_path
  )

  if (!inherits(fit, "gbmt")) {
    stop(
      paste0(
        "Existing model is not a gbmt object for G",
        g,
        ":\n",
        normalizePath(
          fit_path,
          winslash = "/",
          mustWork = TRUE
        )
      )
    )
  }

  required_fit_fields <- c(
    "call",
    "posterior",
    "assign",
    "em",
    "ic",
    "beta",
    "Sigma",
    "logLik"
  )

  missing_fit_fields <- setdiff(
    required_fit_fields,
    names(fit)
  )

  if (length(missing_fit_fields) > 0L) {
    stop(
      paste0(
        "Existing gbmt object is missing exact expected field(s): ",
        paste(
          missing_fit_fields,
          collapse = ", "
        ),
        "\nG",
        g,
        " will not be reinterpreted by guessing object structure."
      )
    )
  }

  if (!file.exists(summary_backup_path)) {
    summary_backup_ok <- file.copy(
      summary_path,
      summary_backup_path,
      overwrite = FALSE
    )

    if (!isTRUE(summary_backup_ok)) {
      stop(
        paste0(
          "Could not create legacy summary backup for G",
          g,
          ":\n",
          normalizePath(
            summary_backup_path,
            winslash = "/",
            mustWork = FALSE
          )
        )
      )
    }
  }

  if (
    file.exists(done_path) &&
    !file.exists(done_backup_path)
  ) {
    done_backup_ok <- file.copy(
      done_path,
      done_backup_path,
      overwrite = FALSE
    )

    if (!isTRUE(done_backup_ok)) {
      stop(
        paste0(
          "Could not create legacy DONE backup for G",
          g,
          ":\n",
          normalizePath(
            done_backup_path,
            winslash = "/",
            mustWork = FALSE
          )
        )
      )
    }
  }

  old_status <- as.character(
    old_summary$status[1L]
  )

  old_iterations <- suppressWarnings(
    as.numeric(
      old_summary$em_iterations[1L]
    )
  )

  em_dt <- extract_em_table(
    fit
  )

  em_status <- extract_em_status(
    fit
  )

  class_obj <- extract_class_quality(
    fit,
    requested_g = g
  )

  q <- class_obj$summary

  package_converged <- isTRUE(
    em_status$converged
  )

  structure_exact <- (
    q$structure_exact[1L] ==
      1L
  )

  all_occupied <- (
    q$all_classes_occupied[1L] ==
      1L
  )

  full_g <- (
    q$full_g[1L] ==
      1L
  )

  status <- if (
    package_converged &&
    full_g
  ) {

    "CONVERGED_VALID_FULL_G"

  } else if (
    package_converged &&
    structure_exact &&
    !all_occupied
  ) {

    "CONVERGED_EMPTY_CLASS"

  } else if (
    package_converged &&
    !structure_exact
  ) {

    "CONVERGED_STRUCTURE_MISMATCH"

  } else if (
    !package_converged &&
    full_g
  ) {

    "NONCONVERGED_FULL_G"

  } else if (
    !package_converged &&
    structure_exact &&
    !all_occupied
  ) {

    "NONCONVERGED_EMPTY_CLASS"

  } else {

    "NONCONVERGED_STRUCTURE_MISMATCH"

  }

  final_logLik <- suppressWarnings(
    as.numeric(
      logLik(fit)
    )
  )

  AIC_value <- safe_ic_value(
    fit,
    "AIC"
  )

  BIC_value <- safe_ic_value(
    fit,
    "BIC"
  )

  CAIC_value <- safe_ic_value(
    fit,
    "CAIC"
  )

  ssBIC_value <- safe_ic_value(
    fit,
    "ssBIC"
  )

  HQIC_value <- safe_ic_value(
    fit,
    "HQIC"
  )

  if (!is.finite(AIC_value)) {
    AIC_value <- tryCatch(
      as.numeric(
        AIC(fit)
      ),
      error = function(e) NA_real_
    )
  }

  if (!is.finite(BIC_value)) {
    BIC_value <- tryCatch(
      as.numeric(
        BIC(fit)
      ),
      error = function(e) NA_real_
    )
  }

  min_class_pct <- q$min_class_pct_nonempty[1L]
  min_appa <- q$min_APPA_nonempty[1L]

  preferred_class_size_ok <- is.finite(
    min_class_pct
  ) &&
    min_class_pct >=
      preferred_min_class_pct

  preferred_appa_ok <- is.finite(
    min_appa
  ) &&
    min_appa >=
      preferred_min_appa

  ic_comparison_eligible <- (
    status ==
      "CONVERGED_VALID_FULL_G"
  )

  # IMPORTANT: construct a brand-new one-row summary rather than mutating
  # old_summary with data.table :=.  In data.table j, a bare symbol that is
  # also a column name resolves to the column before the function-local
  # variable.  The V257 implementation therefore reassigned the old status
  # column to itself (e.g. NONCONVERGED_FULL_G) instead of using the freshly
  # computed local status.  Explicit scalar extraction below removes all
  # data.table non-standard-evaluation ambiguity and preserves the exact
  # historical V251 column order/schema.
  refreshed_summary <- data.table(
    version = as.character(VERSION_TAG),
    requested_g = as.integer(g),
    status = as.character(status),
    package_converged = as.integer(package_converged),
    stored_ng = as.integer(q$stored_ng[1L]),
    posterior_ncol = as.integer(q$posterior_ncol[1L]),
    occupied_groups = as.integer(q$occupied_groups[1L]),
    structure_exact = as.integer(structure_exact),
    all_classes_occupied = as.integer(all_occupied),
    full_g = as.integer(full_g),
    empty_classes = as.character(q$empty_classes[1L]),
    screening_tol = as.numeric(old_summary$screening_tol[1L]),
    nstart = as.integer(old_summary$nstart[1L]),
    delete_empty = as.integer(old_summary$delete_empty[1L]),
    technical_maxit = as.integer(old_summary$technical_maxit[1L]),
    automatic_retry = as.integer(old_summary$automatic_retry[1L]),
    em_iterations = as.numeric(em_status$iterations),
    final_logLik = as.numeric(final_logLik),
    AIC = as.numeric(AIC_value),
    BIC = as.numeric(BIC_value),
    CAIC = as.numeric(CAIC_value),
    ssBIC = as.numeric(ssBIC_value),
    HQIC = as.numeric(HQIC_value),
    min_class_pct_nonempty = as.numeric(min_class_pct),
    min_APPA_nonempty = as.numeric(min_appa),
    mean_max_posterior = as.numeric(q$mean_max_posterior[1L]),
    median_max_posterior = as.numeric(q$median_max_posterior[1L]),
    pct_max_posterior_ge_0_70 = as.numeric(
      q$pct_max_posterior_ge_0_70[1L]
    ),
    pct_max_posterior_ge_0_90 = as.numeric(
      q$pct_max_posterior_ge_0_90[1L]
    ),
    assignment_posterior_mismatch_n = as.integer(
      q$assignment_posterior_mismatch_n[1L]
    ),
    preferred_class_size_ok = as.integer(preferred_class_size_ok),
    preferred_APPA_ok = as.integer(preferred_appa_ok),
    IC_comparison_eligible = as.integer(ic_comparison_eligible),
    seed = as.integer(old_summary$seed[1L]),
    fit_seconds = as.numeric(old_summary$fit_seconds[1L]),
    total_task_seconds = as.numeric(old_summary$total_task_seconds[1L]),
    n_units = as.integer(old_summary$n_units[1L]),
    n_rows = as.integer(old_summary$n_rows[1L]),
    pid = as.integer(old_summary$pid[1L]),
    start_time = as.character(old_summary$start_time[1L]),
    end_time = as.character(old_summary$end_time[1L]),
    model_file = normalizePath(
      fit_path,
      winslash = "/",
      mustWork = TRUE
    ),
    resumed_skip = as.integer(old_summary$resumed_skip[1L]),
    error_message = as.character(old_summary$error_message[1L])
  )

  exact_summary_columns <- c(
    "version",
    "requested_g",
    "status",
    "package_converged",
    "stored_ng",
    "posterior_ncol",
    "occupied_groups",
    "structure_exact",
    "all_classes_occupied",
    "full_g",
    "empty_classes",
    "screening_tol",
    "nstart",
    "delete_empty",
    "technical_maxit",
    "automatic_retry",
    "em_iterations",
    "final_logLik",
    "AIC",
    "BIC",
    "CAIC",
    "ssBIC",
    "HQIC",
    "min_class_pct_nonempty",
    "min_APPA_nonempty",
    "mean_max_posterior",
    "median_max_posterior",
    "pct_max_posterior_ge_0_70",
    "pct_max_posterior_ge_0_90",
    "assignment_posterior_mismatch_n",
    "preferred_class_size_ok",
    "preferred_APPA_ok",
    "IC_comparison_eligible",
    "seed",
    "fit_seconds",
    "total_task_seconds",
    "n_units",
    "n_rows",
    "pid",
    "start_time",
    "end_time",
    "model_file",
    "resumed_skip",
    "error_message"
  )

  if (!identical(names(refreshed_summary), exact_summary_columns)) {
    stop(
      paste0(
        "Internal V258 summary-construction schema mismatch for G",
        g,
        ".\nConstructed columns: ",
        paste(names(refreshed_summary), collapse = ", "),
        "\nExpected columns: ",
        paste(exact_summary_columns, collapse = ", ")
      )
    )
  }

  stale_summary_temps <- Sys.glob(
    paste0(
      summary_path,
      ".tmp_*"
    )
  )

  if (length(stale_summary_temps) > 0L) {
    unlink(
      stale_summary_temps,
      force = TRUE
    )
  }

  tmp_summary <- paste0(
    summary_path,
    ".tmp_",
    Sys.getpid()
  )

  fwrite(
    refreshed_summary,
    tmp_summary
  )

  verify_summary <- read_v251_summary_exact(
    tmp_summary
  )

  verification_problems <- character()

  if (nrow(verify_summary) != 1L) {
    verification_problems <- c(
      verification_problems,
      paste0(
        "row_count=",
        nrow(verify_summary),
        " (expected 1)"
      )
    )
  }

  if (nrow(verify_summary) >= 1L) {

    verified_g <- suppressWarnings(
      as.integer(
        verify_summary$requested_g[1L]
      )
    )

    if (
      length(verified_g) != 1L ||
      is.na(verified_g) ||
      verified_g != as.integer(g)
    ) {
      verification_problems <- c(
        verification_problems,
        paste0(
          "requested_g=",
          as.character(
            verify_summary$requested_g[1L]
          ),
          " (expected ",
          g,
          ")"
        )
      )
    }

    verified_status <- as.character(
      verify_summary$status[1L]
    )

    if (!identical(
      verified_status,
      as.character(status)
    )) {
      verification_problems <- c(
        verification_problems,
        paste0(
          "status=",
          verified_status,
          " (expected ",
          status,
          ")"
        )
      )
    }

    verified_iterations <- suppressWarnings(
      as.numeric(
        verify_summary$em_iterations[1L]
      )
    )

    expected_iterations <- as.numeric(
      em_status$iterations
    )

    if (
      length(verified_iterations) != 1L ||
      length(expected_iterations) != 1L ||
      (
        is.na(verified_iterations) !=
          is.na(expected_iterations)
      ) ||
      (
        is.finite(verified_iterations) &&
        is.finite(expected_iterations) &&
        !isTRUE(all.equal(
          verified_iterations,
          expected_iterations,
          tolerance = 0
        ))
      )
    ) {
      verification_problems <- c(
        verification_problems,
        paste0(
          "em_iterations=",
          as.character(verified_iterations),
          " (expected ",
          as.character(expected_iterations),
          ")"
        )
      )
    }

    verified_empty <- verify_summary$empty_classes[1L]
    if (is.na(verified_empty)) {
      verified_empty <- ""
    }

    expected_empty <- as.character(
      q$empty_classes[1L]
    )
    if (is.na(expected_empty)) {
      expected_empty <- ""
    }

    if (!identical(
      as.character(verified_empty),
      expected_empty
    )) {
      verification_problems <- c(
        verification_problems,
        paste0(
          "empty_classes='",
          as.character(verified_empty),
          "' (expected '",
          expected_empty,
          "')"
        )
      )
    }
  }

  if (length(verification_problems) > 0L) {
    stop(
      paste0(
        "Refreshed summary verification failed for G",
        g,
        ". Details: ",
        paste(
          verification_problems,
          collapse = "; "
        ),
        "\nTemporary file retained for inspection:\n",
        normalizePath(
          tmp_summary,
          winslash = "/",
          mustWork = FALSE
        )
      )
    )
  }

  if (file.exists(summary_path)) {
    file.remove(
      summary_path
    )
  }

  if (!file.rename(
    tmp_summary,
    summary_path
  )) {
    if (!file.copy(
      tmp_summary,
      summary_path,
      overwrite = TRUE
    )) {
      stop(
        paste0(
          "Could not finalize refreshed V251 summary for G",
          g,
          "."
        )
      )
    }
    file.remove(
      tmp_summary
    )
  }

  if (!file.exists(class_path)) {
    fwrite(
      class_obj$class_table,
      class_path
    )
  }

  if (!file.exists(em_runs_path)) {
    fwrite(
      em_dt,
      em_runs_path
    )
  }

  if (file.exists(done_path)) {
    writeLines(
      c(
        paste0("Version: ", VERSION_TAG),
        paste0("Task: G", g),
        paste0("Status: ", status),
        paste0(
          "Package converged: ",
          package_converged
        ),
        paste0(
          "Requested G: ",
          g
        ),
        paste0(
          "Stored ng: ",
          q$stored_ng[1L]
        ),
        paste0(
          "Posterior columns: ",
          q$posterior_ncol[1L]
        ),
        paste0(
          "Occupied groups: ",
          q$occupied_groups[1L]
        ),
        paste0(
          "Empty classes: ",
          ifelse(
            nzchar(
              q$empty_classes[1L]
            ),
            q$empty_classes[1L],
            "none"
          )
        ),
        paste0(
          "EM iterations: ",
          em_status$iterations
        ),
        paste0(
          "Min class %: ",
          min_class_pct
        ),
        paste0(
          "Min APPA: ",
          min_appa
        ),
        paste0(
          "IC comparison eligible: ",
          ic_comparison_eligible
        ),
        paste0(
          "Screening tol: ",
          old_summary$screening_tol[1L]
        ),
        paste0(
          "nstart: ",
          old_summary$nstart[1L]
        ),
        "delete.empty: FALSE",
        "Automatic retry: NONE",
        paste0(
          "Original end: ",
          old_summary$end_time[1L]
        ),
        "Metadata refreshed by V258 from the saved gbmt fit; model was NOT rerun."
      ),
      done_path
    )
  }

  audit_dt <- data.table(
    engine = ENGINE_TAG,
    compatible_output_version = VERSION_TAG,
    requested_g = as.integer(g),
    old_status = old_status,
    refreshed_status = as.character(status),
    old_em_iterations = old_iterations,
    refreshed_em_iterations = as.numeric(
      em_status$iterations
    ),
    model_rerun = 0L,
    fit_file = normalizePath(
      fit_path,
      winslash = "/",
      mustWork = TRUE
    ),
    refreshed_at = format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  )

  audit_exists <- file.exists(
    audit_path
  )

  fwrite(
    audit_dt,
    audit_path,
    append = audit_exists,
    col.names = !audit_exists
  )

  completed_checkpoint_paths <- c(
    file.path(
      output_dir,
      paste0(
        "V252_checkpoint_A_G",
        g,
        ".rds"
      )
    ),
    file.path(
      output_dir,
      paste0(
        "V252_checkpoint_B_G",
        g,
        ".rds"
      )
    )
  )

  for (cp in completed_checkpoint_paths) {
    if (file.exists(cp)) {
      file.remove(cp)
    }
  }

  rm(
    fit,
    old_summary,
    refreshed_summary,
    verify_summary,
    class_obj,
    q,
    em_dt
  )

  gc(
    full = TRUE
  )

  invisible(TRUE)
}


# ============================================================
# 9. One-G worker
# ============================================================

run_one_g_v251 <- function(
  g,
  minimal_input_rds,
  output_dir,
  z_markers,
  model_degree,
  screen_tol,
  nstart,
  stagger_seconds_per_g,
  nonbinding_maxit,
  master_seed,
  expected_patients,
  expected_long_rows,
  preferred_min_class_pct,
  preferred_min_appa,
  minimum_g,
  input_md5,
  checkpoint_every,
  scheduler_state_path,
  scheduler_lock_dir,
  scheduler_event_log,
  scheduler_helper_ps1,
  scheduler_launch_lock_dir,
  scheduler_last_launch_path,
  scheduler_memory_log,
  memory_min_available_gb,
  memory_min_available_fraction,
  memory_max_commit_fraction,
  memory_poll_seconds,
  process_priority_class
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
  })

  data.table::setDTthreads(1L)

  fit_path <- file.path(
    output_dir,
    paste0(
      "V251_screen_fit_G",
      g,
      ".rds"
    )
  )

  partial_fit_path <- file.path(
    output_dir,
    paste0(
      "V251_screen_fit_partial_G",
      g,
      ".rds"
    )
  )

  summary_path <- file.path(
    output_dir,
    paste0(
      "V251_task_summary_G",
      g,
      ".csv"
    )
  )

  class_path <- file.path(
    output_dir,
    paste0(
      "V251_class_quality_G",
      g,
      ".csv"
    )
  )

  em_runs_path <- file.path(
    output_dir,
    paste0(
      "V251_EM_runs_G",
      g,
      ".csv"
    )
  )

  em_progress_path <- file.path(
    output_dir,
    paste0(
      "V251_EM_progress_G",
      g,
      ".log"
    )
  )

  stage_path <- file.path(
    output_dir,
    paste0(
      "V251_stage_G",
      g,
      ".log"
    )
  )

  running_path <- file.path(
    output_dir,
    paste0(
      "V251_RUNNING_G",
      g,
      ".txt"
    )
  )

  done_path <- file.path(
    output_dir,
    paste0(
      "V251_DONE_G",
      g,
      ".txt"
    )
  )

  error_path <- file.path(
    output_dir,
    paste0(
      "V251_ERROR_G",
      g,
      ".txt"
    )
  )


  checkpoint_a <- file.path(
    output_dir,
    paste0(
      "V252_checkpoint_A_G",
      g,
      ".rds"
    )
  )

  checkpoint_b <- file.path(
    output_dir,
    paste0(
      "V252_checkpoint_B_G",
      g,
      ".rds"
    )
  )

  resume_audit_path <- file.path(
    output_dir,
    paste0(
      "V258_resume_audit_G",
      g,
      ".csv"
    )
  )

  # ----------------------------------------------------------
  # Resume protection
  # ----------------------------------------------------------

  if (
    file.exists(summary_path) &&
    file.exists(fit_path)
  ) {

    old <- tryCatch(
      fread(summary_path),
      error = function(e) NULL
    )

    completed_statuses <- c(
      "CONVERGED_VALID_FULL_G",
      "CONVERGED_EMPTY_CLASS",
      "CONVERGED_STRUCTURE_MISMATCH",
      "NONCONVERGED_FULL_G",
      "NONCONVERGED_EMPTY_CLASS",
      "NONCONVERGED_STRUCTURE_MISMATCH"
    )

    if (
      !is.null(old) &&
      nrow(old) >= 1L &&
      "status" %in% names(old) &&
      as.character(old$status[1L]) %in% completed_statuses
    ) {

      old[
        ,
        resumed_skip := 1L
      ]

      return(old)

    }

  }

  cleanup_paths <- c(
    partial_fit_path,
    running_path,
    done_path,
    error_path
  )

  for (f in cleanup_paths) {
    if (file.exists(f)) {
      file.remove(f)
    }
  }

  append_stage(
    stage_path,
    paste0(
      "===== V258 execution session entered; PID=",
      Sys.getpid(),
      " ====="
    )
  )

  append_progress_line(
    em_progress_path,
    paste0(
      "===== V258 execution session entered; PID=",
      Sys.getpid(),
      " ====="
    )
  )

  # Memory gate acts only before this model loads its private input copy.
  scheduler_wait_for_memory(
    g = g,
    helper_path = scheduler_helper_ps1,
    memory_log = scheduler_memory_log,
    min_available_gb = memory_min_available_gb,
    min_available_fraction = memory_min_available_fraction,
    max_commit_fraction = memory_max_commit_fraction,
    poll_seconds = memory_poll_seconds
  )

  launch_meta <- scheduler_wait_for_launch_slot(
    g = g,
    spacing_seconds = stagger_seconds_per_g,
    launch_lock_dir = scheduler_launch_lock_dir,
    last_launch_path = scheduler_last_launch_path
  )

  assignment <- scheduler_register_task(
    g = g,
    pid = Sys.getpid(),
    state_path = scheduler_state_path,
    lock_dir = scheduler_lock_dir,
    event_log = scheduler_event_log,
    helper_path = scheduler_helper_ps1,
    priority_class = process_priority_class
  )

  scheduler_registered <- TRUE

  on.exit(
    {
      if (isTRUE(scheduler_registered)) {
        try(
          scheduler_release_task(
            g = g,
            pid = Sys.getpid(),
            state_path = scheduler_state_path,
            lock_dir = scheduler_lock_dir,
            event_log = scheduler_event_log
          ),
          silent = TRUE
        )
      }
    },
    add = TRUE
  )

  start_time <- Sys.time()
  stagger_seconds <- as.numeric(launch_meta$waited_seconds)

  scheduler_context <- list(
    enabled = TRUE,
    g = as.integer(g),
    pid = as.integer(Sys.getpid()),
    state_path = scheduler_state_path,
    lock_dir = scheduler_lock_dir,
    event_log = scheduler_event_log,
    helper_path = scheduler_helper_ps1,
    priority_class = process_priority_class
  )

  seed_used <- as.integer(
    master_seed +
      g * 100003L
  )


  checkpoint_signature <- make_checkpoint_signature(
    g = g,
    seed_used = seed_used,
    input_md5 = input_md5,
    z_markers = z_markers,
    model_degree = model_degree,
    screen_tol = screen_tol,
    nonbinding_maxit = nonbinding_maxit
  )

  worker_process_identity <- read_ps_csv(
    run_ps_helper(
      mode = "process",
      target_pid = Sys.getpid(),
      helper_path = scheduler_helper_ps1
    )
  )

  required_worker_identity_fields <- c(
    "PID",
    "ProcessName",
    "ExecutablePath",
    "StartTimeLocal",
    "CommandLine",
    "Status"
  )

  if (
    nrow(worker_process_identity) != 1L ||
    !all(
      required_worker_identity_fields %in%
        names(worker_process_identity)
    ) ||
    !identical(
      as.character(worker_process_identity$Status[1L]),
      "FOUND"
    )
  ) {
    stop(
      paste0(
        "Could not record exact V258 worker process identity for G",
        g,
        " PID ",
        Sys.getpid(),
        ". Stop here rather than create an ambiguous RUNNING marker."
      )
    )
  }

  worker_process_name <- as.character(
    worker_process_identity$ProcessName[1L]
  )
  worker_executable_path <- as.character(
    worker_process_identity$ExecutablePath[1L]
  )
  worker_process_start_local <- as.character(
    worker_process_identity$StartTimeLocal[1L]
  )

  if (is.na(worker_process_name)) worker_process_name <- ""
  if (is.na(worker_executable_path)) worker_executable_path <- ""
  if (is.na(worker_process_start_local)) worker_process_start_local <- ""

  if (!nzchar(worker_process_start_local)) {
    stop(
      paste0(
        "Windows did not expose the exact process start time for V258 worker G",
        g,
        " PID ",
        Sys.getpid(),
        ". Stop here rather than create an ambiguous RUNNING marker."
      )
    )
  }

  writeLines(
    c(
      paste0("Compatible output version: V251"),
      paste0("Engine: V258 dynamic scheduler + V252 checkpoint/resume"),
      paste0("Task: G", g),
      "Status: RUNNING",
      paste0(
        "Start: ",
        format(
          start_time,
          "%Y-%m-%d %H:%M:%S"
        )
      ),
      paste0(
        "PID: ",
        Sys.getpid()
      ),
      paste0(
        "Worker process name: ",
        worker_process_name
      ),
      paste0(
        "Worker executable: ",
        worker_executable_path
      ),
      paste0(
        "Worker process start local: ",
        worker_process_start_local
      ),
      paste0(
        "Initial CPU slot: ",
        assignment$slot_id,
        "; tier=",
        assignment$tier_name,
        "; logicalCPU=",
        assignment$logical_index
      ),
      paste0(
        "Random seed: ",
        seed_used
      ),
      paste0(
        "tol: ",
        screen_tol
      ),
      paste0(
        "nstart: ",
        nstart
      ),
      "delete.empty: FALSE",
      paste0(
        "technical maxit: ",
        nonbinding_maxit
      ),
      paste0(
        "global launch spacing: ",
        stagger_seconds_per_g,
        " sec; actual wait in launch gate: ",
        round(stagger_seconds, 3),
        " sec"
      ),
      "automatic retry: NONE",
      "checkpoint/resume: ENABLED",
      "dynamic CPU promotion: ENABLED at checkpointed EM boundaries",
      paste0("checkpoint every accepted iteration: ", checkpoint_every)
    ),
    running_path
  )

  append_stage(
    stage_path,
    paste0(
      "worker assigned CPU slot ",
      assignment$slot_id,
      " tier=",
      assignment$tier_name,
      " logicalCPU=",
      assignment$logical_index
    )
  )

  append_stage(
    stage_path,
    "before readRDS minimal input"
  )

  model_data <- readRDS(
    minimal_input_rds
  )

  append_stage(
    stage_path,
    "after readRDS minimal input"
  )

  if (!is.data.frame(model_data)) {
    model_data <- as.data.frame(model_data)
  }

  model_data$unit_id <- as.character(
    model_data$unit_id
  )

  model_data$time_model <- as.numeric(
    model_data$time_model
  )

  for (v in z_markers) {
    model_data[[v]] <- as.numeric(
      model_data[[v]]
    )
  }

  n_units <- length(
    unique(
      model_data$unit_id
    )
  )

  n_model_rows <- nrow(
    model_data
  )

  if (n_units != expected_patients) {
    stop(
      paste0(
        "G",
        g,
        " worker patient count = ",
        n_units,
        "; expected ",
        expected_patients,
        "."
      )
    )
  }

  if (n_model_rows != expected_long_rows) {
    stop(
      paste0(
        "G",
        g,
        " worker row count = ",
        n_model_rows,
        "; expected ",
        expected_long_rows,
        "."
      )
    )
  }

  append_stage(
    stage_path,
    "worker input validation passed"
  )

  gc(full = TRUE)

  set.seed(
    seed_used
  )

  cat(
    paste0(
      "V258 engine / V251-compatible G",
      g,
      " fitting started at ",
      format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      ),
      "\n"
    ),
    file = em_progress_path,
    append = TRUE
  )

  append_stage(
    stage_path,
    "before gbmt call"
  )

  fit_error <- NA_character_
  fit <- NULL

  fit_start <- Sys.time()

  resume_meta <- list(
    resumed = 0L,
    resumed_from_iteration = 0L,
    final_em_iteration = NA_integer_
  )

  fit_result <- tryCatch(
    {

      gbmt_checkpoint_one_start_v252(
        x.names = z_markers,
        unit = "unit_id",
        time = "time_model",
        ng = as.integer(g),
        d = as.integer(model_degree),
        data = model_data,
        scaling = 0,
        pruning = TRUE,
        delete.empty = FALSE,
        nstart = as.integer(nstart),
        tol = as.numeric(screen_tol),
        maxit = as.integer(nonbinding_maxit),
        quiet = TRUE,
        checkpoint_a = checkpoint_a,
        checkpoint_b = checkpoint_b,
        checkpoint_signature = checkpoint_signature,
        checkpoint_every = as.integer(checkpoint_every),
        em_progress_path = em_progress_path,
        stage_path = stage_path,
        scheduler_context = scheduler_context
      )

    },
    error = function(e) {

      fit_error <<- conditionMessage(e)

      NULL

    }
  )

  if (!is.null(fit_result)) {
    fit <- fit_result$fit
    resume_meta <- fit_result$resume_meta
  }


  fit_end <- Sys.time()

  fit_seconds <- as.numeric(
    difftime(
      fit_end,
      fit_start,
      units = "secs"
    )
  )

  if (is.null(fit)) {

    append_stage(
      stage_path,
      paste0(
        "gbmt ERROR: ",
        fit_error
      )
    )

    error_dt <- data.table(
      requested_g = as.integer(g),
      status = "ERROR",
      error_message = as.character(fit_error),
      start_time = format(
        start_time,
        "%Y-%m-%d %H:%M:%S"
      ),
      end_time = format(
        Sys.time(),
        "%Y-%m-%d %H:%M:%S"
      ),
      fit_seconds = as.numeric(fit_seconds),
      seed = as.integer(seed_used)
    )

    fwrite(
      error_dt,
      summary_path
    )

    if (file.exists(running_path)) {
      file.remove(running_path)
    }

    writeLines(
      c(
        paste0("Version: V251"),
        paste0("Task: G", g),
        "Status: ERROR",
        paste0(
          "Error: ",
          fit_error
        ),
        paste0(
          "End: ",
          format(
            Sys.time(),
            "%Y-%m-%d %H:%M:%S"
          )
        )
      ),
      error_path
    )

    return(error_dt)

  }

  append_stage(
    stage_path,
    "gbmt call returned"
  )

  # ----------------------------------------------------------
  # Save raw fit BEFORE classification checks.
  # This preserves the exact requested-G fit even if an empty
  # modal class is present.
  # ----------------------------------------------------------

  saveRDS(
    fit,
    partial_fit_path
  )

  if (
    !file.exists(partial_fit_path) ||
    is.na(
      file.info(partial_fit_path)$size
    ) ||
    file.info(partial_fit_path)$size <= 0
  ) {
    stop(
      paste0(
        "G",
        g,
        " model partial RDS was not written correctly."
      )
    )
  }

  if (file.exists(fit_path)) {
    file.remove(fit_path)
  }

  renamed <- file.rename(
    partial_fit_path,
    fit_path
  )

  if (!renamed) {

    copied <- file.copy(
      partial_fit_path,
      fit_path,
      overwrite = TRUE
    )

    if (!copied) {
      stop(
        paste0(
          "Could not finalize model RDS for G",
          g,
          "."
        )
      )
    }

    file.remove(
      partial_fit_path
    )

  }

  # ----------------------------------------------------------
  # Extract EM and class diagnostics
  # ----------------------------------------------------------

  em_dt <- extract_em_table(
    fit
  )

  fwrite(
    em_dt,
    em_runs_path
  )

  em_status <- extract_em_status(
    fit
  )

  class_obj <- extract_class_quality(
    fit,
    requested_g = g
  )

  fwrite(
    class_obj$class_table,
    class_path
  )

  q <- class_obj$summary

  package_converged <- isTRUE(
    em_status$converged
  )

  structure_exact <- (
    q$structure_exact[1L] ==
      1L
  )

  all_occupied <- (
    q$all_classes_occupied[1L] ==
      1L
  )

  full_g <- (
    q$full_g[1L] ==
      1L
  )

  # ----------------------------------------------------------
  # Status logic
  # ----------------------------------------------------------

  status <- if (
    package_converged &&
    full_g
  ) {

    "CONVERGED_VALID_FULL_G"

  } else if (
    package_converged &&
    structure_exact &&
    !all_occupied
  ) {

    "CONVERGED_EMPTY_CLASS"

  } else if (
    package_converged &&
    !structure_exact
  ) {

    "CONVERGED_STRUCTURE_MISMATCH"

  } else if (
    !package_converged &&
    full_g
  ) {

    "NONCONVERGED_FULL_G"

  } else if (
    !package_converged &&
    structure_exact &&
    !all_occupied
  ) {

    "NONCONVERGED_EMPTY_CLASS"

  } else {

    "NONCONVERGED_STRUCTURE_MISMATCH"

  }

  # ----------------------------------------------------------
  # Information criteria
  # ----------------------------------------------------------

  final_logLik <- suppressWarnings(
    as.numeric(
      logLik(fit)
    )
  )

  AIC_value <- safe_ic_value(
    fit,
    "AIC"
  )

  BIC_value <- safe_ic_value(
    fit,
    "BIC"
  )

  CAIC_value <- safe_ic_value(
    fit,
    "CAIC"
  )

  ssBIC_value <- safe_ic_value(
    fit,
    "ssBIC"
  )

  HQIC_value <- safe_ic_value(
    fit,
    "HQIC"
  )

  if (!is.finite(AIC_value)) {
    AIC_value <- tryCatch(
      as.numeric(AIC(fit)),
      error = function(e) NA_real_
    )
  }

  if (!is.finite(BIC_value)) {
    BIC_value <- tryCatch(
      as.numeric(BIC(fit)),
      error = function(e) NA_real_
    )
  }

  min_class_pct <- q$min_class_pct_nonempty[1L]
  min_appa <- q$min_APPA_nonempty[1L]

  preferred_class_size_ok <- is.finite(
    min_class_pct
  ) &&
    min_class_pct >=
      preferred_min_class_pct

  preferred_appa_ok <- is.finite(
    min_appa
  ) &&
    min_appa >=
      preferred_min_appa

  ic_comparison_eligible <- (
    status ==
      "CONVERGED_VALID_FULL_G"
  )

  end_time <- Sys.time()

  summary_dt <- data.table(
    version = VERSION_TAG,
    requested_g = as.integer(g),
    status = as.character(status),
    package_converged = as.integer(package_converged),
    stored_ng = as.integer(q$stored_ng[1L]),
    posterior_ncol = as.integer(q$posterior_ncol[1L]),
    occupied_groups = as.integer(q$occupied_groups[1L]),
    structure_exact = as.integer(structure_exact),
    all_classes_occupied = as.integer(all_occupied),
    full_g = as.integer(full_g),
    empty_classes = as.character(q$empty_classes[1L]),
    screening_tol = as.numeric(screen_tol),
    nstart = as.integer(nstart),
    delete_empty = 0L,
    technical_maxit = as.integer(nonbinding_maxit),
    automatic_retry = 0L,
    em_iterations = as.numeric(em_status$iterations),
    final_logLik = as.numeric(final_logLik),
    AIC = as.numeric(AIC_value),
    BIC = as.numeric(BIC_value),
    CAIC = as.numeric(CAIC_value),
    ssBIC = as.numeric(ssBIC_value),
    HQIC = as.numeric(HQIC_value),
    min_class_pct_nonempty = as.numeric(min_class_pct),
    min_APPA_nonempty = as.numeric(min_appa),
    mean_max_posterior = as.numeric(
      q$mean_max_posterior[1L]
    ),
    median_max_posterior = as.numeric(
      q$median_max_posterior[1L]
    ),
    pct_max_posterior_ge_0_70 = as.numeric(
      q$pct_max_posterior_ge_0_70[1L]
    ),
    pct_max_posterior_ge_0_90 = as.numeric(
      q$pct_max_posterior_ge_0_90[1L]
    ),
    assignment_posterior_mismatch_n = as.integer(
      q$assignment_posterior_mismatch_n[1L]
    ),
    preferred_class_size_ok = as.integer(
      preferred_class_size_ok
    ),
    preferred_APPA_ok = as.integer(
      preferred_appa_ok
    ),
    IC_comparison_eligible = as.integer(
      ic_comparison_eligible
    ),
    seed = as.integer(seed_used),
    fit_seconds = as.numeric(fit_seconds),
    total_task_seconds = as.numeric(
      difftime(
        end_time,
        start_time,
        units = "secs"
      )
    ),
    n_units = as.integer(n_units),
    n_rows = as.integer(n_model_rows),
    pid = as.integer(Sys.getpid()),
    start_time = format(
      start_time,
      "%Y-%m-%d %H:%M:%S"
    ),
    end_time = format(
      end_time,
      "%Y-%m-%d %H:%M:%S"
    ),
    model_file = normalizePath(
      fit_path,
      winslash = "/",
      mustWork = TRUE
    ),
    resumed_skip = 0L,
    error_message = NA_character_
  )

  fwrite(
    summary_dt,
    summary_path
  )

  resume_audit_dt <- data.table(
    engine = "V258",
    compatible_output_version = "V251",
    requested_g = as.integer(g),
    resumed = as.integer(resume_meta$resumed),
    resumed_from_iteration = as.integer(resume_meta$resumed_from_iteration),
    final_em_iteration = as.integer(resume_meta$final_em_iteration),
    checkpoint_every_accepted_iteration = as.integer(checkpoint_every),
    checkpoint_signature = as.character(checkpoint_signature),
    input_md5 = as.character(input_md5),
    completed_at = format(
      Sys.time(),
      "%Y-%m-%d %H:%M:%S"
    )
  )

  fwrite(
    resume_audit_dt,
    resume_audit_path
  )

  # Final compatible RDS + summary are safely written. Checkpoints are
  # no longer needed and are removed only at this point.
  for (cp in c(checkpoint_a, checkpoint_b)) {
    if (file.exists(cp)) {
      file.remove(cp)
    }
  }

  if (file.exists(running_path)) {
    file.remove(running_path)
  }

  done_lines <- c(
    paste0("Version: ", VERSION_TAG),
    paste0("Task: G", g),
    paste0("Status: ", status),
    paste0(
      "Package converged: ",
      package_converged
    ),
    paste0(
      "Requested G: ",
      g
    ),
    paste0(
      "Stored ng: ",
      q$stored_ng[1L]
    ),
    paste0(
      "Posterior columns: ",
      q$posterior_ncol[1L]
    ),
    paste0(
      "Occupied groups: ",
      q$occupied_groups[1L]
    ),
    paste0(
      "Empty classes: ",
      ifelse(
        nzchar(q$empty_classes[1L]),
        q$empty_classes[1L],
        "none"
      )
    ),
    paste0(
      "EM iterations: ",
      em_status$iterations
    ),
    paste0(
      "Min class %: ",
      min_class_pct
    ),
    paste0(
      "Min APPA: ",
      min_appa
    ),
    paste0(
      "IC comparison eligible: ",
      ic_comparison_eligible
    ),
    paste0(
      "Screening tol: ",
      screen_tol
    ),
    paste0(
      "nstart: ",
      nstart
    ),
    "delete.empty: FALSE",
    "Automatic retry: NONE",
    paste0(
      "End: ",
      format(
        end_time,
        "%Y-%m-%d %H:%M:%S"
      )
    )
  )

  writeLines(
    done_lines,
    done_path
  )

  append_stage(
    stage_path,
    paste0(
      "task finished: ",
      status
    )
  )

  summary_dt

}


# ============================================================
# 10. Refresh existing V251 metadata, then determine pending G
# ============================================================

cat("Checking existing V251-compatible completed fits...\n")

legacy_refresh_result <- vapply(
  GROUPS_TO_SCREEN,
  refresh_existing_v251_metadata,
  logical(1),
  output_dir = output_dir,
  preferred_min_class_pct = PREFERRED_MIN_CLASS_PCT,
  preferred_min_appa = PREFERRED_MIN_APPA
)

refreshed_groups <- GROUPS_TO_SCREEN[
  legacy_refresh_result
]

if (length(refreshed_groups) > 0L) {
  cat(
    "Existing fit metadata verified/refreshed for: ",
    paste0(
      "G",
      refreshed_groups,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
} else {
  cat("No existing complete fit required metadata refresh.\n")
}

existing_complete_groups <- GROUPS_TO_SCREEN[
  vapply(
    GROUPS_TO_SCREEN,
    task_is_complete,
    logical(1),
    output_dir = output_dir
  )
]

if (length(existing_complete_groups) > 0L) {
  cat(
    "Existing completed G values to skip: ",
    paste0(
      "G",
      existing_complete_groups,
      collapse = ", "
    ),
    "\n",
    sep = ""
  )
} else {
  cat("Existing completed G values to skip: none\n")
}

pending_groups <- GROUPS_TO_SCREEN[
  !vapply(
    GROUPS_TO_SCREEN,
    task_is_complete,
    logical(1),
    output_dir = output_dir
  )
]

cat("Pending G values: ")

if (length(pending_groups) == 0L) {

  cat("none\n\n")

} else {

  cat(
    paste0(
      "G",
      pending_groups,
      collapse = ", "
    ),
    "\n\n",
    sep = ""
  )

}


# ============================================================
# 11. Master-run markers
# ============================================================

master_start <- Sys.time()

master_running_path <- file.path(
  output_dir,
  "V251_MASTER_RUNNING.txt"
)

master_done_path <- file.path(
  output_dir,
  "V251_MASTER_DONE.txt"
)

if (file.exists(master_done_path)) {
  file.remove(master_done_path)
}

writeLines(
  c(
    "V258 dynamic scheduler + V252 checkpoint/resume engine is active (V251-compatible outputs).",
    paste0(
      "Start: ",
      format(
        master_start,
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    "Joint model: MSI + Na + K + Cl + HCO3",
    "Requested G: 5-20",
    paste0("Detected physical cores: ", physical_cores),
    paste0("Detected logical CPUs: ", logical_cores),
    paste0("Maximum active tasks from hardware: ", MAX_ACTIVE_TASKS),
    "Pending models are distributed one G per PSOCK worker until logical CPU capacity is reached.",
    paste0("Global actual-start stagger: ", START_STAGGER_SECONDS, " sec"),
    "CPU order: P > E > LP-E > SMT",
    "Promotion: only after accepted EM iteration + successful checkpoint",
    "Promotion priority: low iteration, then high recent deltaLL, then long iteration time",
    "RAM-aware launch guard: ENABLED",
    "tol: 10",
    "nstart: 1",
    "delete.empty: FALSE",
    "Automatic retry: NONE",
    "EM checkpoint/resume: ENABLED after every accepted iteration",
    "Dual-slot power-safe checkpoints: ENABLED",
    "Final files remain V251-compatible",
    "Wall-clock timeout: NONE",
    "Practical maxit budget: NONE"
  ),
  master_running_path
)


# ============================================================
# 12. Run pending models in parallel
# ============================================================

if (length(pending_groups) > 0L) {

  workers_now <- min(
    MAX_ACTIVE_TASKS,
    length(pending_groups)
  )

  cat(
    "Starting ",
    workers_now,
    " PSOCK workers...\n\n",
    sep = ""
  )

  worker_console_path <- file.path(
    output_dir,
    "V251_PS0CK_worker_console.log"
  )

  cl <- parallel::makePSOCKcluster(
    workers_now,
    outfile = worker_console_path
  )

  cluster_stopped <- FALSE

  on.exit(
    {

      if (
        exists("cl") &&
        !is.null(cl) &&
        !cluster_stopped
      ) {

        try(
          parallel::stopCluster(cl),
          silent = TRUE
        )

      }

    },
    add = TRUE
  )

  parallel::clusterExport(
    cl,
    varlist = c(
      "run_one_g_v251",
      "append_stage",
      "safe_ic_value",
      "extract_em_table",
      "extract_em_status",
      "extract_class_quality",
      "make_checkpoint_signature",
      "append_progress_line",
      "validate_checkpoint_state",
      "read_checkpoint_slot",
      "load_latest_checkpoint",
      "save_checkpoint_dual_slot",
      "make_checkpoint_state",
      "gbmt_fit_checkpoint_v252",
      "gbmt_checkpoint_one_start_v252",
      "scheduler_acquire_lock",
      "scheduler_release_lock",
      "scheduler_save_state_locked",
      "scheduler_log_event_locked",
      "scheduler_apply_affinity",
      "scheduler_wait_for_memory",
      "scheduler_wait_for_launch_slot",
      "scheduler_reserve_free_slots_locked",
      "scheduler_register_task",
      "scheduler_iteration_boundary",
      "scheduler_release_task",
      "query_memory_status",
      "run_ps_helper",
      "read_ps_csv",
      "powershell_exe",
      "PROCESS_PRIORITY_CLASS",
      "scheduler_helper_ps1",
      "CHECKPOINT_STATE_VERSION",
      "VERSION_TAG"
    ),
    envir = environment()
  )

  results_new <- parallel::parLapplyLB(
    cl,
    pending_groups,
    fun = run_one_g_v251,
    minimal_input_rds = minimal_input_rds,
    output_dir = output_dir,
    z_markers = Z_MARKERS,
    model_degree = MODEL_DEGREE,
    screen_tol = SCREEN_TOL,
    nstart = NSTART,
    stagger_seconds_per_g = START_STAGGER_SECONDS,
    nonbinding_maxit = NONBINDING_MAXIT,
    master_seed = MASTER_SEED,
    expected_patients = EXPECTED_PATIENTS,
    expected_long_rows = EXPECTED_LONG_ROWS,
    preferred_min_class_pct = PREFERRED_MIN_CLASS_PCT,
    preferred_min_appa = PREFERRED_MIN_APPA,
    minimum_g = min(GROUPS_TO_SCREEN),
    input_md5 = input_md5,
    checkpoint_every = CHECKPOINT_EVERY_ACCEPTED_ITER,
    scheduler_state_path = scheduler_state_path,
    scheduler_lock_dir = scheduler_lock_dir,
    scheduler_event_log = scheduler_event_log,
    scheduler_helper_ps1 = scheduler_helper_ps1,
    scheduler_launch_lock_dir = scheduler_launch_lock_dir,
    scheduler_last_launch_path = scheduler_last_launch_path,
    scheduler_memory_log = scheduler_memory_log,
    memory_min_available_gb = MEMORY_MIN_AVAILABLE_GB,
    memory_min_available_fraction = MEMORY_MIN_AVAILABLE_FRACTION,
    memory_max_commit_fraction = MEMORY_MAX_COMMIT_FRACTION,
    memory_poll_seconds = MEMORY_POLL_SECONDS,
    process_priority_class = PROCESS_PRIORITY_CLASS
  )

  parallel::stopCluster(
    cl
  )

  cluster_stopped <- TRUE
  cl <- NULL

} else {

  results_new <- list()

}


# ============================================================
# 13. Aggregate per-G summaries
# ============================================================

summary_files <- file.path(
  output_dir,
  paste0(
    "V251_task_summary_G",
    GROUPS_TO_SCREEN,
    ".csv"
  )
)

summary_list <- lapply(
  seq_along(summary_files),
  function(i) {

    f <- summary_files[i]
    g <- GROUPS_TO_SCREEN[i]

    if (!file.exists(f)) {

      return(
        data.table(
          version = VERSION_TAG,
          requested_g = as.integer(g),
          status = "NOT_FINISHED_OR_MISSING"
        )
      )

    }

    tryCatch(
      {

        fread(f)

      },
      error = function(e) {

        data.table(
          version = VERSION_TAG,
          requested_g = as.integer(g),
          status = "SUMMARY_READ_ERROR",
          error_message = conditionMessage(e)
        )

      }
    )

  }
)

comparison_dt <- rbindlist(
  summary_list,
  fill = TRUE
)

setorder(
  comparison_dt,
  requested_g
)

fwrite(
  comparison_dt,
  file.path(
    output_dir,
    "V251_model_comparison.csv"
  )
)


# ============================================================
# 14. Compact run-status table
# ============================================================

desired_status_cols <- c(
  "version",
  "requested_g",
  "status",
  "package_converged",
  "stored_ng",
  "posterior_ncol",
  "occupied_groups",
  "empty_classes",
  "em_iterations",
  "final_logLik",
  "BIC",
  "min_class_pct_nonempty",
  "min_APPA_nonempty",
  "IC_comparison_eligible",
  "fit_seconds",
  "start_time",
  "end_time"
)

existing_status_cols <- intersect(
  desired_status_cols,
  names(comparison_dt)
)

run_status_dt <- comparison_dt[
  ,
  ..existing_status_cols
]

fwrite(
  run_status_dt,
  file.path(
    output_dir,
    "V251_run_status.csv"
  )
)


# ============================================================
# 15. Preliminary ranking
# ============================================================

if (
  all(
    c(
      "status",
      "BIC"
    ) %in%
      names(comparison_dt)
  )
) {

  preliminary_ranking_dt <- comparison_dt[
    status ==
      "CONVERGED_VALID_FULL_G"
  ]

  if (nrow(preliminary_ranking_dt) > 0L) {

    setorder(
      preliminary_ranking_dt,
      BIC,
      requested_g,
      na.last = TRUE
    )

  }

} else {

  preliminary_ranking_dt <- data.table()

}

fwrite(
  preliminary_ranking_dt,
  file.path(
    output_dir,
    "V251_preliminary_ranking.csv"
  )
)


# ============================================================
# 16. Adjacent-G IC changes
# ============================================================

adjacent_dt <- data.table()

if (
  nrow(comparison_dt) > 0L &&
  all(
    c(
      "requested_g",
      "status",
      "BIC",
      "AIC",
      "CAIC",
      "ssBIC",
      "HQIC"
    ) %in%
      names(comparison_dt)
  )
) {

  valid_for_adjacent <- comparison_dt[
    status ==
      "CONVERGED_VALID_FULL_G"
  ][
    order(requested_g)
  ]

  if (nrow(valid_for_adjacent) > 0L) {

    adjacent_dt <- valid_for_adjacent[
      ,
      .(
        requested_g,
        BIC,
        delta_BIC_vs_previous_valid_G = BIC - shift(BIC),
        AIC,
        delta_AIC_vs_previous_valid_G = AIC - shift(AIC),
        CAIC,
        delta_CAIC_vs_previous_valid_G = CAIC - shift(CAIC),
        ssBIC,
        delta_ssBIC_vs_previous_valid_G = ssBIC - shift(ssBIC),
        HQIC,
        delta_HQIC_vs_previous_valid_G = HQIC - shift(HQIC)
      )
    ]

  }

}

fwrite(
  adjacent_dt,
  file.path(
    output_dir,
    "V251_adjacent_IC_change.csv"
  )
)


# ============================================================
# 17. Human-readable report
# ============================================================

master_end <- Sys.time()

elapsed_hours <- as.numeric(
  difftime(
    master_end,
    master_start,
    units = "hours"
  )
)

n_valid <- if (
  "status" %in%
    names(comparison_dt)
) {
  comparison_dt[
    status ==
      "CONVERGED_VALID_FULL_G",
    .N
  ]
} else {
  0L
}

n_empty <- if (
  "status" %in%
    names(comparison_dt)
) {
  comparison_dt[
    grepl(
      "EMPTY_CLASS",
      status,
      fixed = TRUE
    ),
    .N
  ]
} else {
  0L
}

n_structure_mismatch <- if (
  "status" %in%
    names(comparison_dt)
) {
  comparison_dt[
    grepl(
      "STRUCTURE_MISMATCH",
      status,
      fixed = TRUE
    ),
    .N
  ]
} else {
  0L
}

n_error <- if (
  "status" %in%
    names(comparison_dt)
) {
  comparison_dt[
    status ==
      "ERROR",
    .N
  ]
} else {
  0L
}

report_lines <- c(
  "V258 DYNAMIC-SCHEDULER CHECKPOINT/RESUME GBMTM G5-G20 (V251-COMPATIBLE OUTPUTS)",
  "================================================",
  "",
  paste0(
    "Start: ",
    format(
      master_start,
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  paste0(
    "End: ",
    format(
      master_end,
      "%Y-%m-%d %H:%M:%S"
    )
  ),
  paste0(
    "Elapsed hours: ",
    round(
      elapsed_hours,
      3
    )
  ),
  "",
  "MODEL",
  "MSI + Na + K + Cl + HCO3 jointly modeled in one GBMTM.",
  "",
  "SCREEN SETTINGS",
  "Requested G: G5-G20",
  paste0("Hardware-derived max active tasks: ", MAX_ACTIVE_TASKS),
  "Pending models are distributed one G per PSOCK worker.",
  paste0("Global actual-start stagger: ", START_STAGGER_SECONDS, " sec"),
  "CPU order: P > E > LP-E > SMT",
  "Safe-boundary dynamic promotion: ENABLED",
  "RAM-aware launch guard: ENABLED",
  "Polynomial degree: 2",
  "scaling: 0",
  "pruning: TRUE",
  "tol: 10",
  "nstart: 1",
  "delete.empty: FALSE",
  "Automatic seed retry: NONE",
  "EM checkpoint/resume: every accepted iteration",
  "Dual checkpoint slots: A/B",
  "Final model/CSV filenames: V251-compatible",
  "Existing V251_config.csv: preserved; V258 execution metadata is separate",
  "Legacy V251 summary status: re-read from saved gbmt fit before skip/ranking",
  "Wall-clock timeout: NONE",
  "Small/G-specific maxit budget: NONE",
  paste0(
    "Technical maxit: ",
    NONBINDING_MAXIT
  ),
  "",
  "RESULT COUNTS",
  paste0(
    "Converged valid full-G: ",
    n_valid,
    "/",
    length(GROUPS_TO_SCREEN)
  ),
  paste0(
    "Empty-class fits: ",
    n_empty
  ),
  paste0(
    "Structure-mismatch fits: ",
    n_structure_mismatch
  ),
  paste0(
    "Errors: ",
    n_error
  ),
  "",
  "INTERPRETATION RULE",
  "Only CONVERGED_VALID_FULL_G fits are eligible for preliminary IC ranking.",
  "CONVERGED_EMPTY_CLASS means the requested number of components was retained",
  "but one or more modal groups had zero assigned patients.",
  "Such a fit is saved exactly as estimated and is NOT relabeled as G-1.",
  "No automatic random-seed retry is performed.",
  "",
  "IMPORTANT",
  "tol=10 is deliberately relaxed for broad first-stage screening.",
  "Do not use this run alone to lock the final class number.",
  "Shortlisted G values require strict multi-start refinement.",
  "",
  "KEY FILES",
  "V251_model_comparison.csv",
  "V251_run_status.csv",
  "V251_preliminary_ranking.csv",
  "V251_adjacent_IC_change.csv",
  "V251_screen_report.txt"
)

writeLines(
  report_lines,
  file.path(
    output_dir,
    "V251_screen_report.txt"
  )
)


# ============================================================
# 18. Master completion marker
# ============================================================

if (file.exists(master_running_path)) {
  file.remove(master_running_path)
}

master_status <- if (
  n_valid ==
    length(GROUPS_TO_SCREEN)
) {

  "ALL_G5_G20_CONVERGED_VALID_FULL_G"

} else {

  "G5_G20_FINISHED_WITH_ITEMS_TO_REVIEW"

}

writeLines(
  c(
    paste0(
      "Version: ",
      VERSION_TAG
    ),
    paste0(
      "Status: ",
      master_status
    ),
    paste0(
      "Start: ",
      format(
        master_start,
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    paste0(
      "End: ",
      format(
        master_end,
        "%Y-%m-%d %H:%M:%S"
      )
    ),
    paste0(
      "Elapsed hours: ",
      round(
        elapsed_hours,
        3
      )
    ),
    paste0(
      "Converged valid full-G: ",
      n_valid,
      "/",
      length(GROUPS_TO_SCREEN)
    ),
    paste0(
      "Empty-class fits: ",
      n_empty
    ),
    paste0(
      "Structure-mismatch fits: ",
      n_structure_mismatch
    ),
    paste0(
      "Errors: ",
      n_error
    )
  ),
  master_done_path
)


# ============================================================
# 19. Final console output
# ============================================================

cat("\n")
cat("============================================================\n")
cat("V258 DYNAMIC-SCHEDULER CHECKPOINT/RESUME G5-G20 SCREEN FINISHED\n")
cat("============================================================\n\n")

print(
  run_status_dt
)

cat("\nPreliminary ranking among CONVERGED_VALID_FULL_G only:\n\n")

print(
  preliminary_ranking_dt
)

cat("\nAdjacent-G IC changes among valid full-G models:\n\n")

print(
  adjacent_dt
)

cat("\nOutput folder:\n")
cat(
  normalizePath(
    output_dir,
    winslash = "/",
    mustWork = TRUE
  ),
  "\n\n",
  sep = ""
)

cat("Most useful files to send back:\n")
cat("V251_model_comparison.csv\n")
cat("V251_run_status.csv\n")
cat("V251_preliminary_ranking.csv\n")
cat("V251_adjacent_IC_change.csv\n")
cat("V251_screen_report.txt\n\n")

cat("For individual G diagnostics:\n")
cat("V251_task_summary_GX.csv\n")
cat("V251_class_quality_GX.csv\n")
cat("V251_EM_runs_GX.csv\n")
cat("V251_EM_progress_GX.log\n")
cat("V251_stage_GX.log\n\n")

cat("Reminder: scientific settings remain the V251 broad screen with tol=10.\n")
cat("V258 changes execution/audit handling only: checkpoint/resume + dynamic CPU scheduling + RAM-aware launch + automatic completed-G skip + legacy status refresh.\n")
cat("All V251-compatible per-G outputs are written into the existing GBMTM_V251_G5_G20_screen folder.\n")
cat("Final shortlisted models must later use strict multi-start refinement.\n")
cat("============================================================\n")
