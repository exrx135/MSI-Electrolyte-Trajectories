param(
  [Parameter(Mandatory=$true)][ValidateSet("topology","affinity","memory","process")][string]$Mode,
  [int]$TargetPid = 0,
  [int]$LogicalIndex = -1,
  [string]$PriorityClass = "AboveNormal"
)

$ErrorActionPreference = "Stop"

if ($Mode -eq "process") {
  if ($TargetPid -le 0) { throw "TargetPid must be positive." }
  $p = Get-Process -Id $TargetPid -ErrorAction SilentlyContinue
  if ($null -eq $p) {
    [PSCustomObject]@{
      PID = $TargetPid
      ProcessName = ""
      ExecutablePath = ""
      StartTimeLocal = ""
      CommandLine = ""
      Status = "NOT_FOUND"
    } | ConvertTo-Csv -NoTypeInformation
  } else {
    $exePath = ""
    $startLocal = ""
    $commandLine = ""
    try {
      if ($null -ne $p.Path) { $exePath = [string]$p.Path }
    } catch {}
    try {
      $startLocal = $p.StartTime.ToString("yyyy-MM-dd HH:mm:ss.fff")
    } catch {}
    try {
      $cim = Get-CimInstance Win32_Process -Filter ("ProcessId = " + $TargetPid) -ErrorAction SilentlyContinue
      if ($null -ne $cim) {
        if ($cim.ExecutablePath) { $exePath = [string]$cim.ExecutablePath }
        if ($cim.CommandLine) { $commandLine = [string]$cim.CommandLine }
      }
    } catch {}
    [PSCustomObject]@{
      PID = $TargetPid
      ProcessName = $p.ProcessName
      ExecutablePath = $exePath
      StartTimeLocal = $startLocal
      CommandLine = $commandLine
      Status = "FOUND"
    } | ConvertTo-Csv -NoTypeInformation
  }
  exit 0
}

if ($Mode -eq "memory") {
  $os = Get-CimInstance Win32_OperatingSystem
  [PSCustomObject]@{
    TotalPhysicalMB = [math]::Round($os.TotalVisibleMemorySize / 1024.0, 3)
    AvailablePhysicalMB = [math]::Round($os.FreePhysicalMemory / 1024.0, 3)
    CommitLimitMB = [math]::Round($os.TotalVirtualMemorySize / 1024.0, 3)
    CommitUsedMB = [math]::Round(($os.TotalVirtualMemorySize - $os.FreeVirtualMemory) / 1024.0, 3)
  } | ConvertTo-Csv -NoTypeInformation
  exit 0
}

if ($Mode -eq "affinity") {
  if ($TargetPid -le 0) { throw "TargetPid must be positive." }
  if ($LogicalIndex -lt 0 -or $LogicalIndex -gt 62) { throw "LogicalIndex must be between 0 and 62." }
  $p = Get-Process -Id $TargetPid
  $mask = ([Int64]1 -shl $LogicalIndex)
  $p.ProcessorAffinity = [IntPtr]$mask
  if ($PriorityClass -and $PriorityClass -ne "Normal") {
    $p.PriorityClass = $PriorityClass
  }
  [PSCustomObject]@{
    PID = $TargetPid
    LogicalIndex = $LogicalIndex
    AffinityMask = $mask
    PriorityClass = $p.PriorityClass.ToString()
    Status = "OK"
  } | ConvertTo-Csv -NoTypeInformation
  exit 0
}

$src = @"
using System;
using System.Runtime.InteropServices;
public static class V258CpuSetNative {
  [DllImport("kernel32.dll", SetLastError=true)]
  public static extern bool GetSystemCpuSetInformation(
    IntPtr Information,
    uint BufferLength,
    out uint ReturnedLength,
    IntPtr Process,
    uint Flags);
}
"@

Add-Type -TypeDefinition $src
$needed = [uint32]0
[void][V258CpuSetNative]::GetSystemCpuSetInformation([IntPtr]::Zero, 0, [ref]$needed, [IntPtr]::Zero, 0)
if ($needed -le 0) { throw "GetSystemCpuSetInformation returned zero required bytes." }
$buf = [Runtime.InteropServices.Marshal]::AllocHGlobal([int]$needed)
try {
  $returned = [uint32]0
  $ok = [V258CpuSetNative]::GetSystemCpuSetInformation($buf, $needed, [ref]$returned, [IntPtr]::Zero, 0)
  if (-not $ok) { throw "GetSystemCpuSetInformation failed." }
  $cpuName = (Get-CimInstance Win32_Processor | Select-Object -First 1 -ExpandProperty Name)
  $offset = 0
  $rows = @()
  while ($offset -lt $returned) {
    $ptr = [IntPtr]::Add($buf, $offset)
    $size = [Runtime.InteropServices.Marshal]::ReadInt32($ptr, 0)
    $type = [Runtime.InteropServices.Marshal]::ReadInt32($ptr, 4)
    if ($size -le 0) { throw "Invalid CPU set structure size." }
    if ($type -eq 0) {
      $id = [Runtime.InteropServices.Marshal]::ReadInt32($ptr, 8)
      $group = [Runtime.InteropServices.Marshal]::ReadInt16($ptr, 12)
      $logical = [Runtime.InteropServices.Marshal]::ReadByte($ptr, 14)
      $core = [Runtime.InteropServices.Marshal]::ReadByte($ptr, 15)
      $llc = [Runtime.InteropServices.Marshal]::ReadByte($ptr, 16)
      $numa = [Runtime.InteropServices.Marshal]::ReadByte($ptr, 17)
      $eff = [Runtime.InteropServices.Marshal]::ReadByte($ptr, 18)
      $flags = [Runtime.InteropServices.Marshal]::ReadByte($ptr, 19)
      $sched = [Runtime.InteropServices.Marshal]::ReadByte($ptr, 20)
      $rows += [PSCustomObject]@{
        ProcessorName = $cpuName
        CpuSetId = [uint32]$id
        Group = [int]$group
        LogicalProcessorIndex = [int]$logical
        CoreIndex = [int]$core
        LastLevelCacheIndex = [int]$llc
        NumaNodeIndex = [int]$numa
        EfficiencyClass = [int]$eff
        SchedulingClass = [int]$sched
        AllFlags = [int]$flags
      }
    }
    $offset += $size
  }
  $rows | Sort-Object Group,LogicalProcessorIndex | ConvertTo-Csv -NoTypeInformation
} finally {
  [Runtime.InteropServices.Marshal]::FreeHGlobal($buf)
}
