param(
  [Parameter(Mandatory=$true)][ValidateSet("topology","affinity","memory","process_memory")][string]$Mode,
  [int]$TargetPid = 0,
  [int]$LogicalIndex = -1,
  [string]$PriorityClass = "AboveNormal",
  [string]$PidList = ""
)
$ErrorActionPreference = "Stop"

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

if ($Mode -eq "process_memory") {
  $ids = @()
  if ($PidList) {
    $ids = $PidList.Split(",") | ForEach-Object { [int]$_.Trim() }
  }
  $procs = @()
  foreach ($id in $ids) {
    $p = Get-Process -Id $id -ErrorAction SilentlyContinue
    if ($null -ne $p) { $procs += $p }
  }
  $ws = ($procs | Measure-Object -Property WorkingSet64 -Sum).Sum
  $pm = ($procs | Measure-Object -Property PrivateMemorySize64 -Sum).Sum
  if ($null -eq $ws) { $ws = 0 }
  if ($null -eq $pm) { $pm = 0 }
  [PSCustomObject]@{
    ProcessCount = $procs.Count
    WorkingSetMB = [math]::Round($ws / 1MB, 3)
    PrivateMemoryMB = [math]::Round($pm / 1MB, 3)
  } | ConvertTo-Csv -NoTypeInformation
  exit 0
}

if ($Mode -eq "affinity") {
  if ($TargetPid -le 0) { throw "TargetPid must be positive." }
  if ($LogicalIndex -lt 0 -or $LogicalIndex -gt 62) { throw "LogicalIndex must be 0-62." }
  $p = Get-Process -Id $TargetPid
  $mask = ([Int64]1 -shl $LogicalIndex)
  $p.ProcessorAffinity = [IntPtr]$mask
  if ($PriorityClass -and $PriorityClass -ne "Normal") {
    $p.PriorityClass = $PriorityClass
  }
  [PSCustomObject]@{
    PID = $TargetPid
    LogicalIndex = $LogicalIndex
    PriorityClass = $p.PriorityClass.ToString()
    Status = "OK"
  } | ConvertTo-Csv -NoTypeInformation
  exit 0
}

$src = @"
using System;
using System.Runtime.InteropServices;
public static class V295CpuSetNative {
  [DllImport("kernel32.dll", SetLastError=true)]
  public static extern bool GetSystemCpuSetInformation(
    IntPtr Information, uint BufferLength, out uint ReturnedLength,
    IntPtr Process, uint Flags);
}
"@
Add-Type -TypeDefinition $src
$needed = [uint32]0
[void][V295CpuSetNative]::GetSystemCpuSetInformation([IntPtr]::Zero, 0, [ref]$needed, [IntPtr]::Zero, 0)
if ($needed -le 0) { throw "GetSystemCpuSetInformation returned zero bytes." }
$buf = [Runtime.InteropServices.Marshal]::AllocHGlobal([int]$needed)
try {
  $returned = [uint32]0
  $ok = [V295CpuSetNative]::GetSystemCpuSetInformation($buf, $needed, [ref]$returned, [IntPtr]::Zero, 0)
  if (-not $ok) { throw "GetSystemCpuSetInformation failed." }
  $cpuName = (Get-CimInstance Win32_Processor | Select-Object -First 1 -ExpandProperty Name)
  $offset = 0
  $rows = @()
  while ($offset -lt $returned) {
    $ptr = [IntPtr]::Add($buf, $offset)
    $size = [Runtime.InteropServices.Marshal]::ReadInt32($ptr, 0)
    $type = [Runtime.InteropServices.Marshal]::ReadInt32($ptr, 4)
    if ($size -le 0) { throw "Invalid CPU-set structure size." }
    if ($type -eq 0) {
      $rows += [PSCustomObject]@{
        ProcessorName = $cpuName
        Group = [int][Runtime.InteropServices.Marshal]::ReadInt16($ptr, 12)
        LogicalProcessorIndex = [int][Runtime.InteropServices.Marshal]::ReadByte($ptr, 14)
        CoreIndex = [int][Runtime.InteropServices.Marshal]::ReadByte($ptr, 15)
        LastLevelCacheIndex = [int][Runtime.InteropServices.Marshal]::ReadByte($ptr, 16)
        EfficiencyClass = [int][Runtime.InteropServices.Marshal]::ReadByte($ptr, 18)
        SchedulingClass = [int][Runtime.InteropServices.Marshal]::ReadByte($ptr, 20)
      }
    }
    $offset += $size
  }
  $rows | Sort-Object Group,LogicalProcessorIndex | ConvertTo-Csv -NoTypeInformation
} finally {
  [Runtime.InteropServices.Marshal]::FreeHGlobal($buf)
}
