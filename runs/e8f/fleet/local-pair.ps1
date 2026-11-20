# local-pair.ps1 : run ONE frozen-orbital pair that the fleet could not finish (killed on memory at about
# 15 GB RAM with spilling off), on a Windows machine with room to spill the whole table to disk.
# It measures what the fleet never could: how long a pair takes with the table spilled, and the peak
# memory the process reaches.
#
# From the repository root (E:\Projects\phys), after `git fetch origin claude/e8f-orchestrator` and
# `git checkout claude/e8f-orchestrator`:
#   cargo build --release --bin phys-overlap
#   powershell -ExecutionPolicy Bypass -File runs\e8f\fleet\local-pair.ps1
# Optional: -Name butane -Part 2 -Pair 1 -SpillDir E:\phys-spill -SpillGB 300
#
# Writes runs\e8f\species\<Name>\local-<Part>-<Pair>.log (the engine's output), local-<Part>-<Pair>.txt (the
# pair's result line) and local-<Part>-<Pair>.mem (elapsed time and peak working set, sampled every 30 s).
param(
  [string]$Name = "butane",
  [int]$Part = 2,
  [int]$Pair = 1,
  [string]$SpillDir = "E:\phys-spill",
  [double]$SpillGB = 300
)
$ErrorActionPreference = "Stop"
$root = (git rev-parse --show-toplevel).Trim()
$exe = Join-Path $root "target\release\phys-overlap.exe"
if (-not (Test-Path $exe)) { throw "build it first: cargo build --release --bin phys-overlap" }
$dir = Join-Path $root "runs\e8f\species\$Name"
Set-Location $dir

$sel = "sel-part-$Part-$Pair.txt"
if (-not (Test-Path $sel)) {
  # the fleet's 3-pair parts: pair J of part K is line (K-1)*3+J of the selection
  $line = ($Part - 1) * 3 + $Pair
  (Get-Content "pairs-$Name-sel.txt" | Where-Object { $_ -ne "" })[$line - 1] | Set-Content $sel
}

New-Item -ItemType Directory -Force -Path $SpillDir | Out-Null
Get-ChildItem $SpillDir -Force | Remove-Item -Recurse -Force
$free = (Get-PSDrive ($SpillDir.Substring(0,1))).Free / 1e9
Write-Host ("spill to {0}: {1:N0} GB free, cap {2} GB" -f $SpillDir, $free, $SpillGB)
if ($free -lt 110) { Write-Warning "the largest table seen for butane is 103 GB; this disk may not hold it" }

$env:PHYS_SPILL_DIR = $SpillDir
$env:PHYS_SPILL_MAX_GB = "$SpillGB"
$tag = "local-$Part-$Pair"
$out = "$tag.txt"; $log = "$tag.log"; $mem = "$tag.mem"
if (Test-Path $out) { Remove-Item $out }

$start = Get-Date
$p = Start-Process -FilePath $exe -ArgumentList @($Name, $sel, "--frozen", "--out", $out) `
  -RedirectStandardOutput $log -RedirectStandardError "$tag.err" -NoNewWindow -PassThru
$peak = 0
while (-not $p.HasExited) {
  Start-Sleep -Seconds 30
  try { $p.Refresh(); if ($p.PeakWorkingSet64 -gt $peak) { $peak = $p.PeakWorkingSet64 } } catch {}
  $h = ((Get-Date) - $start).TotalHours
  "{0:N2} h  peak {1:N1} GB  spill {2:N1} GB" -f $h, ($peak / 1e9), `
    ((Get-ChildItem $SpillDir -Recurse -File -ErrorAction SilentlyContinue | Measure-Object Length -Sum).Sum / 1e9) |
    Tee-Object -FilePath $mem -Append | Out-Null
}
$h = ((Get-Date) - $start).TotalHours
$lines = if (Test-Path $out) { (Get-Content $out | Where-Object { $_ -ne "" }).Count } else { 0 }
$summary = "exit {0}; {1:N2} h; peak working set {2:N1} GB; result lines {3}" -f $p.ExitCode, $h, ($peak / 1e9), $lines
$summary | Tee-Object -FilePath $mem -Append
Get-Content "$tag.err" -ErrorAction SilentlyContinue | Add-Content $log
Write-Host "paste the last lines of $mem and $log to the orchestrator, or commit runs\e8f\species\$Name\$tag.* to a branch claude/e8f-local and push it."
