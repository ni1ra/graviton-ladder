<#
.SYNOPSIS
  Refresh the Pages site: public page from origin/main, nightly preview from origin/nightly.

.DESCRIPTION
  Fetches origin, then runs site/publish.py --push. Appends everything to the log file
  (outside the repository) and exits non-zero if the fetch or the publish failed.
  Meant for the Windows Task Scheduler (see site/README.md); also fine to run by hand.
  -DryRun prepares the commit without pushing. -Nightly picks another ref for the preview.
#>
param(
  [string]$Stable  = 'origin/main',
  [string]$Nightly = 'origin/nightly',
  [string]$Log     = 'C:\Users\nira\work\artifacts\graviton-site\nightly-log.txt',
  [switch]$DryRun
)

$ErrorActionPreference = 'Continue'   # native commands write progress to stderr; the exit codes decide
$repo = Split-Path -Parent $PSScriptRoot
$python = if ($env:GRAVITON_PYTHON) { $env:GRAVITON_PYTHON }
          elseif (Test-Path 'C:\Python314\python.exe') { 'C:\Python314\python.exe' }
          else { 'python' }
$env:PYTHONUTF8 = '1'

New-Item -ItemType Directory -Force -Path (Split-Path -Parent $Log) | Out-Null

function Write-Log([string]$text) {
  $line = "$((Get-Date).ToUniversalTime().ToString('u')) $text"
  Add-Content -Path $Log -Value $line -Encoding utf8
  Write-Host $line
}

Write-Log "--- nightly publish start (stable=$Stable nightly=$Nightly dryrun=$($DryRun.IsPresent))"

& git -C $repo fetch --quiet origin 2>&1 | ForEach-Object { Write-Log "$_" }
if ($LASTEXITCODE -ne 0) {
  Write-Log "FAILED: git fetch exited with $LASTEXITCODE"
  exit 1
}

# --allow-dirty: the refs are named explicitly and read from git, so local edits cannot leak in.
$pyArgs = @("$repo\site\publish.py", '--stable', $Stable, '--nightly', $Nightly, '--allow-dirty')
if (-not $DryRun) { $pyArgs += '--push' }
& $python @pyArgs 2>&1 | ForEach-Object { Write-Log "$_" }
$code = $LASTEXITCODE
if ($code -ne 0) {
  Write-Log "FAILED: publish.py exited with $code"
  exit $code
}

Write-Log "--- nightly publish done"
exit 0
