# wait.ps1 <logfile> <pattern> : poll up to 165 s, one status line. Stops by itself (bounded loop).
param([string]$log, [string]$pat)
Set-Location $PSScriptRoot
$t = 0; $minfree = 99
while ($t -lt 33) {
    $f = (Get-CimInstance Win32_OperatingSystem).FreePhysicalMemory / 1MB
    if ($f -lt $minfree) { $minfree = $f }
    if ((Test-Path $log) -and (Select-String -Path $log -Pattern $pat -Quiet)) { break }
    Start-Sleep 5; $t++
}
$hit = (Test-Path $log) -and (Select-String -Path $log -Pattern $pat -Quiet)
"{0} hit={1} rows={2} minfree={3:N1}GB {4}" -f (Get-Date -Format HH:mm), $hit, (Get-ChildItem results\raw\E_*.json).Count, $minfree, ((Get-ChildItem logs\1*,logs\2* -ErrorAction SilentlyContinue | % { $_.Name.Substring(0,2) + ':' + (Get-Content $_.FullName | Measure-Object -Line).Lines }) -join ' ')
