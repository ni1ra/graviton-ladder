$root = Split-Path -Parent $PSScriptRoot
$n = (Get-ChildItem "$root\logs" -Filter *.log | Measure-Object).Count + 1
$log = "{0}\logs\{1:D3}-run.log" -f $root, $n
& C:\Python314\python.exe "$PSScriptRoot\second.py" 2>&1 | Tee-Object -FilePath $log
"--- flint check ---" | Tee-Object -FilePath $log -Append
& record\wave73\even-nonzero\venv\Scripts\python.exe "$PSScriptRoot\check_flint.py" 2>&1 | Tee-Object -FilePath $log -Append
