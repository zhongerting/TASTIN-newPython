$ErrorActionPreference = "Stop"

$RepoRoot = Split-Path -Parent $PSScriptRoot
Set-Location $RepoRoot

$Python = "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe"
$OutputDir = "testModule\v11_caseA_closed_loop_ref_overnight_20260616"
$CasePrefix = "v11_caseA_closed_loop_ref_overnight_20260616"

New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null
"wrapper start $(Get-Date -Format o)" | Out-File -Encoding utf8 "$OutputDir\wrapper.log"

try {
  & $Python testModule\run_v11_caseA_closed_loop.py `
    --restart-in testModule\v11_caseA_closed_loop_ref_smoke_10s\v11_caseA_closed_loop_ref_smoke_10s_latest_restart.npz `
    --duration 50000 `
    --record-interval 50 `
    --restart-interval 50 `
    --max-dt 0.05 `
    --enable-pump-head-control `
    --pump-control-interval 50 `
    --output-dir $OutputDir `
    --case-prefix $CasePrefix `
    1> "$OutputDir\run.out" `
    2> "$OutputDir\run.err"
} catch {
  $_ | Out-File -Encoding utf8 "$OutputDir\wrapper.err"
  throw
}

"wrapper done $(Get-Date -Format o)" | Out-File -Encoding utf8 -Append "$OutputDir\wrapper.log"
