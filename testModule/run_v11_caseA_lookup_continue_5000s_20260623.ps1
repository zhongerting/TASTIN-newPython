$ErrorActionPreference = "Stop"

$RepoRoot = Split-Path -Parent $PSScriptRoot
Set-Location $RepoRoot

$Python = "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe"
$OutputDir = "testModule\v11_lookup_continue_5000s_from22240_20260623"
$CasePrefix = "v11_lookup_continue_5000s_from22240_20260623"

$env:THERMOCALC_PYD_DIR = (Resolve-Path "ThermoCalc\build_cp312\Release").Path
$env:THERMOCALC_ENABLE_LOOKUP = "1"
$env:THERMOCALC_LOOKUP_DB = (Resolve-Path "ThermoCalc\emission_runtime_db_v2").Path
$env:THERMOCALC_LOOKUP_REGIONS = "core"

New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null
"wrapper start $(Get-Date -Format o)" | Out-File -Encoding utf8 "$OutputDir\wrapper.log"
"THERMOCALC_PYD_DIR=$env:THERMOCALC_PYD_DIR" | Out-File -Encoding utf8 -Append "$OutputDir\wrapper.log"
"THERMOCALC_LOOKUP_DB=$env:THERMOCALC_LOOKUP_DB" | Out-File -Encoding utf8 -Append "$OutputDir\wrapper.log"
"THERMOCALC_LOOKUP_REGIONS=$env:THERMOCALC_LOOKUP_REGIONS" | Out-File -Encoding utf8 -Append "$OutputDir\wrapper.log"

try {
  & $Python testModule\run_v11_caseA_closed_loop.py `
    --restart-in testModule\v11_caseA_closed_loop_ref_ring024_outer030_20260616\v11_caseA_closed_loop_ref_ring024_outer030_20260616_latest_restart.npz `
    --duration 5000 `
    --record-interval 50 `
    --restart-interval 50 `
    --max-dt 0.05 `
    --ring-emissivity 0.24 `
    --outer-header-emissivity 0.30 `
    --pump-total-head-pa 6483.548313292204 `
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
