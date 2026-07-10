$ErrorActionPreference = "Stop"

$RepoRoot = Split-Path -Parent $PSScriptRoot
Set-Location $RepoRoot

$Python = "E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe"
$OutputDir = "testModule\v11_lookup_pcs_0p02_5torr_12000s_from23500_20260623"
$CasePrefix = "v11_lookup_pcs_0p02_5torr_12000s_from23500_20260623"
$RestartIn = "testModule\v11_lookup_pcs_0p02_5torr_long_from22450_20260623\v11_lookup_pcs_0p02_5torr_long_from22450_20260623_latest_restart.npz"

$env:THERMOCALC_PYD_DIR = (Resolve-Path "ThermoCalc\build_cp312\Release").Path
$env:THERMOCALC_ENABLE_LOOKUP = "1"
$env:THERMOCALC_LOOKUP_DB = (Resolve-Path "ThermoCalc\emission_runtime_db_v2\pcs_0p02_5torr").Path
$env:THERMOCALC_LOOKUP_REGIONS = "core,startup,high_power,accident"

New-Item -ItemType Directory -Force -Path $OutputDir | Out-Null
"wrapper start $(Get-Date -Format o)" | Out-File -Encoding utf8 "$OutputDir\wrapper.log"
"THERMOCALC_PYD_DIR=$env:THERMOCALC_PYD_DIR" | Out-File -Encoding utf8 -Append "$OutputDir\wrapper.log"
"THERMOCALC_ENABLE_LOOKUP=$env:THERMOCALC_ENABLE_LOOKUP" | Out-File -Encoding utf8 -Append "$OutputDir\wrapper.log"
"THERMOCALC_LOOKUP_DB=$env:THERMOCALC_LOOKUP_DB" | Out-File -Encoding utf8 -Append "$OutputDir\wrapper.log"
"THERMOCALC_LOOKUP_REGIONS=$env:THERMOCALC_LOOKUP_REGIONS" | Out-File -Encoding utf8 -Append "$OutputDir\wrapper.log"
"restart_in=$RestartIn" | Out-File -Encoding utf8 -Append "$OutputDir\wrapper.log"

& $Python testModule\probe_v11_lookup_hit_rate_20260623.py `
  --restart-in $RestartIn `
  --output-json "$OutputDir\lookup_preflight.json" `
  --lookup-db $env:THERMOCALC_LOOKUP_DB `
  --lookup-regions $env:THERMOCALC_LOOKUP_REGIONS `
  --pump-total-head-pa 6483.548313292204 `
  --ring-emissivity 0.24 `
  --outer-header-emissivity 0.30 `
  1> "$OutputDir\lookup_preflight.out" `
  2> "$OutputDir\lookup_preflight.err"

try {
  & $Python testModule\run_v11_caseA_closed_loop.py `
    --restart-in $RestartIn `
    --duration 12000 `
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
