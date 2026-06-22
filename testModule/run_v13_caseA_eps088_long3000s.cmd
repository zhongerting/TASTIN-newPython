@echo off
set OUTDIR=testModule\v13_caseA_closed_loop_eps088_long3000s
if not exist "%OUTDIR%" mkdir "%OUTDIR%"
"E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe" testModule\run_v13_caseA_closed_loop.py ^
  --restart-in testModule\v13_caseA_closed_loop_eps088_100s_probe\v13_caseA_closed_loop_eps088_100s_probe_latest_restart.npz ^
  --tube-emissivity 0.88 ^
  --fin-emissivity 0.88 ^
  --duration 3000 ^
  --record-interval 300 ^
  --restart-interval 300 ^
  --max-dt 0.5 ^
  --output-dir "%OUTDIR%" ^
  --case-prefix v13_caseA_closed_loop_eps088_long3000s ^
  > "%OUTDIR%\run.out" 2> "%OUTDIR%\run.err"
