# V14 Fixed-I Low-Power Tuning Execution Plan

> **For agentic workers:** use subagent review for implementation changes; numerical candidates run as isolated Python processes. The main thread owns all parameter decisions.

**Goal:** Find fixed `Q_low/W_low` and a smooth descent that reaches `Pe/Pe0=38%-42%`, preserves the thermal system, and remains valid for a frozen full external-heat orbit.

**Architecture:** One runner produces isolated, restartable candidates from the same immutable checkpoint. Endpoint tuning uses a standard provisional ramp and frozen hold; trajectory tuning starts only after endpoint selection. Parallelism is across processes and candidates, never inside one ThermoCalc process.

**Runtime:** `E:\Users\HC Zhao\anaconda3\envs\tastin-python\python.exe`

## Global Constraints

- Initial checkpoint: `checkpoint_t019265s.npz`, SHA256 `3EC31DDB1538DE299C346317BBA7F644A6995327A3BC3493FB07D0E26C65D88D`.
- Main TEC is series `fixed_i` with runtime-refreshed `I0`; fallback/open-circuit and non-convergence invalidate a candidate.
- Thermal power is primary; flow is auxiliary. Both change continuously and have independent endpoints.
- At freeze: `Pe/Pe0=0.38..0.42`. Frozen orbit: every record `0.35..0.45`, orbit mean `0.38..0.42`.
- `Tout <= Tout0 + 0.5 K`; fuel 2700 K, collector 1500 K, emitter 3000 K, coolant 1058 K, moderator 930 K, reflector 1000 K.
- Prefer the endpoint whose heat-pipe, collector-ring and radiator temperatures are closest to the same-phase baseline; break ties by smaller flow reduction.
- `I0/Q_low/W_low` never change after freeze.
- Every candidate starts from the same original checkpoint and writes to a unique directory.
- Numerical candidates use independent processes with `MKL_NUM_THREADS=OMP_NUM_THREADS=OPENBLAS_NUM_THREADS=NUMEXPR_NUM_THREADS=1`.

## Task 1: Qualify the runner for unattended runs

**Files:** `run_v14_low_power_fixed_i.py`, `test_low_power_fixed_i.py`

- [ ] Add failing tests for exact endpoint application, staged recording, hard trips, incremental history, checkpoint rotation and non-zero failure exit.
- [ ] Verify the tests fail for the missing behavior.
- [ ] Reuse existing V14/LOCA temperature collectors and restart helpers where possible; implement only the missing runner behavior.
- [ ] Apply setpoints at the step end, explicitly apply exact `Q_low/W_low` before freeze/save, then verify the stored values.
- [ ] Write each record immediately and flush; save periodic restart through a temporary file and `os.replace`.
- [ ] Run focused tests and a 0.1 s integration smoke.

## Task 2: Add selection diagnostics

**Files:** `run_v14_low_power_fixed_i.py`, focused test file, `DEBUG_LOG.md`

- [ ] Add a failing synthetic diagnostic test.
- [ ] Record heat-pipe evaporator/adiabatic/condenser min/mean/max temperature, axial conductivity min/mean/max, axial temperature difference and rejection.
- [ ] Record collector-ring and radiator temperature/rejection metrics plus the existing electrical, outlet and limit metrics.
- [ ] Record same-phase baseline deltas used to rank `W_low`.
- [ ] Verify finite values on the checkpoint smoke and document field definitions.

## Task 3: Establish safe process concurrency

- [ ] Launch two identical 20 s read-only candidates with unique output/log directories and staggered starts.
- [ ] Measure per-process wall time, aggregate throughput and available memory.
- [ ] Accept two-way concurrency only if throughput is at least 1.6x and per-process slowdown is at most 25% with at least 12 GB memory free.
- [ ] Test three-way concurrency only if two-way passes; cap production at three unless a four-way probe also meets the same gate.
- [ ] Record the selected concurrency in `DEBUG_LOG.md`.

## Task 4: Local response batch

- [ ] Run in parallel from the original checkpoint: `(q,w)=(0.95,1.00)`, `(1.00,0.95)`, `(0.95,0.95)`.
- [ ] For each: identical 20 s fixed-I baseline, 100 s provisional smootherstep, 200 s frozen hold.
- [ ] Reject hard failures; compare three consecutive hold windows to detect continuing drift.
- [ ] Estimate local `Pe/Tout/heat-pipe-temperature` response to `Q`, `W`, and their interaction.

## Task 5: Endpoint search

- [ ] Generate only three candidates per batch: predicted center and two bracketing points.
- [ ] Limit one-batch movement to 5 percentage points in `Q/Q0` and 3 percentage points in `W/W0`.
- [ ] Use the same provisional curve and duration for every candidate in a batch.
- [ ] A candidate is not a hit until exact `Q/W` are frozen for 300-600 s and drift projection remains inside `Pe/Pe0=0.35..0.45`.
- [ ] Rank valid hits by same-phase heat-pipe/collector/radiator temperature error, then by smaller flow reduction.
- [ ] Freeze the selected `Q_low/W_low` in the log; do not tune them during trajectory search.

## Task 6: Smooth trajectory search

- [ ] Run three 1000 s candidates in parallel with fixed endpoints: cubic synchronous, quintic synchronous, and quintic with flow delayed by `0.1T`.
- [ ] Reject hard failures, endpoint values outside `0.38..0.42`, or projected frozen drift outside `0.35..0.45`.
- [ ] Rank survivors lexicographically by endpoint compliance, final-100 s drift, `TV(Pe)/abs(Pe_end-Pe_start)`, peak 60 s electrical slope, and outlet-temperature slope.
- [ ] If the best path is still moving toward the target, rerun only the best one or two from the original checkpoint at 1500 s, then 2000 s if needed.
- [ ] If a settled endpoint misses the target, return to Task 5 instead of extending time.

## Task 7: Final qualification

- [ ] Recompute the selected endpoint and trajectory once in a single process from the original checkpoint.
- [ ] Freeze for 300-600 s, then run at least one complete `5668.144369 s` external-heat period without changing `I0/Q_low/W_low`.
- [ ] Verify instantaneous/mean power bands, all thermal limits, convergence and same-phase start/end drift.
- [ ] If key temperature drift exceeds 2 K or electrical drift exceeds `0.01 Pe0`, run one additional orbit.
- [ ] Write final parameters, trajectory, hashes, histories, failures and acceptance evidence to `DEBUG_LOG.md` and the case README.
