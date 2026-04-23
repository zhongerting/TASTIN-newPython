"""
Standalone verification case for a long-lasting small dU/dt tail.

The purpose is to isolate the energy-accounting question from the full v4_2
coolant-loop model.  The case represents a finite thermal mass connected to a
steady 96 kW heat-transfer path:

    C * dT/dt = Q_loop - Q_rej
    Q_rej     = Q_loop + G * (T - T_ss)

Therefore:

    dU/dt   = C * dT/dt = -G * (T - T_ss)
    residual = Q_loop - Q_rej - dU/dt = 0

With a large C/G time constant, dU/dt can remain close to -0.11 kW for
thousands of seconds while releasing only a few hundred kJ in total.
"""

import csv
import os

import matplotlib.pyplot as plt
import numpy as np


Q_LOOP_KW = 96.37
T_SS = 765.0
C_TOTAL = 10.0e6          # J/K, order of the finite-domain heat capacity
TAU = 60_000.0            # s, long tail time constant
G = C_TOTAL / TAU         # W/K
T_REF = 1000.0            # s
DUDT_REF_KW = -0.112      # kW, prints as -0.11 kW
T_END = 5000.0
DT = 1.0

OUT_DIR = os.path.dirname(os.path.abspath(__file__))
CSV_PATH = os.path.join(OUT_DIR, "verify_dudt_long_tail_case.csv")
PLOT_PATH = os.path.join(OUT_DIR, "verify_dudt_long_tail_case.png")


def analytic_state(t):
    """Return theta, q_rej_kw, dudt_kw, residual_kw at time t."""
    theta_ref = abs(DUDT_REF_KW) * 1000.0 / G
    theta = theta_ref * np.exp(-(t - T_REF) / TAU)
    q_rej_kw = Q_LOOP_KW + (G * theta) / 1000.0
    dudt_kw = -(G * theta) / 1000.0
    residual_kw = Q_LOOP_KW - q_rej_kw - dudt_kw
    return theta, q_rej_kw, dudt_kw, residual_kw


def integrate_storage(t1, t2):
    """Analytic released storage energy between t1 and t2, positive in J."""
    theta1 = analytic_state(t1)[0]
    theta2 = analytic_state(t2)[0]
    return C_TOTAL * (theta1 - theta2)


def main():
    times = np.arange(0.0, T_END + DT, DT)
    rows = []

    for t in times:
        theta, q_rej_kw, dudt_kw, residual_kw = analytic_state(t)
        rows.append({
            "t": t,
            "T_storage": T_SS + theta,
            "theta_to_steady": theta,
            "Q_loop_kW": Q_LOOP_KW,
            "Q_rej_kW": q_rej_kw,
            "dUdt_kW": dudt_kw,
            "residual_kW": residual_kw,
            "U_excess_MJ": C_TOTAL * theta / 1.0e6,
        })

    with open(CSV_PATH, "w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0].keys()))
        writer.writeheader()
        writer.writerows(rows)

    fig, axes = plt.subplots(3, 1, figsize=(10, 9), sharex=True)
    t_arr = np.array([r["t"] for r in rows])
    dudt_arr = np.array([r["dUdt_kW"] for r in rows])
    q_loop_arr = np.array([r["Q_loop_kW"] for r in rows])
    q_rej_arr = np.array([r["Q_rej_kW"] for r in rows])
    theta_arr = np.array([r["theta_to_steady"] for r in rows])
    u_arr = np.array([r["U_excess_MJ"] for r in rows])

    axes[0].plot(t_arr, q_loop_arr, label="Q_loop", lw=2)
    axes[0].plot(t_arr, q_rej_arr, label="Q_rej = Q_loop - dUdt", lw=2)
    axes[0].set_ylabel("Heat rate [kW]")
    axes[0].grid(True, linestyle="--", alpha=0.45)
    axes[0].legend(loc="best")

    axes[1].plot(t_arr, dudt_arr, color="tab:red", lw=2, label="dUdt")
    axes[1].axhline(-0.11, color="gray", linestyle=":", label="-0.11 kW")
    axes[1].set_ylabel("dU/dt [kW]")
    axes[1].grid(True, linestyle="--", alpha=0.45)
    axes[1].legend(loc="best")

    axes[2].plot(t_arr, theta_arr, color="tab:green", lw=2, label="T - T_ss")
    axes[2].plot(t_arr, u_arr, color="tab:purple", lw=2, label="U_excess [MJ]")
    axes[2].set_xlabel("Time [s]")
    axes[2].set_ylabel("Tail state")
    axes[2].grid(True, linestyle="--", alpha=0.45)
    axes[2].legend(loc="best")

    plt.tight_layout()
    plt.savefig(PLOT_PATH, dpi=130)
    plt.close()

    print("=" * 72)
    print("Standalone dU/dt long-tail verification")
    print("=" * 72)
    print(f"Q_loop             = {Q_LOOP_KW:.3f} kW")
    print(f"C_total            = {C_TOTAL/1e6:.3f} MJ/K")
    print(f"tau = C/G          = {TAU:.1f} s")
    print(f"G                  = {G:.3f} W/K")
    print(f"dUdt at {T_REF:g}s       = {analytic_state(T_REF)[2]:.6f} kW")
    print()

    for t in [1000.0, 2000.0, 3000.0, 4000.0, 4274.87]:
        theta, q_rej_kw, dudt_kw, residual_kw = analytic_state(t)
        print(
            f"t={t:7.2f}s | T_storage={T_SS + theta:8.5f}K "
            f"| theta={theta:8.5f}K | Q_rej={q_rej_kw:8.4f}kW "
            f"| dUdt={dudt_kw:8.5f}kW | printed={dudt_kw:5.2f}kW "
            f"| residual={residual_kw: .3e}kW"
        )

    released_1000_4000 = integrate_storage(1000.0, 4000.0)
    main_heat_1000_4000 = Q_LOOP_KW * 1000.0 * (4000.0 - 1000.0)
    theta_1000 = analytic_state(1000.0)[0]
    theta_4000 = analytic_state(4000.0)[0]

    print()
    print("[1000s, 4000s] integral check")
    print(f"released storage   = {released_1000_4000/1e6:.6f} MJ")
    print(f"main heat transfer = {main_heat_1000_4000/1e6:.3f} MJ")
    print(f"storage/main ratio = {released_1000_4000/main_heat_1000_4000*100:.4f} %")
    print(f"theta drop         = {theta_1000 - theta_4000:.6f} K")
    print()
    print(f"CSV saved:  {CSV_PATH}")
    print(f"Plot saved: {PLOT_PATH}")


if __name__ == "__main__":
    main()
