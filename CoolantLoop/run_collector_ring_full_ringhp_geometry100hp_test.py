import argparse
import argparse
import math
import os
import sys

if hasattr(sys.stdout, "reconfigure"):
    sys.stdout.reconfigure(encoding="utf-8", line_buffering=True, write_through=True)
if hasattr(sys.stderr, "reconfigure"):
    sys.stderr.reconfigure(encoding="utf-8", line_buffering=True, write_through=True)

current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

import CoolantLoop.model_collector_ring_full_ringhp as model


DEFAULT_CASE_NAME = "collector_ring_full_ringhp_geometry100hp_test"

RECT_LENGTH = 0.110
RECT_WIDTH = 0.040
RECT_AREA = RECT_LENGTH * RECT_WIDTH
RECT_WETTED_PERIMETER = 2.0 * (RECT_LENGTH + RECT_WIDTH)
RECT_DH = 4.0 * RECT_AREA / RECT_WETTED_PERIMETER

HP_TOTAL_COUNT = 100
HP_MULTIPLIERS_RING = [5 if idx in (0, 6, 12, 18) else 4 for idx in range(model.N_RING)]

HP_L_EVA = 0.100
HP_L_ABA = 0.0
HP_L_CON = 0.500
HP_R_VAPOR = 0.0075
HP_WICK_THICKNESS = 0.0005
HP_WALL_THICKNESS = 0.0010
HP_R_IN_WALL = HP_R_VAPOR + HP_WICK_THICKNESS
HP_R_OUT = HP_R_IN_WALL + HP_WALL_THICKNESS
HP_POROSITY = 0.6

FIN_AVERAGE_WIDTH = 0.020
FIN_THICKNESS = 0.0004
FIN_HEIGHT = FIN_AVERAGE_WIDTH
N_FIN_HEIGHT = 15

HEADER_WALL_THICKNESS = 0.002
HEADER_EQUIV_R_IN = RECT_WETTED_PERIMETER / (2.0 * math.pi)
HEADER_EQUIV_R_OUT = HEADER_EQUIV_R_IN + HEADER_WALL_THICKNESS


def apply_geometry_overrides(args):
    overrides = {
        "W_TOTAL": float(args.total_flow),
        "W_BRANCH_TOTAL": float(args.total_flow) / 3.0,
        "W_INLET_LEG_INIT": float(args.total_flow) / 3.0,
        "T_INIT": float(args.init_temp),
        "HP_INITIAL_TEMP": float(args.hp_init_temp),
        "AREA_RING": RECT_AREA,
        "DH_RING": RECT_DH,
        "PERIM_HEADER": RECT_WETTED_PERIMETER,
        "R_IN_RING": HEADER_EQUIV_R_IN,
        "R_OUT_RING": HEADER_EQUIV_R_OUT,
        "HP_COUNT_PER_SECTOR": HP_TOTAL_COUNT / 6.0,
        "HP_MULTIPLIERS_RING": list(HP_MULTIPLIERS_RING),
        "HP_MULTIPLIERS_SECTOR": list(HP_MULTIPLIERS_RING[: model.N_SECTOR]),
        "R_OUT_HP": HP_R_OUT,
        "R_IN_HP": HP_R_IN_WALL,
        "R_VAPOR_HP": HP_R_VAPOR,
        "L_EVA": HP_L_EVA,
        "L_ABA": HP_L_ABA,
        "L_CON": HP_L_CON,
        "HP_N_EVA": 1,
        "HP_N_ABA": 0,
        "HP_N_CON": 12,
        "HP_N_WICK": 1,
        "HP_N_WALL": 2,
        "POROSITY": HP_POROSITY,
        "THIN_FIN": FIN_THICKNESS,
        "FIN_HEIGHT": FIN_HEIGHT,
        "N_FIN_HEIGHT": N_FIN_HEIGHT,
        "HP_EMISSIVITY": float(args.hp_emissivity),
        "FIN_EMISSIVITY": float(args.fin_emissivity),
    }
    old_values = {name: getattr(model, name) for name in overrides}
    for name, value in overrides.items():
        setattr(model, name, value)
    return old_values


def restore_overrides(old_values):
    for name, value in old_values.items():
        setattr(model, name, value)


def print_geometry_summary(args):
    print("=" * 78)
    print("Geometry variant: rectangular collector ring + 100 heat pipes")
    print(f"  Ring rectangular section : {RECT_LENGTH * 1000:.1f} mm x {RECT_WIDTH * 1000:.1f} mm")
    print(f"  Ring flow area           : {RECT_AREA:.9e} m2")
    print(f"  Ring wetted perimeter    : {RECT_WETTED_PERIMETER:.9e} m")
    print(f"  Ring hydraulic diameter  : {RECT_DH:.9e} m")
    print(f"  Header equivalent radii  : {HEADER_EQUIV_R_IN:.9e} / {HEADER_EQUIV_R_OUT:.9e} m")
    print(f"  HP multipliers sum       : {sum(HP_MULTIPLIERS_RING)}")
    print(f"  HP multipliers           : {HP_MULTIPLIERS_RING}")
    print(f"  HP lengths eva/aba/con   : {HP_L_EVA:.3f} / {HP_L_ABA:.3f} / {HP_L_CON:.3f} m")
    print(f"  HP radii vapor/in/out    : {HP_R_VAPOR:.4f} / {HP_R_IN_WALL:.4f} / {HP_R_OUT:.4f} m")
    print(f"  HP porosity              : {HP_POROSITY:.3f}")
    print(f"  Fin width/thickness      : {FIN_AVERAGE_WIDTH:.4f} / {FIN_THICKNESS:.4f} m")
    print(f"  Emissivity HP/fin        : {args.hp_emissivity:.3f} / {args.fin_emissivity:.3f}")
    print("=" * 78)


def run_case(args):
    old_values = apply_geometry_overrides(args)
    try:
        print_geometry_summary(args)
        return model.run_case(
            case_name=args.case_name,
            t_end=args.t_end,
            min_dt=args.min_dt,
            max_dt=args.max_dt,
            safety_factor=args.safety_factor,
            inner_iter=args.inner_iter,
            print_every_time=args.print_every,
            csv_path=args.csv_path,
            restart_from=args.restart_from,
            restart_save_path=args.restart_save_path,
            restart_save_every=args.restart_save_every,
        )
    finally:
        restore_overrides(old_values)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Build and test a geometry variant with a rectangular collector ring and 100 heat pipes."
    )
    parser.add_argument("--total-flow", type=float, default=1.3)
    parser.add_argument("--init-temp", type=float, default=model.T_INLET)
    parser.add_argument("--hp-init-temp", type=float, default=800.0)
    parser.add_argument("--hp-emissivity", type=float, default=0.85)
    parser.add_argument("--fin-emissivity", type=float, default=0.85)
    parser.add_argument("--t-end", type=float, default=0.05)
    parser.add_argument("--case-name", default=DEFAULT_CASE_NAME)
    parser.add_argument("--min-dt", type=float, default=2.0e-4)
    parser.add_argument("--max-dt", type=float, default=5.0e-3)
    parser.add_argument("--safety-factor", type=float, default=0.1)
    parser.add_argument("--inner-iter", type=int, default=2)
    parser.add_argument("--print-every", type=float, default=0.05)
    parser.add_argument("--restart-save-every", type=float, default=0.0)
    parser.add_argument("--restart-from", default=None)
    parser.add_argument(
        "--csv-path",
        default=os.path.join(current_dir, f"{DEFAULT_CASE_NAME}_history.csv"),
    )
    parser.add_argument(
        "--restart-save-path",
        default=os.path.join(current_dir, f"{DEFAULT_CASE_NAME}_restart.npz"),
    )
    return parser.parse_args()


def main():
    run_case(parse_args())


if __name__ == "__main__":
    main()
