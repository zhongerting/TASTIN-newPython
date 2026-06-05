import argparse
import os
import sys

current_dir = os.path.dirname(os.path.abspath(__file__))
root_dir = os.path.abspath(os.path.join(current_dir, ".."))
if root_dir not in sys.path:
    sys.path.insert(0, root_dir)

from Materials.Solids.KHP import PotassiumHP
from Materials.Solids.WickMaterial import WickMaterial

import CoolantLoop.run_collector_ring_full_ringhp_geometry100hp_test as geom_case


model = geom_case.model

DEFAULT_CASE_NAME = "collector_ring_full_ringhp_geometry100hp_potassium_test"


def apply_material_overrides():
    old_values = {
        "mat_hp_fluid": model.mat_hp_fluid,
        "mat_wick": model.mat_wick,
    }

    hp_fluid = PotassiumHP(name="HP_Fluid_K")
    model.mat_hp_fluid = hp_fluid
    model.mat_wick = WickMaterial(
        name="WickMaterial_K",
        solid_mat=model.mat_wall,
        fluid_mat=hp_fluid,
        porosity=model.POROSITY,
        r_vapor=model.R_VAPOR_HP,
        r_in_wall=model.R_IN_HP,
    )
    return old_values


def restore_material_overrides(old_values):
    for name, value in old_values.items():
        setattr(model, name, value)


def print_case_summary():
    print("Heat-pipe working fluid : PotassiumHP")
    print("Heat-pipe fluid name    : HP_Fluid_K")
    print("Wick material name      : WickMaterial_K")
    print(f"Wick porosity/radii     : {model.POROSITY:.4f} / {model.R_VAPOR_HP:.6f} / {model.R_IN_HP:.6f} m")


def run_case(args):
    geometry_old_values = geom_case.apply_geometry_overrides(args)
    material_old_values = None
    try:
        material_old_values = apply_material_overrides()
        geom_case.print_geometry_summary(args)
        print_case_summary()
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
        if material_old_values is not None:
            restore_material_overrides(material_old_values)
        geom_case.restore_overrides(geometry_old_values)


def parse_args():
    parser = argparse.ArgumentParser(
        description=(
            "Build and test a potassium heat-pipe variant of the rectangular "
            "collector-ring + 100 heat-pipe case."
        )
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
