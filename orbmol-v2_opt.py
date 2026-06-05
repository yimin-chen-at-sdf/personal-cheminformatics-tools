"""
Comments to be added
"""

import os
import argparse
from pathlib import Path
import ase
from ase.io import read, write
from ase.optimize import BFGS
from orb_models.forcefield import pretrained
from orb_models.forcefield.inference.calculator import ORBCalculator

def parse_args():
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("--device", "-d", required=True, choices=["cpu", "cuda"], help="Device to use: 'cpu' or 'cuda'")
    parser.add_argument("--precision", "-p", default="float32-high", choices=["float32-high", "float32-highest", "float64"], help="Calculation precision: 'float32-high', 'float32-highest', or 'float64' (default: float32-high)")
    parser.add_argument("--weights", "-w", type=str, help="Path to the predownloaded check point file")
    parser.add_argument("--input", "-i", required=True, type=str, help="Path of the input xyz file")
    parser.add_argument("--charge", "-c", default=0, type=int, help="Net charge of the molecule with the default value being zero")
    parser.add_argument("--multiplicity", "-m", default=1, type=int, help="Multiplicity of the molecule with the default value being one")
    parser.add_argument("--trajectory", "-t", action="store_true", help="Once this option is specified, trajectory of geometry optimization will be outputted")
    return parser.parse_args()

def resolve_weights(weights_arg):
    if weights_arg is None:
        print("No --weights or -w argument provided. The program might need to download check point file.")
        return None

    weights_path = Path(weights_arg)

    if not weights_path.is_file():
        print(f"Error: invalid checkpoint file path: {weights_path}", file=sys.stderr)
        sys.exit(1)

    print(f"Predownloaded check point file will be used: {weights_path}")
    return weights_path

def notify_user(args):
    print(f"Using device: {args.device}")
    print(f"Using precision: {args.precision}")
    weights_path = resolve_weights(args.weights)
    input_path = Path(args.input).resolve()
    output_dir = Path.cwd()
    if not input_path.exists():
        parser.error(f"Input file does not exist: {input_path}")
    if input_path.suffix.lower() != ".xyz":
        parser.error("Input file must have a .xyz extension")
    opt_filename = f"{input_path.stem}_opt{input_path.suffix}"
    opt_path = output_dir / opt_filename
    return input_path, opt_path, weights_path

def main():
    args = parse_args()
    input_path, opt_path, weights_path = notify_user(args)
    device = args.device
    if weights_path is None:
        orbff, atoms_adapter = pretrained.orbmol_v2(device=device, precision=args.precision)
    else:
        orbff, atoms_adapter = pretrained.orbmol_v2(weights_path=weights_path, device=device, precision=args.precision)
    calc = ORBCalculator(orbff, atoms_adapter=atoms_adapter, device=device)

    atoms=read(input_path)
    atoms.info["charge"] = args.charge
    atoms.info["spin"] = args.multiplicity
    atoms.calc = calc

    if args.trajectory:
        trj_path = opt_path.with_name(opt_path.name[:-len("_opt.xyz")] + "_trj.xyz")
        intermediate_path = opt_path.with_name("orb_optimized.traj")
        dyn = BFGS(atoms, trajectory=intermediate_path)
        dyn.run(fmax=0.01)
        write(opt_path, atoms)
        images = read(intermediate_path, index=":")
        write(trj_path, images)
        intermediate_path.unlink(missing_ok=True)
    else:
        dyn = dyn = BFGS(atoms)
        dyn.run(fmax=0.01)
        write(opt_path, atoms)

if __name__ == "__main__":
    main()
