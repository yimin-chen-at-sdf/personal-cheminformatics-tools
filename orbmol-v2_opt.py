"""
This script performs geometry optimization with OrbMol-v2. This script relies 
on ase (Atomic Simulation Environment) and orb-models. The user supplies an 
input.xyz and the result of geometry optimization with be input_opt.xyz. If you
choose to preserve the process of geometry optimization with the "-t" argument,
an input_trj.xyz will be produced as well. If you use CPU to run the 
calculations, do not forget to set OMP_NUM_THREADS and MKL_NUM_THREADS 
environment variables prior to running this script, otherwise the optimization 
process can be very slow. The check point file can be obtained from 
https://huggingface.co/Orbital-Materials/orbmol-v2/tree/main prior to running 
calculations.

Optimize a neutral closed-shell molecule:
  python orbmol-v2_opt.py -d cpu -i input.xyz
Use a predownloaded check point file to optimize a neutral closed-shell 
molecule:
  python orbmol-v2_opt.py -d cpu -w /path/to/check/point/file -i input.xyz
Optimize a closed-shell molecule with +1 charge
  python orbmol-v2_opt.py -d cpu -i input.xyz -c 1
Optimize a neutral radical species (S = 1/2, 2S + 1 = 2):
  python orbmol-v2_opt.py -d cpu -i input.xyz -m 2
Optimize a neutral closed-shell molecule and preserve the process of geometry 
optimization:
  python orbmol-v2_opt.py -d cpu -i input.xyz -t
Use NVIDIA GPU rather than CPU to optimize a neutral closed-shell molecule:
  python orbmol-v2_opt.py -d cuda -i input.xyz
Optimize a neutral closed-shell molecule with the precision being float32-
highest instead of the default float32-high:
  python orbmol-v2_opt.py -d cpu -p float32-highest -i input.xyz
"""

import os
import argparse
from pathlib import Path
from ase.io import read, write
from sella import Sella
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
        print("No --weights or -w argument is provided. The program might need to download check point file.")
        return None

    weights_path = Path(weights_arg).resolve()

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

def set_calculator(device, precision, weights_path):
    if weights_path is None:
        orbff, atoms_adapter = pretrained.orbmol_v2(device=device, precision=precision)
    else:
        orbff, atoms_adapter = pretrained.orbmol_v2(weights_path=weights_path, device=device, precision=precision)
    return ORBCalculator(orbff, atoms_adapter=atoms_adapter, device=device)

def set_atoms(input_path, charge, multiplicity, calc):
    atoms = read(input_path)
    atoms.info["charge"] = charge
    atoms.info["spin"] = multiplicity
    atoms.calc = calc
    return atoms

def main():
    args = parse_args()
    input_path, opt_path, weights_path = notify_user(args)
    calc = set_calculator(args.device, args.precision, weights_path)
    atoms = set_atoms(input_path, args.charge, args.multiplicity, calc)

    if args.trajectory:
        trj_path = opt_path.with_name(opt_path.name[:-len("_opt.xyz")] + "_trj.xyz")
        intermediate_path = opt_path.with_name(opt_path.name[:-len("_opt.xyz")] + "_opt.traj")
        dyn = Sella(atoms, order=0, internal=True, trajectory=os.fspath(intermediate_path))
        dyn.run(1e-3, 1000)
        write(opt_path, atoms)
        images = read(intermediate_path, index=":")
        write(trj_path, images)
        intermediate_path.unlink(missing_ok=True)
    else:
        dyn = Sella(atoms, order=0, internal=True)
        dyn.run(1e-3, 1000)
        write(opt_path, atoms)

if __name__ == "__main__":
    main()
