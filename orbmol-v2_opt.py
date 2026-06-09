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
import csv
import numpy as np
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
    parser.add_argument("--export_csv", "-e", action="store_true", help="Once this option is specified, a csv file about geometry optimization will be outputted")
    parser.add_argument("--fmax", type=float, default=0.01, help="Threshold for maximum force acting on any atom of the system under investigation with the default value being 0.01")
    parser.add_argument("--maxcycles", type=int, default=1000, help="Maximum steps of geometry optimization with the default value being 1000")
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

def extract_energies_and_fmax(traj_path):
    frames = read(traj_path, index=":")
    energies = []
    fmax_list = []

    for atoms in frames:
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        fmax = float(np.linalg.norm(forces, axis=1).max())

        energies.append(float(energy))
        fmax_list.append(fmax)

    return energies, fmax_list

def run_continuous_optimization(atoms, input_path, opt_path, output_trajectory, fmax, maxcycles):
    intermediate_path = opt_path.with_name(opt_path.name[:-len("_opt.xyz")] + "_opt.traj")
    opt = Sella(atoms, order=0, internal=True, trajectory=os.fspath(intermediate_path))
    opt.run(fmax=fmax, steps=maxcycles)
    write(opt_path, atoms)
    energies, fmax_list = extract_energies_and_fmax(intermediate_path)
    last_energy = None
    dE = []
    for energy in energies:
        if last_energy is None:
            energy_change = "N/A"
        else:
            energy_change = float(energy - last_energy)
        dE.append(energy_change)
        last_energy = energy
    if output_trajectory:
        trj_path = opt_path.with_name(opt_path.name[:-len("_opt.xyz")] + "_trj.xyz")
        images = read(intermediate_path, index=":")
        write(trj_path, images)
    intermediate_path.unlink(missing_ok=True)
    return dE, fmax_list

def initialize_csv(csv_path):
    with open(csv_path, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["number", "dE", "fmax"])


def write_csv(start_step, dE_block, fmax_block, csv_path):
    with open(csv_path, "a", newline="") as file:
        writer = csv.writer(file)
        for i, (de, fmax) in enumerate(zip(dE_block, fmax_block)):
            de_out = de if de == "N/A" else f"{de:.8f}"
            fmax_out = f"{fmax:.8f}"
            writer.writerow([start_step + i, de_out, fmax_out])

def main():
    args = parse_args()
    input_path, opt_path, weights_path = notify_user(args)
    calc = set_calculator(args.device, args.precision, weights_path)
    atoms = set_atoms(input_path, args.charge, args.multiplicity, calc)

    dE, fmax_list = run_continuous_optimization(atoms, input_path, opt_path, args.trajectory, args.fmax, args.maxcycles)
    if args.export_csv:
        csv_path = opt_path.with_name(opt_path.name[:-len("_opt.xyz")] + "_opt.csv")
        initialize_csv(csv_path)
        write_csv(0, dE, fmax_list, csv_path)

if __name__ == "__main__":
    main()
