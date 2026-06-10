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

def build_parser():
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("--device", "-d", required=True, choices=["cpu", "cuda"], help="Device to use: 'cpu' or 'cuda'")
    parser.add_argument("--precision", "-p", default="float32-high", choices=["float32-high", "float32-highest", "float64"], help="Calculation precision: 'float32-high', 'float32-highest', or 'float64' (default: float32-high)")
    parser.add_argument("--weights", "-w", type=str, help="Path to the predownloaded check point file")
    parser.add_argument("--input", "-i", required=True, type=str, help="Path of the input xyz file")
    parser.add_argument("--charge", "-c", default=0, type=int, help="Net charge of the molecule with the default value being zero")
    parser.add_argument("--multiplicity", "-m", default=1, type=int, help="Multiplicity of the molecule with the default value being one")
    parser.add_argument("--trajectory", "-t", action="store_true", help="Once this option is specified, trajectory of geometry optimization will be outputted")
    parser.add_argument("--export_csv", "-e", action="store_true", help="Once this option is specified, a csv file about geometry optimization will be outputted")
    parser.add_argument("--fmax_threshold", type=float, default=0.01, help="Threshold for maximum force acting on any atom of the system under investigation with the default value being 0.01")
    parser.add_argument("--maxcycles", type=int, default=1000, help="Maximum steps of geometry optimization with the default value being 1000")
    parser.add_argument("--run_consecutive_optimization", action="store_true", help="Once this option is specified, consecutive geometry optimization will be performed, which consists of several blocks")
    parser.add_argument("--steps_per_block", type=int, default=argparse.SUPPRESS, help="The number of steps per block for running consecutive geometry optimization with the default value being 10")
    parser.add_argument("--energy_change_threshold", type=float, default=argparse.SUPPRESS, help="Threshold for energy change in consecutive geometry optimization with the default value being 3e-5")
    parser.add_argument("--nde_check", type=int, default=argparse.SUPPRESS, help="Number of steps involved in checking energy change in one block of calculation in consecutive geometry optimization with the default value being 3")
    return parser

def validate_dependencies_in_consecutive_optimization(args, parser):
    steps_specified = hasattr(args, "steps_per_block")
    energy_specified = hasattr(args, "energy_change_threshold")
    nde_specified = hasattr(args, "nde_check")

    if steps_specified and not args.run_consecutive_optimization:
        parser.error("--steps_per_block requires --run_consecutive_optimization")
    if energy_specified and not args.run_consecutive_optimization:
        parser.error("--energy_change_threshold requires --run_consecutive_optimization")
    if nde_specified and not args.run_consecutive_optimization:
        parser.error("--nde_check requires --run_consecutive_optimization")

def apply_defaults_in_consecutive_optimization(args):
    if not hasattr(args, "steps_per_block"):
        args.steps_per_block = 10
    if not hasattr(args, "energy_change_threshold"):
        args.energy_change_threshold = 3e-5
    if not hasattr(args, "nde_check"):
        args.nde_check = 3

def validate_values(args, parser):
    if args.maxcycles <= 0:
        parser.error("--maxcycles must be greater than 0")
    if args.steps_per_block <= 1:
        parser.error("--steps_per_block must be greater than 1")
    if args.steps_per_block >= args.maxcycles:
        parser.error("--steps_per_block must be smaller than --maxcycles")
    if args.energy_change_threshold <= 0:
        parser.error("--energy_change_threshold must be a positive float")
    if args.nde_check <= 0:
        parser.error("--nde_check must be greater than 0")
    if args.nde_check > args.steps_per_block:
        parser.error("--nde_check must be smaller than or equal to --steps_per_block")

def parse_args(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)

    validate_dependencies_in_consecutive_optimization(args, parser)
    apply_defaults_in_consecutive_optimization(args)
    validate_values(args, parser)

    return args

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
    atoms = read(input_path, format='xyz')
    atoms.info["charge"] = charge
    atoms.info["spin"] = multiplicity
    atoms.calc = calc
    return atoms

def extract_energies_and_fmax(traj_path, iblock):
    frames = read(traj_path, index=":")
    energies = []
    fmax_list = []

    for atoms in frames:
        energy = atoms.get_potential_energy()
        forces = atoms.get_forces()
        fmax = float(np.linalg.norm(forces, axis=1).max())

        energies.append(float(energy))
        fmax_list.append(fmax)

    if iblock > 0:
        del energies[0]
        del fmax_list[0]

    return energies, fmax_list

def initialize_csv(csv_path):
    with open(csv_path, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["step", "dE", "fmax"])


def write_csv(start_step, dE_block, fmax_block, csv_path):
    with open(csv_path, "a", newline="") as file:
        writer = csv.writer(file)
        for i, (de, fmax) in enumerate(zip(dE_block, fmax_block)):
            de_out = de if de == "N/A" else f"{de:.8f}"
            fmax_out = f"{fmax:.8f}"
            writer.writerow([start_step + i, de_out, fmax_out])

def check_energy_force_convergence(dE_block, fmax_block, energy_change_threshold, fmax_threshold, nde_check):
    last_dE_values = dE_block[-nde_check:]

    return (
        all(abs(x) < energy_change_threshold for x in last_dE_values)
        and fmax_block[-1] < fmax_threshold
    )

def combine_xyz_files(opt_path, iblock):
    base = opt_path.stem.removesuffix("_opt")
    merged = opt_path.with_name(f"{base}_trj_0000.xyz")
    final_path = opt_path.with_name(f"{base}_trj.xyz")

    if iblock == 0:
        merged.rename(final_path)
        return final_path

    for i in range(1, iblock+1):
        xyzpath = opt_path.with_name(f"{base}_trj_{i:04d}.xyz")
        with merged.open("a", encoding="utf-8") as fout, xyzpath.open("r", encoding="utf-8") as fin:
            fout.write(fin.read())
        xyzpath.unlink(missing_ok=True)

    merged.rename(final_path)
    return final_path

def perform_consecutive_optimization(atoms, opt_path, output_trajectory, fmax_threshold, steps_per_block, nblocks, export_csv, energy_change_threshold, nde_check):
    dE = []
    fmax_history = []
    last_energy = None

    if export_csv:
        csv_path = opt_path.with_name(opt_path.name[:-len("_opt.xyz")] + "_opt.csv")
        initialize_csv(csv_path)

    for iblock in range(nblocks):
        intermediate_path = opt_path.with_name(f"{opt_path.stem}_{iblock:04d}.traj")
        tighter_fmax = fmax_threshold / 10.0
        opt = Sella(atoms, order=0, internal=True, trajectory=os.fspath(intermediate_path))
        opt.run(fmax=tighter_fmax, steps=steps_per_block)

        energies_block, fmax_block = extract_energies_and_fmax(intermediate_path, iblock)
        dE_block = []
        for energy in energies_block:
            if last_energy is None:
                de = "N/A"
            else:
                de = float(energy - last_energy)
            dE.append(de)
            dE_block.append(de)
            last_energy = energy
        fmax_history.extend(fmax_block)
        
        if output_trajectory:
            base = opt_path.stem.removesuffix("_opt")
            trj_path = opt_path.with_name(f"{base}_trj_{iblock:04d}.xyz")
            images = read(intermediate_path, index=":")
            images = images[1:] if iblock > 0 else images
            write(trj_path, images)
        intermediate_path.unlink(missing_ok=True)

        if export_csv:
            start_step = len(fmax_history) - len(fmax_block)
            write_csv(start_step, dE_block, fmax_block, csv_path)

        # Stop condition 1:
        # This block produced fewer frames than the requested block size.
        if len(fmax_block) < steps_per_block:
            stop_reason = (
                f"Consecutive geometry optimization ended: block {iblock:04d} produced "
                f"{len(fmax_block)} steps, fewer than steps_per_block={steps_per_block}."
            )
            print(stop_reason)
            break

        # Stop condition 2:
        # The last nde_check dE values are all small enough and the last fmax is small enough.
        if check_energy_force_convergence(dE_block, fmax_block, energy_change_threshold, fmax_threshold, nde_check):
            stop_reason = (
                f"Consecutive geometry optimization ended: the absolute values of the last "
                f"{nde_check} dE values are all smaller than {energy_change_threshold} "
                f"and the last fmax is smaller than {fmax_threshold}."
            )
            print(stop_reason)
            break
    if output_trajectory:
        combine_xyz_files(opt_path, iblock)
    write(opt_path, atoms, format='xyz')
    return dE, fmax_history, last_energy

def perform_continuous_optimization(atoms, opt_path, output_trajectory, fmax_threshold, maxcycles):
    intermediate_path = opt_path.with_name(opt_path.name[:-len("_opt.xyz")] + "_opt.traj")
    opt = Sella(atoms, order=0, internal=True, trajectory=os.fspath(intermediate_path))
    opt.run(fmax=fmax_threshold, steps=maxcycles)

    energies, fmax_list = extract_energies_and_fmax(intermediate_path, 1)
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
    write(opt_path, atoms, format='xyz')
    return dE, fmax_list, last_energy

def main():
    args = parse_args()
    input_path, opt_path, weights_path = notify_user(args)
    calc = set_calculator(args.device, args.precision, weights_path)
    atoms = set_atoms(input_path, args.charge, args.multiplicity, calc)

    if args.run_consecutive_optimization:
        nblocks = args.maxcycles // args.steps_per_block
        if args.maxcycles % args.steps_per_block > 0:
            nblocks += 1
        dE, fmax_history, last_energy = perform_consecutive_optimization(atoms, opt_path, args.trajectory, args.fmax_threshold, args.steps_per_block, nblocks, args.export_csv, args.energy_change_threshold, args.nde_check)
        steps = len(fmax_history) - 1
        print(f"The final energy is {last_energy:.8f} after {steps} steps of consecutive geometry optimization.")
    else:
        dE, fmax_list, last_energy = perform_continuous_optimization(atoms, opt_path, args.trajectory, args.fmax_threshold, args.maxcycles)
        steps = len(fmax_list) - 1
        print(f"The final energy is {last_energy:.8f} after {steps} steps of continuous geometry optimization.")
        if args.export_csv:
            csv_path = opt_path.with_name(opt_path.name[:-len("_opt.xyz")] + "_opt.csv")
            initialize_csv(csv_path)
            write_csv(0, dE, fmax_list, csv_path)

if __name__ == "__main__":
    main()
