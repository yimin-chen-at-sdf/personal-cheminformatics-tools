#!/usr/bin/env python3
"""
This script performs geometry optimization with OrbMol-v2. This script relies 
on ase (Atomic Simulation Environment), orb-models, and sella. The user 
supplies an input.xyz and the result of geometry optimization will be 
input_opt.xyz. If the user chooses to preserve the process of geometry 
optimization with the "-t" argument, an input_trj.xyz will be produced as well.
If the user uses CPU to run the calculations, OMP_NUM_THREADS and 
MKL_NUM_THREADS environment variables should be set prior to running this 
script, otherwise the optimization process can be very slow. The check point 
file can be obtained from 
https://huggingface.co/Orbital-Materials/orbmol-v2/tree/main prior to running 
calculations.

Optimize a neutral closed-shell molecule:
  python orbmol-v2_opt.py -d cpu -i input.xyz
Optimize a neutral closed-shell molecule and preserve the process of geometry 
optimization:
  python orbmol-v2_opt.py -d cpu -i input.xyz -t
Optimize a neutral closed-shell molecule and export a csv file about the 
process of geometry optimization:
  python orbmol-v2_opt.py -d cpu -i input.xyz -e
Use a predownloaded check point file to optimize a neutral closed-shell 
molecule:
  python orbmol-v2_opt.py -d cpu -w /path/to/check/point/file -i input.xyz
Optimize a closed-shell molecule with +1 charge
  python orbmol-v2_opt.py -d cpu -i input.xyz -c 1
Optimize a neutral radical species (S = 1/2, 2S + 1 = 2):
  python orbmol-v2_opt.py -d cpu -i input.xyz -m 2
Use NVIDIA GPU rather than CPU to optimize a neutral closed-shell molecule:
  python orbmol-v2_opt.py -d cuda -i input.xyz
Optimize a neutral closed-shell molecule with the precision being float32-
highest instead of the default float32-high:
  python orbmol-v2_opt.py -d cpu -p float32-highest -i input.xyz
Optimize a neutral closed-shell molecule and change the threshold for maximum 
force acting on any atom to 0.02:
  python orbmol-v2_opt.py -d cpu -i input.xyz --fmax_threshold 0.02
Optimize a neutral closed-shell molecule and change maximum steps of geometry 
optimization to 250:
  python orbmol-v2_opt.py -d cpu -i input.xyz --maxcycles 250
Use an experimental feature to optimize a neutral closed-shell molecule with 
several blocks of geometry optimization:
  python orbmol-v2_opt.py -d cpu -i input.xyz --run_consecutive_optimization

Advanced usage with bond constraints:
When any bond constraint is specified, the process of geometry optimization 
will be preserved by default. Only one-based indexing of atoms is allowed when 
sepecifying bond constraints. The threshold for maximum force acting on any 
atom in this case is 1e-3 instead of the default 0.01. The maximum steps of 
geometry optimization in this case is 250 instead of the default 1000. 
Constrain the bond between atom 1 and atom 2:
  python orbmol-v2_opt.py -d cpu -i input.xyz --fix_bond 1 2
Constrain the bond between atom 1 and atom 2 to 1.5 Angstrom:
  python orbmol-v2_opt.py -d cpu -i input.xyz --fix_bond 1 2 --target 1.5
NOTE: The command above does not change the geometry immediately.
Constrain the bond between atom 1 and atom 2 along with the bond between atom 7
and atom 8:
  python orbmol-v2_opt.py -d cpu -i input.xyz --fix_bond 1 2 7 8
Constrain the bond between atom 1 and atom 2 while constrain the bond between 
atom 7 and atomm 8 to 1.5 Angstrom:
  python orbmol-v2_opt.py -d cpu -i input.xyz --fix_bond 1 2 7 8 --target C 1.5
NOTE: The command above does not change the geometry immediately.

The author only tested this script on python 3.13.13. The author used ase 
3.28.0, orb-models 0.7.0, and sella 2.4.2.
"""

import os
import platform
import argparse
from pathlib import Path
import csv
import numpy as np
from ase.io import read, write
from sella import Sella, Constraints
from orb_models.forcefield import pretrained
from orb_models.forcefield.inference.calculator import ORBCalculator

def positive_int(value):
    """
    Custom type function to validate positive integers.

    Args:
      value (str): One-based indexing of atom

    Returns:
      int
    """
    try:
        value = int(value)
    except ValueError:
        raise argparse.ArgumentTypeError(f"{value!r} is not an integer")

    if value <= 0:
        raise argparse.ArgumentTypeError(f"{value} is not a positive integer")

    return value


def target_value(value):
    """
    Accept either 'C' or a floating-point target value.

    Args:
      value (str): 'C' or bond length

    Returns:
      'C' for a bond length retaining current value
      float for an explicitly specified bond length
    """
    if value == "C":
        return "C"

    try:
        return float(value)
    except ValueError:
        raise argparse.ArgumentTypeError(
            f"{value!r} must be 'C' or a floating-point number"
        )

def build_parser():
    parser = argparse.ArgumentParser(usage=__doc__)
    parser.add_argument("--device", "-d", required=True, choices=["cpu", "cuda"], help="Device to use: 'cpu' or 'cuda'")
    parser.add_argument("--precision", "-p", default="float32-high", choices=["float32-high", "float32-highest", "float64"], help="Calculation precision: 'float32-high', 'float32-highest', or 'float64' (default: float32-high)")
    parser.add_argument("--weights", "-w", type=str, help="Path to the predownloaded check point file")
    parser.add_argument("--input", "-i", required=True, type=str, help="Path to the input xyz file")
    parser.add_argument("--charge", "-c", type=int, default=0, help="Net charge of the molecule with the default value being zero")
    parser.add_argument("--multiplicity", "-m", type=int, default=1, help="Multiplicity of the molecule with the default value being one")
    parser.add_argument("--trajectory", "-t", action="store_true", help="Once this option is specified, trajectory of geometry optimization will be outputted. If the user enforces constraint on some bond, trajectory of geometry optimization will be outputted by default even without specifying this option.")
    parser.add_argument("--export_csv", "-e", action="store_true", help="Once this option is specified, a csv file about geometry optimization will be outputted")
    parser.add_argument("--fmax_threshold", type=float, default=0.01, help="Threshold for maximum force acting on any atom of the system under investigation with the default value being 0.01. However, during geometry optimization with bond constraint, the default value will be 1e-3 and cannot be changed without editing the codes.")
    parser.add_argument("--maxcycles", type=int, default=1000, help="Maximum steps of geometry optimization with the default value being 1000. However, during geometry optimization with bond constraint, the default value will be 300 and cannot be changed without editing the codes.")
    parser.add_argument("--run_consecutive_optimization", action="store_true", help="Once this option is specified, consecutive geometry optimization will be performed, which consists of several blocks")
    parser.add_argument("--steps_per_block", type=int, default=argparse.SUPPRESS, help="The number of steps per block for running consecutive geometry optimization with the default value being 10")
    parser.add_argument("--energy_change_threshold", type=float, default=argparse.SUPPRESS, help="Threshold for energy change in consecutive geometry optimization with the default value being 3e-5")
    parser.add_argument("--nde_check", type=int, default=argparse.SUPPRESS, help="Number of steps involved in checking energy change in one block of calculation in consecutive geometry optimization with the default value being 3")
    parser.add_argument("--fix_bond", nargs="+", type=positive_int, help="Atom-index pairs for fixed bonds, e.g. --fix_bond 1 2 8 7. One-based indexing should be used.")
    parser.add_argument("--target", nargs="+", type=target_value, default=argparse.SUPPRESS, help="One target per fixed bond: C or a floating-point value. Here C means constant value. The unit or dimension here is Angstrom.")
    return parser

def validate_dependencies_in_consecutive_optimization(args, parser):
    """
    This function forbids the user from abusing some setting of the 
    experimental feature named consecutive optimization.
    """
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
    """
    Set the number of steps per block for running consecutive geometry 
    optimization to 10. Set the threshold for energy change in consecutive 
    geometry optimization to 3e-5. Set the number of steps used for checking 
    energy change convergence to 3.
    """
    if not hasattr(args, "steps_per_block"):
        args.steps_per_block = 10
    if not hasattr(args, "energy_change_threshold"):
        args.energy_change_threshold = 3e-5
    if not hasattr(args, "nde_check"):
        args.nde_check = 3

def validate_values(args, parser):
    """
    This function forbids the user from wrongly setting some parameters for 
    geometry optimization.
    """
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

    return parser, args

def check_cpu_environment():
    """Report Linux status and OpenMP/MKL thread-variable status."""
    if platform.system() != "Linux":
        print("The operating system is not Linux.")
        return
    print("The operating system is Linux.")
    omp_num_threads = os.environ.get("OMP_NUM_THREADS")
    mkl_num_threads = os.environ.get("MKL_NUM_THREADS")
    if omp_num_threads is None:
        print("The OMP_NUM_THREADS environment variable has not been set.")
    else:
        print("The OMP_NUM_THREADS environment variable has been set.")
    if mkl_num_threads is None:
        print("The MKL_NUM_THREADS environment variable has not been set.")
    else:
        print("The MKL_NUM_THREADS environment variable has been set.")
    if omp_num_threads is None and mkl_num_threads is None:
        print("This can slow down the calculations in the next step. You are supposed to stop the program, set the two environment variables, and rerun the program.")
    elif omp_num_threads is not None and mkl_num_threads is not None:
        if omp_num_threads == mkl_num_threads:
            print("They are equal to each other.")
        else:
            print("They are not equal to each other.")

def resolve_weights(weights_arg):
    """
    This function determines whether predownloaded check point file will be 
    used.

    Args:
      weights_arg (str): path to the predownloaded check point file

    Returns:
      path object
    """
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
    """
    This function notifies the user about some settings of the calculation to 
    be performed.
    """
    print(f"Using device: {args.device}")
    if args.device == "cpu":
        check_cpu_environment()
    print(f"Using precision: {args.precision}")
    if args.fix_bond is not None:
        print("Geometry optimization will be performed with the constraints specified by the user.")

    input_path = Path(args.input).resolve()
    if not input_path.exists():
        parser.error(f"Input file does not exist: {input_path}")
    if input_path.suffix.lower() != ".xyz":
        parser.error("Input file must have a .xyz extension")
    if len(input_path.suffixes) > 1:
        parser.error("Input file has more than one extension")
    with input_path.open("r", encoding="utf-8") as xyz_file:
        number_of_atoms = int(xyz_file.readline().strip())
    print(f"The system under investigation has {number_of_atoms} atoms")

    output_dir = Path.cwd()
    opt_filename = f"{input_path.stem}_opt{input_path.suffix}"
    opt_path = output_dir / opt_filename

    weights_path = resolve_weights(args.weights)

    return input_path, number_of_atoms, opt_path, weights_path

def validate_constraints(args, parser):
    """Validate constraint-related arguments."""
    has_fix_bond = args.fix_bond is not None
    has_target = hasattr(args, "target")

    # The user cannot supply "--target" argument without supplying 
    # "--fix_bond" argument.
    if has_target and not has_fix_bond:
        parser.error("--target can only be specified together with --fix_bond")

    # The user does not supply any "--fix_bond" argument.
    if not has_fix_bond:
        return args

    # --fix_bond must contain complete atom pairs.
    if len(args.fix_bond) % 2 != 0:
        parser.error("The last atom in --fix_bond does not have bond specified")

    number_of_bonds = len(args.fix_bond) // 2

    # One --target input is required for each atom pair.
    if has_target and len(args.target) != number_of_bonds:
        parser.error(
            "The number of --target values must equal the number "
            "of fixed bonds"
        )

    return args

def validate_bond_atom_indices(number_of_atoms, atom_indices, argument_name, parser):
    """
    Ensure that atom indices supplied by the user exist in atoms.

    Args:
      number_of_atoms (int): The number of atoms in the system under 
      investigation.
      atom_indices (list[int]): One-based atom indices supplied by the user.
      None means that the corresponding argument was not supplied.
      argument_name (str): Name shown in error messages.
      parser (argparse.ArgumentParser): Used for standard argparse-style error 
      messages.
    """
    if atom_indices is None:
        return

    largest_index = max(atom_indices)

    if largest_index > number_of_atoms:
        parser.error(
            f"Atom index {largest_index} in {argument_name} exceeds "
            f"the number of atoms in the system ({number_of_atoms})"
        )

def set_calculator(device, precision, weights_path):
    """
    This functions sets a calculator compatible with ASE.

    Args:
      device (str): "cpu" or "cuda"
      precision (str): "float32-high", "float32-highest", or "float64"
      weights_path (path object): path to the predownloaded check point file

    Returns:
      an ASE-compatible calculator
    """
    if weights_path is None:
        orbff, atoms_adapter = pretrained.orbmol_v2(device=device, precision=precision)
    else:
        orbff, atoms_adapter = pretrained.orbmol_v2(weights_path=weights_path, device=device, precision=precision)
    return ORBCalculator(orbff, atoms_adapter=atoms_adapter, device=device)

def set_atoms(input_path, charge, multiplicity, calc):
    """
    This function reads the user-specified xyz file and sets the charge and 
    multiplicity.

    Args:
      input_path (path object): A path to the xyz file
      charge (int): Net charge of the molecule
      multiplicity (int): Multiplicity of the molecule
      calc: Calculator

    Returns:
      The molecule to be optimized
    """
    atoms = read(input_path, format='xyz')
    atoms.info["charge"] = charge
    atoms.info["spin"] = multiplicity
    atoms.calc = calc
    return atoms

def initialize_csv(csv_path):
    with open(csv_path, "w", newline="") as file:
        writer = csv.writer(file)
        writer.writerow(["step", "dE", "fmax"])

def one_based_to_zero_based(atom_indices):
    """
    Convert one-based atom indices to zero-based atom indices.

    Args:
      indices (list[int]): One-based atom indices

    Returns:
      zero-based atom indices
    """
    if atom_indices is None:
        return None

    return [atom_index - 1 for atom_index in atom_indices]

def set_sella_optimizer(atoms, traj_path, fixed_bond_pairs=None, target_list=None):
    """
    Create a Sella geometry optimizer.

    Args:
      atoms (ase.Atoms): The molecule to be optimized
      traj_path (path object): A file with its extension being traj which 
      stores the process of geometry optimization
      fixed_bond_pairs (list[int]): atom pairs of fixed bond
      target_list (list[float | str]): If the element is a float number, the 
      bond distance will be changed to that value in Anstrom during geometry
      optimization. If the element is 'C', the bond distance will remain the 
      same value.

    Returns:
      Sella optimizer
    """
    # The user does not supply any "--fix_bond" argument. Geometry optimization
    # will be performed without constraints.
    if fixed_bond_pairs is None:
        return Sella(atoms, order=0, internal=True, trajectory=traj_path)

    cons = Constraints(atoms)
    # Sella requires zero-based indexing of atom pairs while the atom pairs 
    # supplied by the user has one-based indexing.
    zero_based_indices = one_based_to_zero_based(fixed_bond_pairs)
    constraint_pairs = None
    if zero_based_indices is not None:
        constraint_pairs = list(
            zip(zero_based_indices[::2], zero_based_indices[1::2])
        )
    # The user supplies "--fix_bond" argument without "--target" argument.
    if target_list is None:
        for bond in constraint_pairs:
            cons.fix_bond(bond)

    # The user supplies "--fix_bond" argument along with "--target" argument.
    else:
        for bond, target in zip(constraint_pairs, target_list):
            if target == "C":
                cons.fix_bond(bond)
            else:
                cons.fix_bond(bond, target=target)

    # Geometry optimization will be performed with constraints.
    return Sella(atoms, order=0, constraints=cons, trajectory=traj_path)

def extract_energies_and_fmax(traj_path, iblock):
    """
    This function read a file with its extension being traj and outputs 
    energies along with maximum force acting on any atom of the system under 
    investigation.

    Args:
      traj_path (path object): A file with its extension being traj which 
      stores the process of geometry optimization
      iblock (int): If it is zero, then the first frame will be processed. If 
      it is greater than zero, then the first frame will be discarded.

    Returns:
      energies (list[float]): list of energy
      fmax_list (list[float]): list of maximum force acting on any atom
    """
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

def write_csv(start_step, dE_block, fmax_block, csv_path):
    """
    This function exports the process of geometry optimization to a csv file.

    Args:
      start_step (int): Zero-based index for the step of geometry optimization
      dE_block (list[float | str]): If the element is "N/A", then the 
      corresponding step should be the first step in geometry optimization. If 
      the element is a float number, it should be the energy change in geometry
      optimization.
      fmax_block (list[float]): List of maximum force acting on any atom
      csv_path (path object): path to the csv file being written
    """
    with open(csv_path, "a", newline="") as file:
        writer = csv.writer(file)
        for i, (de, fmax) in enumerate(zip(dE_block, fmax_block)):
            de_out = de if de == "N/A" else f"{de:.8f}"
            fmax_out = f"{fmax:.8f}"
            writer.writerow([start_step + i, de_out, fmax_out])

def check_energy_force_convergence(dE_block, fmax_block, energy_change_threshold, fmax_threshold, nde_check):
    """
    This function tells whether the consecutive geometry optimization has 
    reached convergence. The absolute values of the energy changes in the last 
    three steps should be smaller than certain value and the maximum force 
    acting on any atom in the last step should be smaller than certain value.

    Args:
      dE_block (list[float | str]): If the element is "N/A", then the 
      corresponding step should be the first step in geometry optimization. If 
      the element is a float number, it should be the energy change in geometry
      optimization.
      fmax_block (list[float]): List of maximum force acting on any atom
      energy_change_threshold (float): The threshold for energy changes in the last 
      several steps. The default value is 3e-5 in consecutive geometry 
      optimization.
      fmax_threshold (float): Threshold for maximum force acting on any atom. 
      The default value is 0.01 in consecutive geometry optimization.
      nde_check (int): The number of steps used for checking energy change 
      convergence in geometry optimization. The default value is 3.
    """
    last_dE_values = dE_block[-nde_check:]

    return (
        all(abs(x) < energy_change_threshold for x in last_dE_values)
        and fmax_block[-1] < fmax_threshold
    )

def combine_xyz_files(opt_path, iblock):
    """
    This functions combines multiple xyz files into one xyz file. This function
    is used in consecutive geometry optimization where the process of geometry 
    optimization is generated block by block, which results in multiple xyz 
    files.

    Args:
      opt_path (path object): Path to the xyz file storing the optimized 
      geometry
      iblock (int): number of xyz files to be combined minus one

    Returns:
      a path object of the generated xyz file
    """
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
        opt = set_sella_optimizer(atoms, traj_path=os.fspath(intermediate_path), fixed_bond_pairs=None, target_list=None)
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

def perform_continuous_optimization(atoms, opt_path, output_trajectory, fmax_threshold, maxcycles, fixed_bond_pairs, target_list):
    intermediate_path = opt_path.with_name(opt_path.name[:-len("_opt.xyz")] + "_opt.traj")
    opt = set_sella_optimizer(atoms=atoms, traj_path=os.fspath(intermediate_path), fixed_bond_pairs=fixed_bond_pairs, target_list=target_list)
    if fixed_bond_pairs is None:
        opt.run(fmax=fmax_threshold, steps=maxcycles)
    else:
        opt.run(fmax=1e-3, steps=250)

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

    if output_trajectory or fixed_bond_pairs is not None:
        trj_path = opt_path.with_name(opt_path.name[:-len("_opt.xyz")] + "_trj.xyz")
        images = read(intermediate_path, index=":")
        write(trj_path, images)
    intermediate_path.unlink(missing_ok=True)
    write(opt_path, atoms, format='xyz')
    return dE, fmax_list, last_energy

def main():
    parser, args = parse_args()
    input_path, number_of_atoms, opt_path, weights_path = notify_user(args)
    args = validate_constraints(args, parser)
    validate_bond_atom_indices(number_of_atoms=number_of_atoms, atom_indices=args.fix_bond, argument_name="--fix_bond", parser=parser)
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
        dE, fmax_list, last_energy = perform_continuous_optimization(atoms, opt_path, args.trajectory, args.fmax_threshold, args.maxcycles, args.fix_bond, getattr(args, "target", None))
        steps = len(fmax_list) - 1
        print(f"The final energy is {last_energy:.8f} after {steps} steps of continuous geometry optimization.")
        if args.export_csv:
            csv_path = opt_path.with_name(opt_path.name[:-len("_opt.xyz")] + "_opt.csv")
            initialize_csv(csv_path)
            write_csv(0, dE, fmax_list, csv_path)

if __name__ == "__main__":
    main()
