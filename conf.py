#!/usr/bin/env python3
"""
Generate conformers from a one-record SDF file.

The input SDF must contain exactly one valid molecular record with the intended
atom connectivity, bond orders, and formal charges.

Usage:
    python conf.py -i <filename.sdf>
"""

import argparse
import sys
from pathlib import Path

from rdkit import Chem
from openconf import ConformerConfig, generate_conformers


def get_parser():
    """Builds and returns the ArgumentParser."""
    parser = argparse.ArgumentParser(usage=__doc__)

    parser.add_argument(
        '--input',
        '-i',
        required=True,
        type=str,
        help='Path to a one-record input SDF file',
    )

    return parser


def load_one_sdf_molecule(input_path):
    """
    Load exactly one valid molecule from an SDF file.

    Explicit hydrogens are retained while reading so that the input graph can
    be checked before generating a hydrogen-suppressed SMILES for openconf.
    """
    supplier = Chem.SDMolSupplier(
        str(input_path),
        sanitize=True,
        removeHs=False,
        strictParsing=True,
    )

    mols = [mol for mol in supplier if mol is not None]

    if not mols:
        raise ValueError(
            f"RDKit could not read a valid molecular record from: {input_path}"
        )

    if len(mols) != 1:
        raise ValueError(
            f"Expected exactly one valid molecule in {input_path}, "
            f"but found {len(mols)}."
        )

    return mols[0]


def main(argv=None):
    if argv is None:
        argv = sys.argv[1:]

    parser = get_parser()
    args = parser.parse_args(argv)

    input_path = Path(args.input)

    if not input_path.is_file():
        parser.error(
            f"Input file does not exist or is not a regular file: {input_path}"
        )

    if input_path.suffix.lower() != ".sdf":
        parser.error(f"Input file must have a .sdf extension: {input_path}")

    print("Reading SDF file")

    try:
        mol = load_one_sdf_molecule(input_path)
    except ValueError as exc:
        parser.error(str(exc))

    print(f"Atoms: {mol.GetNumAtoms()}")
    print(f"Bonds: {mol.GetNumBonds()}")
    print(f"Conformers in input: {mol.GetNumConformers()}")
    print(f"Formal charge: {Chem.GetFormalCharge(mol)}")

    rawsmiles = Chem.MolToSmiles(mol, isomericSmiles=True)
    print("SMILES from input SDF:")
    print(rawsmiles)

    mol_no_h = Chem.RemoveHs(mol)
    smiles = Chem.MolToSmiles(mol_no_h, isomericSmiles=True)

    print("SMILES passed to openconf:")
    print(smiles)

    newmol = Chem.MolFromSmiles(smiles)

    if newmol is None:
        parser.error(
            "RDKit could not parse the SMILES generated from the input SDF."
        )

    config = ConformerConfig(
        max_out=100,
        pool_max=1000,
        n_steps=400,
        energy_window_kcal=12.0,
        seed_n_per_rotor=5,
        seed_prune_rms_thresh=0.5,
        do_final_refine=True,
        minimize_batch_size=8,
        parent_strategy="softmax",
        final_select="energy",
        use_low_mode_following=True
    )

    ensemble = generate_conformers(newmol, config=config)

    print(f"Generated {ensemble.n_conformers} conformers")
    print(ensemble.summary())

    output_path = input_path.with_name(
        f"{input_path.stem}_conformers.xyz"
    )

    ensemble.to_xyz(str(output_path))
    print(f"Wrote conformers to: {output_path}")


if __name__ == "__main__":
    main()
