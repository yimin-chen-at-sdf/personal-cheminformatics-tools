"""
This script employs xyz2mol developed by Jan Halborg Jensen group at University
of Copenhagen to transform the first frame of an xyz file to a SMILES string. 
This script requires RDKit which should be at least version 2022.09.1. Cite the
following paper if you use this script in any published work.
Y. Kim, W. Y. Kim, Bull. Korean Chem. Soc. 2015, 36(7), 1769–1777

For a neutral molecule, generate a canonical SMILES string with the following 
command:
  python xyz2smiles.py -i <filename>
For a molecue whose charge is -1, generate a canonical SMILES string with the 
following command:
  python xyz2smiles.py -i <filename> -c -1
For an aromatic hydrocarbon, generate a kekule SMILES string with the following
command:
  python xyz2smiles.py -i <filename> -k
For an aromatic hydrocarbon whose charge is -1, generate a kekule SMILES string
with the following command:
  python xyz2smiles.py -i <filename> -c -1 -k
"""

import argparse
from rdkit import Chem
from rdkit.Chem import rdDetermineBonds

def remove_explicit_hydrogens_from_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    mol = Chem.RemoveHs(mol)
    return Chem.MolToSmiles(mol)

def kekulize_smiles(smiles):
    mol = Chem.MolFromSmiles(smiles, sanitize=False)
    Chem.Kekulize(mol)
    return Chem.MolToSmiles(mol, kekuleSmiles=True)

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--input', '-i', required=True, type=str, help='Path of the input xyz file')
parser.add_argument('--charge', '-c', default=argparse.SUPPRESS, type=int, help='Net charge of the molecule')
parser.add_argument('--kekule', '-k', action='store_false', help='Once this option is specified, kekule SMILES instead of the default canonical SMILES will be outputted')
args = parser.parse_args()

rawmol = Chem.MolFromXYZFile(args.input)
mol = Chem.Mol(rawmol)
if 'charge' in args:
    c = args.charge
else:
    c = 0
rdDetermineBonds.DetermineBonds(mol, charge=c)
rawsmiles = Chem.MolToSmiles(mol)
smiles = remove_explicit_hydrogens_from_smiles(rawsmiles)
if args.kekule:
    print(smiles)
else:
    print(kekulize_smiles(smiles))
