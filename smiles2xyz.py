'''
This script generates xyz files from a SMILES string. This script requires
RDKit which should be at least version 2020.03. Cite the following paper 
if you use this script in any published work.
S. Wang, J. Witek, G. A. Landrum, S. Riniker, J. Chem. Inf. Model. 2020, 60(4), 2044–2058

Generate one conformer for a molecule with the following command:
  python smiles2xyz.py -s "C[C@@H]([N+]1=NOC(/N=C([O-])/NC2=CC=CC=C2)=C1)CC3=CC=CC=C3" -c mesocarb
Generate ten conformers for a molecule with the following command:
  python smiles2xyz.py -s "C[C@@H]([N+]1=NOC(/N=C([O-])/NC2=CC=CC=C2)=C1)CC3=CC=CC=C3" -c mesocarb -n 10
'''

import sys
import argparse
from rdkit import Chem
from rdkit.Chem import AllChem

def write_xyz_file(fragment, fragment_name):
    number_of_atoms = fragment.GetNumAtoms()
    symbols = [a.GetSymbol() for a in fragment.GetAtoms()]
    fNames = []
    for i,conf in enumerate(fragment.GetConformers()):
        file_name = fragment_name+"_"+str(i)+".xyz"
        fNames.append(file_name)
        with open(file_name, "w") as file:
            file.write(str(number_of_atoms)+"\n")
            file.write("\n")
            for atom,symbol in enumerate(symbols):
                p = conf.GetAtomPosition(atom)
                line = " ".join((symbol,str(p.x),str(p.y),str(p.z),"\n"))
                file.write(line)
    return fNames

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--string', '-s', required=True, type=str, help='SMILES string')
parser.add_argument('--compound', '-c', required=True, type=str, help='Compound name, which will be used as the prefix for xyz file')
parser.add_argument('--number', '-n', default=argparse.SUPPRESS, type=int, help='Number of conformers for a compound')
args = parser.parse_args()

mol = Chem.MolFromSmiles(args.string)
mol_h = Chem.AddHs(mol)

params = Chem.rdDistGeom.srETKDGv3()
params.randomSeed = 0xf00d
params.clearConfs = True
if 'number' in args:
    if args.number > 0:
        number_of_conformers = args.number
    else:
        print("The argument for number of conformers in the input cannot be used.")
        sys.exit(1)
else:
    number_of_conformers = 1
cids = AllChem.EmbedMultipleConfs(mol_h, number_of_conformers, params)
fn = write_xyz_file(mol_h, args.compound)
