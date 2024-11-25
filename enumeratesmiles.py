'''
This script generates multiple possible SMILES strings from a given SMILES string.

Perform ten trials for generating SMILES strings for mesocarb:
  python enumeratesmiles.py -s "C[C@@H]([N+]1=NOC(/N=C([O-])/NC2=CC=CC=C2)=C1)CC3=CC=CC=C3" -n 10
'''

import argparse
import numpy as np
from rdkit import Chem

parser = argparse.ArgumentParser(usage=__doc__)
parser.add_argument('--string', '-s', required=True, type=str, help='SMILES string')
parser.add_argument('--number', '-n', required=True, type=int, help='Number of trials for generating different SMILES strings')
args = parser.parse_args()

mol = Chem.MolFromSmiles(args.string)
smileslist = Chem.MolToRandomSmilesVect(mol, args.number, isomericSmiles=True, kekuleSmiles=True, randomSeed=np.random.randint(0, 1000))
print(*smileslist, sep="\n")
