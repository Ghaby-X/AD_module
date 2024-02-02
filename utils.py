import pandas as pd
import numpy as np
import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem

def compute_morganfps(smile_array):
    mols = get_mols_from_array(smile_array)
    fps = [AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048) for mol in mols]
    return fps


def get_mols_from_array(smile_array):
    mols = [Chem.MolFromSmiles(smi) for smi in smile_array]
    return mols

