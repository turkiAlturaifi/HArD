import csv
import os
from rdkit import Chem
from rdkit.Chem import AllChem
from openbabel import openbabel as ob
import numpy as np

def read_xyz(filename):
    mol = ob.OBMol()
    conv = ob.OBConversion()
    conv.SetInFormat("xyz")
    conv.ReadFile(mol, filename)
    return mol

def get_connectivity_matrix(mol):
    n_atoms = mol.NumAtoms()
    matrix = np.zeros((n_atoms, n_atoms), dtype=int)
    for bond in ob.OBMolBondIter(mol):
        i, j = bond.GetBeginAtomIdx() - 1, bond.GetEndAtomIdx() - 1
        matrix[i, j] = matrix[j, i] = 1
    return matrix

def compare_connectivity(smiles, xyz_file):
    mol_smiles = Chem.MolFromSmiles(smiles)
    mol_smiles = Chem.AddHs(mol_smiles)
    AllChem.EmbedMolecule(mol_smiles)
    conv = ob.OBConversion()
    conv.SetInAndOutFormats("sdf", "sdf")
    mol_ob = ob.OBMol()
    conv.ReadString(mol_ob, Chem.MolToMolBlock(mol_smiles))
    matrix_smiles = get_connectivity_matrix(mol_ob)
    mol_xyz = read_xyz(xyz_file)
    matrix_xyz = get_connectivity_matrix(mol_xyz)
    return np.array_equal(matrix_smiles, matrix_xyz)

# Read CSV file
with open('final.csv', 'r') as f:
    reader = csv.DictReader(f)
    for row in reader:
        id = row['ID']
        smiles = row['SMILES']
        xyz_file = f"{id}.xyz"
        if os.path.exists(xyz_file):
            connectivity_preserved = compare_connectivity(smiles, xyz_file)
            if not connectivity_preserved:
                if id.endswith('1_3_1_c'):
                    print(f"ID: {id}, Connectivity not preserved (this can be ignored)")
                else:
                    print(f"ID: {id}, Connectivity not preserved")
        else:
            print(f"ID: {id}, XYZ file not found")