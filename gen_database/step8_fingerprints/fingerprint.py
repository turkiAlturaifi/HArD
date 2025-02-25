import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
import csv

csv_file_path = 'final.csv'
df = pd.read_csv(csv_file_path)
output_csv = 'fingerprints.csv'


with open(output_csv, mode='w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)
    csv_writer.writerow([
        'ID', 'SMILES', 'parent_name', 'number_of_rings', 'number_of_heteroatoms',
        'number_of_hydrogen_bond_donors', 'number_of_hydrogen_bond_acceptors', 'topological_polar_surface_area_tpsa',
        'molecular_weight', 'logp', 'rotatable_bonds_count',
        'number_of_heavy_atoms', 'number_of_carbon_atoms', 'number_of_nitrogen_atoms',
        'number_of_oxygen_atoms', 'number_of_sulfur_atoms',
        'fraction_of_aromatic_bonds', 'fraction_of_sp3_hybridized_carbons', 'fraction_of_sp2_hybridized_carbons'
    ])

    for index, row in df.iterrows():
        try:
            filename = row['ID']
            smiles = row['SMILES']
            parent=row['parent_name']
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                raise ValueError(f"Could not parse SMILES from {filename}")

            #  descriptors
            num_rings = Descriptors.RingCount(mol)
            num_heteroatoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() not in {1, 6})  # Not H or C
            num_hbond_donors = Descriptors.NumHDonors(mol)
            num_hbond_acceptors = Descriptors.NumHAcceptors(mol)
            tpsa = Descriptors.TPSA(mol)
            mol_weight = Descriptors.MolWt(mol)
            logp = Descriptors.MolLogP(mol)
            rotatable_bonds = Descriptors.NumRotatableBonds(mol)
            #num_aliphatic_chains = Descriptors.NumAliphaticCarbocycles(mol)
            num_heavy_atoms = Descriptors.HeavyAtomCount(mol)
            num_carbon_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
            num_nitrogen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7)
            num_oxygen_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
            num_sulfur_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 16)
            #num_phosphorus_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 15)
            fraction_aromatic_bonds = sum(1 for bond in mol.GetBonds() if bond.GetIsAromatic()) / mol.GetNumBonds()
            fraction_sp3_carbons = Descriptors.FractionCSP3(mol)
            fraction_sp2_carbons = 1 - fraction_sp3_carbons
            
            # csv
            csv_writer.writerow([
                filename, smiles, parent,num_rings, num_heteroatoms,
                num_hbond_donors, num_hbond_acceptors, tpsa,
                mol_weight, logp, rotatable_bonds,
                num_heavy_atoms, num_carbon_atoms, num_nitrogen_atoms,
                num_oxygen_atoms, num_sulfur_atoms,
                fraction_aromatic_bonds, fraction_sp3_carbons, fraction_sp2_carbons
            ])

        except Exception as e:
            print(f"Error processing {filename}: {e}")

output_csv
