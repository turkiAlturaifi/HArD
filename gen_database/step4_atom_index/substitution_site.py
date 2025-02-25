import pandas as pd
import numpy as np
import natsort
from rdkit import Chem
from rdkit.Chem import AllChem
import os
# reformat csv
combined_df = pd.read_csv('final.csv')
combined_df['BaseFilename'] = combined_df['ID'].str[:-2]
base_filenames = combined_df['BaseFilename'].unique()
new_columns = []
# Triple the amount of headers, adding _a, _b, _c to each column except 'ID'
for column in combined_df.columns:
    if column not in ['ID', 'BaseFilename']:
        new_columns.extend([f'{column}_a', f'{column}_b', f'{column}_c'])

new_columns = ['ID'] + new_columns
combined_rows = []

for base_filename in base_filenames:
    rows = combined_df[combined_df['BaseFilename'] == base_filename]
    combined_row = {'ID': base_filename}
    for suffix in ['a', 'b', 'c']:
        row = rows[rows['ID'].str.endswith(suffix)]
        if not row.empty:
            row = row.iloc[0]
            for column in combined_df.columns:
                if column not in ['ID', 'BaseFilename']:
                    combined_row[f'{column}_{suffix}'] = row[column]
        else:
            for column in combined_df.columns:
                if column not in ['ID', 'BaseFilename']:
                    combined_row[f'{column}_{suffix}'] = np.nan

    combined_rows.append(combined_row)
combined_df_final = pd.DataFrame(combined_rows, columns=new_columns)
combined_df_final = combined_df_final.sort_values(by='ID', key=natsort.natsort_keygen())
#output_path = 'final_combined.csv'
#combined_df_final.to_csv(output_path, index=False)

#print(f"Combined final CSV file created: {output_path}")


# get the atom number on the parent that correspond to the atom number on the acid/base pair
formylation_reaction = AllChem.ReactionFromSmarts('[cH:1]>>[c:1]C(=O)O')

def find_formylation_site(smiles_a, smiles_b):
    mol_a = Chem.MolFromSmiles(smiles_a)
    mol_b = Chem.MolFromSmiles(smiles_b)
    smiles_b = Chem.MolToSmiles(mol_b)
    for atom in mol_a.GetAtoms():
        if atom.GetAtomicNum() == 6 and atom.GetTotalNumHs() > 0:
            test_mol = Chem.Mol(mol_a)
            test_mol.UpdatePropertyCache(strict=False)
            Chem.SanitizeMol(test_mol, catchErrors=True)
            idx = atom.GetIdx()
            test_mol.GetAtomWithIdx(idx).SetAtomMapNum(1)
            ps = formylation_reaction.RunReactants((test_mol,))
            for p in ps:
                for product in p:
                    Chem.SanitizeMol(product, catchErrors=True)
                    product_smiles = Chem.MolToSmiles(product)
                    if product_smiles == smiles_b:
                        return idx + 1
    return None

df = combined_df_final
df['SMILES_a'] = df['SMILES_a'].fillna(method='ffill')
df['Formylation_Site'] = df.apply(lambda row: find_formylation_site(row['SMILES_a'], row['SMILES_b']) if pd.notna(row['SMILES_b']) else None, axis=1)
# df.to_csv('final_combined_modified.csv', index=False)
# print("Modified CSV file saved as 'final_combined_modified.csv'.")

# Now you get the atom number of H atoms
def read_xyz(xyz_file):
    """Read an XYZ file and return element symbols and coordinates."""
    with open(xyz_file, 'r') as file:
        lines = file.readlines()[2:]  # Skip the first two metadata lines
    symbols = []
    coordinates = []
    for line in lines:
        parts = line.split()
        symbols.append(parts[0])
        coordinates.append([float(x) for x in parts[1:4]])
    return symbols, np.array(coordinates)

def find_bonded_hydrogen(symbols, coords, carbon_atom_index):
    """Finds the index of the hydrogen atom bonded to a specified carbon atom."""
    carbon_coords = coords[carbon_atom_index]
    bonded_hydrogen_index = None
    min_distance = float('inf')  # Initialize with a large number

    for i, (symbol, coord) in enumerate(zip(symbols, coords)):
        if symbol == 'H':
            distance = np.linalg.norm(coord - carbon_coords)
            if distance < 1.2 and distance < min_distance:  # C-H bond typical distance
                min_distance = distance
                bonded_hydrogen_index = i

    return bonded_hydrogen_index + 1 if bonded_hydrogen_index is not None else None

# df = pd.read_csv('final_combined_modified.csv')
df['z_axis_atom'] = None
xyz_directory = '.'
for index, row in df.iterrows():
    if pd.notna(row['Formylation_Site']):  # Check if Formylation_Site is not NaN
        xyz_file = os.path.join(xyz_directory, f"{row['ID']}_a.xyz")
        if not os.path.exists(xyz_file):
            # Handle duplicates in SMILES_a
            duplicate_ids = df[df['SMILES_a'] == row['SMILES_a']]['ID']
            for dup_id in duplicate_ids:
                xyz_file = os.path.join(xyz_directory, f"{dup_id}_a.xyz")
                if os.path.exists(xyz_file):
                    break
        
        if os.path.exists(xyz_file):
            symbols, coords = read_xyz(xyz_file)
            hydrogen_index = find_bonded_hydrogen(symbols, coords, int(row['Formylation_Site']) - 1)
            df.at[index, 'z_axis_atom'] = hydrogen_index
    else:
        print(f"Skipping row {index+1}: Missing 'Formylation_Site'.")

df.to_csv('ipso_a.csv', index=False)
