import pandas as pd
import os
import glob
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdmolops

# def process_files(ring_file, smiles_file, output_file):
#     # Load data
#     ring_df = pd.read_csv(ring_file)
#     smiles_df = pd.read_csv(smiles_file)
    
#     # Check initial data
#     # print("Initial data from ring file:", ring_df.head())
#     # print("Initial data from smiles file:", smiles_df.head())

#     # Correct ID extraction from File
#     ring_df['ID'] = ring_df['File'].str.extract(r'(\d+_\d+_\d+_\d+)')[0]
    
#     # Rename 'SMILES_b' column to 'SMILES' for uniformity before merging
#     smiles_df.rename(columns={'SMILES_b': 'SMILES'}, inplace=True)

#     # Merge the DataFrames on 'ID'
#     combined_df = pd.merge(ring_df, smiles_df[['ID', 'SMILES']], on='ID', how='left')
    
#     # Check merged data
#     # print("Merged data:", combined_df.head())

#     # Process only entries where 'File' contains 'b'
#     # Safely convert 'excluded_atoms' to list, checking its type first
#     combined_df['excluded_atoms'] = combined_df['excluded_atoms'].apply(lambda x: eval(x) if isinstance(x, str) else x)
#     processed_data = []
#     for index, row in combined_df.iterrows():
#         if 'b' in row['File']:
#             # Check if SMILES is not NaN and is a string
#             if pd.notna(row['SMILES']) and isinstance(row['SMILES'], str):
#                 # Extract the last three numbers from 'excluded_atoms'
#                 excluded_atoms = row['excluded_atoms'][-3:]
#                 # Load the molecule
#                 mol = Chem.MolFromSmiles(row['SMILES'])
#                 mol = Chem.AddHs(mol)
#                 # Find the hydrogen atom bonded to either of the last two atoms in excluded_atoms
#                 co2h_atoms = excluded_atoms.copy()
#                 for atom_idx in excluded_atoms[-2:]:
#                     atom = mol.GetAtomWithIdx(atom_idx - 1)  # zero-based index in RDKit
#                     for neighbor in atom.GetNeighbors():
#                         if neighbor.GetAtomicNum() == 1:  # hydrogen atom
#                             co2h_atoms.append(neighbor.GetIdx() + 1)  # convert to one-based index
#                             break
#                 row['CO2H_atoms'] = co2h_atoms
#             else:
#                 row['CO2H_atoms'] = 'Invalid SMILES'
#                 print(f"Invalid or missing SMILES for ID {row['ID']}")
#         else:
#             row['CO2H_atoms'] = []
#         processed_data.append(row)

#     # Convert list back to DataFrame
#     new_df = pd.DataFrame(processed_data)

#     # Write to CSV
#     new_df.to_csv(output_file, index=False)

# # Usage
# process_files('../sp/ring_atoms7_updated.csv', '../sp/final_combined_modified1_with_hydrogens.csv', 'ring_atoms7_b.csv')


def extract_charges(log_file_path):
    try:
        with open(log_file_path, 'r') as file:
            lines = file.readlines()
        start = False
        charges = []
        for line in lines:
            if "Hirshfeld charges, spin densities" in line:
                start = True
                continue
            if start:
                if 'Tot' in line or "Q-H" in line:  # Skip processing if it's the ending line or a header line
                    if 'Tot' in line:  # Stop processing altogether if it's the ending line
                        break
                    continue
                parts = line.split()
                if len(parts) > 6:  # Ensure there are enough parts to avoid IndexError
                    q_h = float(parts[2])  # Q-H value
                    q_cm5 = float(parts[-1])  # Q-CM5 value
                    charges.append((q_h, q_cm5))
                else:
                    print(f"Unexpected format in line: {line}")
    except Exception as e:
        print(f"Failed to parse file {log_file_path}: {e}")
        return []
    return charges

def calculate_hetaryl_chrg(log_file, charges, row, b_data):
    file_type = log_file.split('_')[-2]
    hetaryl_chrg_qh = 0  # Charge on the Q-H axis
    hetaryl_chrg_qcm5 = 0  # Charge on the Q-CM5 axis

    if file_type == 'a':
        z_axis_atom = int(row['z_axis_atom']) - 1
        if z_axis_atom < len(charges):
            hetaryl_chrg_qh, hetaryl_chrg_qcm5 = -1* charges[z_axis_atom][0],-1* charges[z_axis_atom][-1]
    elif file_type == 'b':
        excluded_atoms = [int(atom) - 1 for atom in eval(b_data['CO2H_atoms'])]
        hetaryl_chrg_qh = -1 * sum(charges[i][0] for i in excluded_atoms if i < len(charges))
        hetaryl_chrg_qcm5 = -1 * sum(charges[i][1] for i in excluded_atoms if i < len(charges))
    elif file_type == 'c':
        excluded_atoms = [int(atom) - 1 for atom in eval(row['excluded_atoms'])[1:]]
        hetaryl_chrg_qh = -1 - sum(charges[i][0] for i in excluded_atoms if i < len(charges))
        hetaryl_chrg_qcm5 = -1 -sum(charges[i][1] for i in excluded_atoms if i < len(charges))

    return hetaryl_chrg_qh, hetaryl_chrg_qcm5

def calculate_properties(row, charges, b_data=None):
    try:
        ring_atoms = [int(atom) - 1 for atom in eval(row['ring_atoms'])]
        atom_with_fg1 = int(row['Atom with FG1']) - 1 if not pd.isna(row['Atom with FG1']) else None
        atom_with_fg2 = int(row['Atom with FG2']) - 1 if not pd.isna(row['Atom with FG2']) else None
        metal_index = int(row['metal_index']) - 1 if not pd.isna(row['metal_index']) else None

        # Ensure we choose the FG atom that is not the same as the metal index
        fg_index = None
        if atom_with_fg1 is not None and atom_with_fg1 != metal_index:
            fg_index = atom_with_fg1
        elif atom_with_fg2 is not None and atom_with_fg2 != metal_index:
            fg_index = atom_with_fg2

        hetaryl_chrg_qh, hetaryl_chrg_qcm5 = calculate_hetaryl_chrg(row['File'], charges, row, b_data)

        properties = {
            'File': row['File'],
            'Q-H_ipso_acid': charges[metal_index][0] if metal_index is not None and metal_index < len(charges) else None,
            'Q-H_ipso_fg': charges[fg_index][0] if fg_index is not None and fg_index < len(charges) else None,
            'Q-H_ring_ave': np.mean([charges[i][0] for i in ring_atoms if i < len(charges)]),
            'Q-CM5_ipso_acid': charges[metal_index][1] if metal_index is not None and metal_index < len(charges) else None,
            'Q-CM5_ipso_fg': charges[fg_index][1] if fg_index is not None and fg_index < len(charges) else None,
            'Q-CM5_ring_ave': np.mean([charges[i][1] for i in ring_atoms if i < len(charges)]),
            'hetaryl_chrg_qh': hetaryl_chrg_qh,
            'hetaryl_chrg_qcm5': hetaryl_chrg_qcm5
        }

        return properties
    except Exception as e:
        print(f"Error calculating properties for {row['File']}: {e}")
        return None

# Load CSV data
csv_data = pd.read_csv('ring_atoms_final.csv')
b_type_csv_data = pd.read_csv('ring_atoms_final.csv')

# Process each log file in the directory
log_files = glob.glob('*_SP_hirsh.log')
all_results = []

for log_file in log_files:
    # print(f"Processing {log_file}")
    charges = extract_charges(log_file)
    if not charges:
        print(f"Skipping {log_file} due to parsing error.")
        continue

    filtered_data = csv_data[csv_data['File'].apply(lambda x: x.split('_SP')[0]) == os.path.basename(log_file).split('_SP')[0]]
    if filtered_data.empty:
        print(f"No matching entries for {log_file} in CSV.")
        continue

    b_data = b_type_csv_data[b_type_csv_data['File'].apply(lambda x: x.split('_SP')[0]) == os.path.basename(log_file).split('_SP')[0]]
    b_data = b_data.iloc[0] if not b_data.empty else None

    try:
        for index, row in filtered_data.iterrows():
            # print(f"Calculating properties for row: {index}")
            properties = calculate_properties(row, charges, b_data)
            if properties:
                all_results.append(properties)
    except Exception as e:
        print(f"Error processing results for {log_file}: {e}")

# Convert all results to a DataFrame and save to CSV
results_df = pd.DataFrame(all_results)
results_df.to_csv('charges_hirsh.csv', index=False)


# clean up
def modify_csv():
    file_path = 'charges_hirsh.csv'
    df = pd.read_csv(file_path)
    df.drop(columns=['Q-H_ipso_fg', 'Q-CM5_ipso_fg'], inplace=True)
    df.rename(columns={
        'Q-H_ipso_acid': 'hirshfeld_ipso_ipso',
        'Q-H_ring_ave': 'hirshfeld_ring_average',
        'Q-CM5_ipso_acid': 'cm5_ipso',
        'Q-CM5_ring_ave': 'cm5_ring_average',
        'hetaryl_chrg_qh': 'hirshfeld_hetaryl',
        'hetaryl_chrg_qcm5': 'cm5_hetaryl'
    }, inplace=True)
    df.to_csv(file_path, index=False)

if __name__ == "__main__":
    modify_csv()

print("charges_hirsh.csv")