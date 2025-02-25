import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdmolops
import cclib
import numpy as np
import glob
import os

# Read the CSV file
csv_file_path = 'ring_atoms_final.csv'
csv_data = pd.read_csv(csv_file_path)

b_type_csv_path = 'ring_atoms_final.csv'
b_type_csv_data = pd.read_csv(b_type_csv_path)

# Function to parse a log file and extract Mulliken and Natural charges
def parse_log_file(logfile_path):
    try:
        logfile = cclib.io.ccread(logfile_path)
        mulliken_charges = logfile.atomcharges['mulliken'] if 'mulliken' in logfile.atomcharges else None
        natural_charges = logfile.atomcharges['natural'] if 'natural' in logfile.atomcharges else None
        return mulliken_charges, natural_charges
    except Exception as e:
        print(f"Error parsing {logfile_path}: {e}")
        return None, None

def calculate_hetaryl_chrg(log_file, natural_charges, mulliken_charges, row, b_data):
    file_type = log_file.split('_')[-2]
    hetaryl_chrg_natural = hetaryl_chrg_mulliken = 0
    # print(log_file)

    if file_type == 'a':
        z_axis_atom = int(row['z_axis_atom']) - 1
        hetaryl_chrg_natural = -1 * natural_charges[z_axis_atom]
        hetaryl_chrg_mulliken = -1 * mulliken_charges[z_axis_atom]
    elif file_type == 'b':
        excluded_atoms = [int(atom) - 1 for atom in eval(b_data['CO2H_atoms'].iloc[0])]
        # print(excluded_atoms)
        sum_natural = sum(natural_charges[i] for i in excluded_atoms)
        sum_mulliken = sum(mulliken_charges[i] for i in excluded_atoms)
        # base_charge = 1
        hetaryl_chrg_natural = -1 * sum_natural
        hetaryl_chrg_mulliken = -1 * sum_mulliken
    elif file_type == 'c':
        excluded_atoms = [int(atom) - 1 for atom in eval(row['excluded_atoms'])[1:]]
        sum_natural = sum(natural_charges[i] for i in excluded_atoms)
        sum_mulliken = sum(mulliken_charges[i] for i in excluded_atoms)
        base_charge = -1
        hetaryl_chrg_natural = base_charge - sum_natural
        hetaryl_chrg_mulliken = base_charge - sum_mulliken
    # print('hetaryl_chrg_natural',hetaryl_chrg_natural)
    # print(' hetaryl_chrg_mulliken',hetaryl_chrg_natural)

    return hetaryl_chrg_natural, hetaryl_chrg_mulliken

# Function to calculate properties for a given row and log file charges
def calculate_properties(row, mulliken_charges, natural_charges):
    try:
        ring_atoms = [int(atom) - 1 for atom in eval(row['ring_atoms'])]
        atom_with_fg1 = int(row['Atom with FG1']) - 1 if not pd.isna(row['Atom with FG1']) else None
        atom_with_fg2 = int(row['Atom with FG2']) - 1 if not pd.isna(row['Atom with FG2']) else None
        metal_index = int(row['metal_index']) - 1 if not pd.isna(row['metal_index']) else None

        # NPA charges
        npa_ipso_acid = natural_charges[metal_index] if metal_index is not None else np.nan
        if atom_with_fg2 is not None:
            npa_ipso_fg = natural_charges[atom_with_fg2]
        elif atom_with_fg1 is not None:
            npa_ipso_fg = natural_charges[atom_with_fg1]
        else:
            npa_ipso_fg = np.nan
        npa_ring_ave = sum(natural_charges[i] for i in ring_atoms) / len(ring_atoms)

        # Mulliken charges
        mulliken_ipso_acid = mulliken_charges[metal_index] if metal_index is not None else np.nan
        if atom_with_fg2 is not None:
            mulliken_ipso_fg = mulliken_charges[atom_with_fg2]
        elif atom_with_fg1 is not None:
            mulliken_ipso_fg = mulliken_charges[atom_with_fg1]
        else:
            mulliken_ipso_fg = np.nan
        mulliken_ring_ave = sum(mulliken_charges[i] for i in ring_atoms) / len(ring_atoms)

        properties = {
        'File': row['File'],
        'npa_ipso_acid': natural_charges[int(row['metal_index']) - 1] if pd.notna(row['metal_index']) else np.nan,
        'npa_ipso_fg': natural_charges[int(row['Atom with FG1']) - 1] if pd.notna(row['Atom with FG1']) else np.nan,
        'npa_ring_ave': np.mean([natural_charges[int(atom) - 1] for atom in eval(row['ring_atoms'])]),
        'mulliken_ipso_acid': mulliken_charges[int(row['metal_index']) - 1] if pd.notna(row['metal_index']) else np.nan,
        'mulliken_ipso_fg': mulliken_charges[int(row['Atom with FG1']) - 1] if pd.notna(row['Atom with FG1']) else np.nan,
        'mulliken_ring_ave': np.mean([mulliken_charges[int(atom) - 1] for atom in eval(row['ring_atoms'])])
        }

          # Calculate hetaryl charges and update the properties dictionary
        hetaryl_chrg_natural, hetaryl_chrg_mulliken = calculate_hetaryl_chrg(log_file, natural_charges, mulliken_charges, row, b_data)
        properties.update({
            'hetaryl_chrg_natural': hetaryl_chrg_natural,
            'hetaryl_chrg_mulliken': hetaryl_chrg_mulliken
        })

        return properties
    
    except Exception as e:
        print(f"Error w {log_file}")
        print(f"Error calculating properties for {row['File']}: {e}")
        return None

# Get all log files in the current directory ending with b_SP.log or c_SP.log
# log_files = glob.glob('*[bc]_SP.log')
log_files = glob.glob('*_SP.log')
all_results = []

# Process each log file
for log_file in log_files:
    mulliken_charges, natural_charges = parse_log_file(log_file)
    if mulliken_charges is None or natural_charges is None:
        print(f"Skipping {log_file} due to parsing error.")
        continue
    csv_file_name = log_file.replace('.log', '.xyz')
    filtered_data = csv_data[csv_data['File'] == csv_file_name]
    b_data = b_type_csv_data[b_type_csv_data['File'] == csv_file_name] if 'b' in log_file else None
    try:
        results = filtered_data.apply(lambda row: calculate_properties(row, mulliken_charges, natural_charges), axis=1)
        all_results.extend([res for res in results.tolist() if res is not None])
    except Exception as e:
        print(f"Error processing results for {log_file}: {e}")

# Convert all results to a DataFrame and save to CSV
results_df = pd.DataFrame(all_results)
results_df.to_csv('charges-npa.csv', index=False)

# clean up
def modify_csv():
    file_path = 'charges-npa.csv'
    df = pd.read_csv(file_path)
    df.drop(columns=['npa_ipso_fg', 'mulliken_ipso_fg'], inplace=True)
    df.rename(columns={
        'npa_ipso_acid': 'npa_ipso',
        'npa_ring_ave': 'npa_ring_average',
        'mulliken_ipso_acid': 'mulliken_ipso',
        'mulliken_ring_ave': 'mulliken_ring_average',
        'hetaryl_chrg_natural': 'npa_hetaryl',
        'hetaryl_chrg_mulliken': 'mulliken_hetaryl'
    }, inplace=True)
    df.to_csv(file_path, index=False)

if __name__ == "__main__":
    modify_csv()