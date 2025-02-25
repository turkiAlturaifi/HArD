import os
import glob
import pandas as pd
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdDetermineBonds

def xyz_to_mol_and_count_rings(xyz_path):
    mol = AllChem.MolFromXYZFile(xyz_path)
    if mol is None:
        print(f"Warning: Could not convert {xyz_path} to a valid molecule.")
        return None
    rdDetermineBonds.DetermineConnectivity(mol)
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    ssr = Chem.GetSymmSSSR(mol)
    num_rings = len(ssr)
    return num_rings

def get_base_id(file_name):
    return '_'.join(file_name.split('_')[:4])

def get_suffix(file_name):
    return file_name.split('_')[4]

output_csv_file = 'ring_counts_with_warnings.csv'
xyz_files = glob.glob("*.xyz")
results = []

for xyz_file in xyz_files:
    print(xyz_file)
    base_id = get_base_id(os.path.basename(xyz_file))
    suffix = get_suffix(os.path.basename(xyz_file))
    num_rings = xyz_to_mol_and_count_rings(xyz_file)

    if num_rings is not None:
        first_field = int(base_id.split('_')[0])
        expected_rings = 1 if 1 <= first_field <= 33 else 2   
        warning = ''
        if num_rings != expected_rings:
            warning = 'bad geometry'
        results.append({
            'filename': f"{base_id}_{suffix}.xyz",
            'id': base_id,
            'expected_number_of_rings': expected_rings,
            'computed_number_of_rings': num_rings,
            'warning': warning
        })
    else:
        print(f"Warning: Could not count rings for {xyz_file}.")

results_df = pd.DataFrame(results)
results_df.to_csv(output_csv_file, index=False)
print(f"Results saved to {output_csv_file}")
