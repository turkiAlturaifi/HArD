import os
import csv
import pandas as pd
import subprocess
import glob
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolTransforms, rdDetermineBonds
df = pd.read_csv('ring_atoms_final.csv')
output_df = pd.DataFrame(columns=['File', 'ring_atoms', 'atom_1', 'atom_2', 'atom_3', 'atom_4'])

rows = []

for index, row in df.iterrows():
    if '_c' in row['File']:
        file = row['File']
        ring_atoms = row['ring_atoms']
        atom_1 = row['xz_atom']
        atom_2 = row['metal_index']
        atom_3 = row['z_axis_atom']
        atom_4 = eval(row['excluded_atoms'])[2] if len(eval(row['excluded_atoms'])) > 1 else ''
        
        rows.append({
            'File': file,
            'ring_atoms': ring_atoms,
            'atom_1': atom_1,
            'atom_2': atom_2,
            'atom_3': atom_3,
            'atom_4': atom_4
        })

output_df = pd.concat([output_df, pd.DataFrame(rows)], ignore_index=True)
output_df.to_csv('output.csv', index=False)

def read_xyz_as_mol(file_path):
    with open(file_path, 'r') as f:
        xyz_data = f.read()
    mol = Chem.MolFromXYZBlock(xyz_data)
    return mol

def get_dihedral_deg(mol, indices):
    return Chem.rdMolTransforms.GetDihedralDeg(mol.GetConformer(), *indices)

df_atoms = pd.read_csv('output.csv')
output_df = pd.DataFrame(columns=['File', 'dihedral_angle', 'is_coplanar'])
for index, row in df_atoms.iterrows():
    file_name = row['File']
    file_path = os.path.join(os.getcwd(), file_name)
    if os.path.exists(file_path):
        print(file_name)
        mol = read_xyz_as_mol(file_path)
        if mol:
            try:
                #atom_indices = [int(row['atom_' + str(i)]) for i in range(1, 5) if pd.notna(row['atom_' + str(i)])]
                atom_indices = [int(row['atom_' + str(i)]) - 1 for i in range(1, 5) if pd.notna(row['atom_' + str(i)])]
                if len(atom_indices) == 4:  # Check if all four indices are present
                    dihedral = get_dihedral_deg(mol, atom_indices)
                    is_coplanar = abs(dihedral) < 10 or abs(abs(dihedral) - 180) < 10
                    temp_df = pd.DataFrame({'File': [row['File']], 'dihedral_angle': [dihedral], 'is_coplanar': [is_coplanar]})
                else:
                    raise ValueError("Missing atom indices")
            except ValueError as e:
                print(f"Warning: {e} for {row['File']}")
                temp_df = pd.DataFrame({'File': [row['File']], 'dihedral_angle': [''], 'is_coplanar': ['']})
        else:
            print(f"Warning: Could not load molecule from {file_name}")
            temp_df = pd.DataFrame({'File': [row['File']], 'dihedral_angle': [''], 'is_coplanar': ['']})
    else:
        print(f"Warning: File {file_name} not found")
        temp_df = pd.DataFrame({'File': [row['File']], 'dihedral_angle': [''], 'is_coplanar': ['']})
    
    output_df = pd.concat([output_df, temp_df], ignore_index=True)
output_df.to_csv('dihedral_c_temp.csv', index=False)

df1 = pd.read_csv('dihedral_c_temp.csv')
df2 = pd.read_csv('output.csv')
combined_df = pd.merge(df1, df2, on='File')
combined_df.to_csv('combined_output.csv', index=False)

def set_dihedral(mol, idx1, idx2, idx3, idx4, angle_deg):
    """Set the dihedral angle to the specified value."""
    conf = mol.GetConformer()
    rdMolTransforms.SetDihedralDeg(conf, idx1-1, idx2-1, idx3-1, idx4-1, angle_deg)
    return mol

def save_xyz(mol, filename):
    """Save the molecule to an XYZ file."""
    with open(filename, 'w') as f:
        conf = mol.GetConformer()
        num_atoms = mol.GetNumAtoms()
        f.write(f"{num_atoms}\n\n")
        for i in range(num_atoms):
            atom = mol.GetAtomWithIdx(i)
            pos = conf.GetAtomPosition(i)
            f.write(f"{atom.GetSymbol()} {pos.x:.6f} {pos.y:.6f} {pos.z:.6f}\n")

def process_file(input_file, ring_atoms, idx1, idx2, idx3, idx4, angle_deg):
    """Process an individual XYZ file and adjust the dihedral angle."""
    try:
        mol = AllChem.MolFromXYZFile(input_file)
        if mol is None:
            return False, "Failed to load molecule"

        rdDetermineBonds.DetermineConnectivity(mol)
        Chem.SanitizeMol(mol)
        mol.UpdatePropertyCache()

        mol = set_dihedral(mol, idx1, idx2, idx3, idx4, angle_deg)
        output_file = input_file.replace('.xyz', '_conf2.xyz')
        save_xyz(mol, output_file)
        return True, None
    except Exception as e:
        return False, str(e)

def main():
    csv_input = "combined_output.csv"
    result_csv = "processing_results.csv"
    angle_deg = 90  # or -90

    with open(csv_input, mode='r') as infile:
        reader = csv.DictReader(infile)
        rows = list(reader)

    results = []
    for row in rows:
        file_id = row['File']
        input_file = file_id
        
        if not os.path.isfile(input_file):
            results.append({'File': file_id, 'status': 'skipped', 'reason': 'File not found'})
            continue

        try:
            ring_atoms = eval(row['ring_atoms'])
            idx1 = int(float(row['atom_1']))
            idx2 = int(float(row['atom_2']))
            idx3 = int(float(row['atom_3']))
            idx4 = int(float(row['atom_4']))
        except (ValueError, SyntaxError) as e:
            results.append({'File': file_id, 'status': 'skipped', 'reason': f'Invalid data: {str(e)}'})
            continue

        success, reason = process_file(input_file, ring_atoms, idx1, idx2, idx3, idx4, angle_deg)
        if success:
            results.append({'File': file_id, 'status': 'processed', 'reason': None})
        else:
            results.append({'File': file_id, 'status': 'skipped', 'reason': reason})

    with open(result_csv, mode='w', newline='') as outfile:
        writer = csv.DictWriter(outfile, fieldnames=['File', 'status', 'reason'])
        writer.writeheader()
        writer.writerows(results)

    print(f"Processing complete. Results saved to dihedral_c.csv")

if __name__ == "__main__":
    main()


df1 = pd.read_csv('combined_output.csv')
df2 = pd.read_csv('processing_results.csv')
combined_df = pd.merge(df1, df2, on='File')
combined_df.to_csv('dihedral_c.csv', index=False)
files_to_delete = ['output.csv','dihedral_c_temp.csv','combined_output.csv','processing_results.csv']
command4 = ["rm"] + files_to_delete
subprocess.run(command4, check=True)