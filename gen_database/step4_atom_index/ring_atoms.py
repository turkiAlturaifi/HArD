import os
import csv
import subprocess
import glob
from openbabel import openbabel as ob
from rdkit import Chem
from rdkit.Chem import AllChem, Draw, rdDetermineBonds, rdmolops
import pandas as pd

print('-----section 1----')

def convert_xyz_to_smiles(xyz_file):
    """Convert xyz file to SMILES using Open Babel."""
    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats("xyz", "smi")

    mol = ob.OBMol()
    if obConversion.ReadFile(mol, xyz_file):
        return obConversion.WriteString(mol).strip()
    return None

def convert_log_to_smiles(log_file):
    """Convert Gaussian log file to SMILES using Open Babel."""
    obConversion = ob.OBConversion()
    obConversion.SetInAndOutFormats("g16", "smi")

    mol = ob.OBMol()
    if obConversion.ReadFile(mol, log_file):
        return obConversion.WriteString(mol).strip()
    return None

def find_ring_atom_numbers(molecule):
    """Identify atom numbers in rings and count non-ring atoms attached to ring atoms."""

    ssr = Chem.GetSymmSSSR(molecule) 
    ring_atoms = []
    attached_non_ring_atoms = {}

    for ring in ssr:
       ring_atoms.extend([atom_idx + 1 for atom_idx in list(ring)])
            
    for atom in molecule.GetAtoms():
        atom_idx = atom.GetIdx() + 1
        if atom.IsInRing():
            non_ring_count = 0
            for neighbor in atom.GetNeighbors():
                neighbor_idx = neighbor.GetIdx() + 1
                if neighbor_idx not in ring_atoms and neighbor.GetAtomicNum() != 1:
                    non_ring_count += 1
                    # print('neighbor '+str(neighbor_idx)+': '+str(non_ring_count))
            attached_non_ring_atoms[atom_idx] = non_ring_count

    return ring_atoms, attached_non_ring_atoms

def check_attached_to_carboxyl_group(atom_idx, molecule, ring_atoms):
    """Check if the atom is connected to a non-ring atom that is further connected to two oxygens."""
    atom = molecule.GetAtomWithIdx(atom_idx - 1)  # Convert to 0-based index
    for neighbor in atom.GetNeighbors():
        neighbor_idx = neighbor.GetIdx() + 1
        # Ensure the neighbor is a non-ring carbon and not a hydrogen atom
        if neighbor_idx not in ring_atoms and neighbor.GetAtomicNum() == 6:
            oxygen_count = sum(1 for n in neighbor.GetNeighbors() if n.GetAtomicNum() == 8)  # Count oxygen atoms
            if oxygen_count == 2:
                return True
    return False

def find_oxygen_atoms(z_axis_atom, molecule):
    """Find the two oxygen atoms connected to the z_axis_atom."""
    oxygen_atoms = []
    for neighbor in molecule.GetAtomWithIdx(z_axis_atom - 1).GetNeighbors():
        if neighbor.GetAtomicNum() == 8:  # Oxygen atom
            oxygen_atoms.append(neighbor.GetIdx() + 1)
    return oxygen_atoms if len(oxygen_atoms) == 2 else []

def main():
#    log_files = [file for file in os.listdir() if file.endswith('.log')]
    xyz_files = [file for file in os.listdir() if file.endswith('.xyz')]
    output_data = []
    mols = []
    labels = []

    for xyz_file in xyz_files:
        
        try:
# read xyz into a mol object
            # print('reading '+xyz_file)
            mol = AllChem.MolFromXYZFile(xyz_file)
            rdDetermineBonds.DetermineConnectivity(mol)
            mols.append(mol)
            labels.append(xyz_file)

            smiles = convert_xyz_to_smiles(xyz_file)
            
            ring_atom_numbers, attached_non_ring_atoms = find_ring_atom_numbers(mol)
            non_ring_attached_atoms = [atom for atom, count in attached_non_ring_atoms.items() if count > 0]
            row = [xyz_file, ring_atom_numbers]
            if len(non_ring_attached_atoms) >= 1:
                row.append(non_ring_attached_atoms[0])
            if len(non_ring_attached_atoms) >= 2:
                row.append(non_ring_attached_atoms[1])
            output_data.append(row)
        except Exception as e:
            print(f"Error processing {xyz_file}: {e}")

    with open('ring_atoms7.csv', 'w', newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerow(['File', 'ring_atoms', 'Atom with FG1', 'Atom with FG2', 'Attached to Carboxyl Group', 'metal_index', 'xz_atom', 'z_axis_atom', 'excluded_atoms'])
        for row in output_data:
            xyz_file, ring_atoms, *fg_atoms = row
            if len(fg_atoms) < 2:
                fg_atoms.extend([None] * (2 - len(fg_atoms)))
            atom_with_fg1, atom_with_fg2 = fg_atoms
            # mol = Chem.MolFromSmiles(convert_xyz_to_smiles(xyz_file))
            mol = AllChem.MolFromXYZFile(xyz_file)
            rdDetermineBonds.DetermineConnectivity(mol)
            attached_to_carboxyl = []
            metal_index = None
            z_axis_atom = None
            for atom_idx in fg_atoms:
                if atom_idx:
                    try:
                        is_attached = check_attached_to_carboxyl_group(atom_idx, mol, ring_atoms)
                        attached_to_carboxyl.append(is_attached)
                        if is_attached:
                            metal_index = atom_idx
                            for neighbor in mol.GetAtomWithIdx(metal_index - 1).GetNeighbors():
                                neighbor_idx = neighbor.GetIdx() + 1
                                if neighbor_idx not in ring_atoms and neighbor.GetAtomicNum() != 1:
                                    z_axis_atom = neighbor_idx
                                    break
                    except Exception as e:
                        print(f"Error checking carboxyl group attachment for atom {atom_idx} in {xyz_file}: {e}")
            xz_atom = None
            if metal_index:
                try:
                    idx = ring_atoms.index(metal_index)
                    if idx + 1 < len(ring_atoms):
                        xz_atom = ring_atoms[idx + 1]
                    else:
                        xz_atom = ring_atoms[idx - 1]  
                except ValueError:
                    xz_atom = None
            oxygen_atoms = find_oxygen_atoms(z_axis_atom, mol) if z_axis_atom else []
            excluded_atoms = [metal_index, z_axis_atom] + oxygen_atoms
            writer.writerow([xyz_file, ring_atoms, atom_with_fg1, atom_with_fg2, attached_to_carboxyl, metal_index, xz_atom, z_axis_atom, excluded_atoms])
if __name__ == '__main__':
    main()

print('-----section 2----')
import pandas as pd

ring_df = pd.read_csv('ring_atoms7.csv')
combined_df = pd.read_csv('ipso_a.csv')

ring_df['ring_atoms'] = ring_df['ring_atoms'].apply(eval)

mask = ring_df['File'].str.endswith('_a.xyz')
filtered_ring_df = ring_df[mask]

def update_row(row, combined_df):
    file_id = row['File'][:-6]
    combined_row = combined_df[combined_df['ID'] == file_id]
    if not combined_row.empty:
        # formylation_site = int(combined_row['Formylation_Site'].values[0])
        formylation_site = combined_row['Formylation_Site'].values[0]
        if pd.notna(formylation_site):  # Check if not NaN
            formylation_site = int(formylation_site)  # Safely convert to integer
            # z_axis_atom = combined_row['z_axis_atom'].values[0] 
            z_axis_atom = combined_row['z_axis_atom'].values[0]
            if pd.notna(z_axis_atom):  # Check if not NaN
                z_axis_atom = int(z_axis_atom)  # Safely convert to integer
            row['Atom with FG2'] = formylation_site
            row['metal_index'] = formylation_site
            row['z_axis_atom'] = z_axis_atom 

            row['Attached to Carboxyl Group'] = [False, True]
            row['excluded_atoms'] = ""

            ring_atoms = row['ring_atoms']
            if formylation_site in ring_atoms:
                index = ring_atoms.index(formylation_site)
                next_index = (index + 1) % len(ring_atoms)
                row['xz_atom'] = ring_atoms[next_index]
    return row
updated_ring_df = filtered_ring_df.apply(lambda row: update_row(row, combined_df), axis=1)
final_df = pd.concat([ring_df[~mask], updated_ring_df], ignore_index=True)
final_df.to_csv('ring_atoms7_updated.csv', index=False)
df = pd.read_csv('ring_atoms7_updated.csv', header=None)
def replace_and_copy(row):
    if '_a_SP.log' in row[0]:
        base_row = df[df[0] == row[0].replace('_a_SP.log', '_a.log')].iloc[0]
        row[1:] = base_row[1:]

df.apply(replace_and_copy, axis=1)
df.to_csv('ring_atoms7_updated.csv', index=False, header=False)


def process_files(ring_file, smiles_file, output_file):
    ring_df = pd.read_csv(ring_file)
    smiles_df = pd.read_csv(smiles_file)
    ring_df['ID'] = ring_df['File'].str.extract(r'(\d+_\d+_\d+_\d+)')[0]
    smiles_df.rename(columns={'SMILES_b': 'SMILES'}, inplace=True)
    combined_df = pd.merge(ring_df, smiles_df[['ID', 'SMILES']], on='ID', how='left')
    combined_df['excluded_atoms'] = combined_df['excluded_atoms'].apply(lambda x: eval(x) if isinstance(x, str) else x)
    processed_data = []
    for index, row in combined_df.iterrows():
        if 'b' in row['File']:
            if pd.notna(row['SMILES']) and isinstance(row['SMILES'], str):
                excluded_atoms = row['excluded_atoms'][-3:]
                mol = Chem.MolFromSmiles(row['SMILES'])
                mol = Chem.AddHs(mol)
                co2h_atoms = excluded_atoms.copy()
                for atom_idx in excluded_atoms[-2:]:
                    if atom_idx is not None:  # Check if atom_idx is not None
                        atom = mol.GetAtomWithIdx(atom_idx - 1)  # zero-based index in RDKit
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetAtomicNum() == 1:  # hydrogen atom
                                co2h_atoms.append(neighbor.GetIdx() + 1)  # convert to one-based index
                                break
                row['CO2H_atoms'] = co2h_atoms
            else:
                row['CO2H_atoms'] = 'Invalid SMILES'
                print(f"Invalid or missing SMILES for ID {row['ID']}")
        else:
            row['CO2H_atoms'] = []
        processed_data.append(row)
    new_df = pd.DataFrame(processed_data)
    new_df.to_csv(output_file, index=False)

process_files('ring_atoms7_updated.csv', 'ipso_a.csv', 'ring_atoms_final.csv')

files_to_delete = ['ring_atoms7_updated.csv','ring_atoms7.csv']
command4 = ["rm"] + files_to_delete
subprocess.run(command4, check=True)