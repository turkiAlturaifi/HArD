from rdkit import Chem
from rdkit.Chem import AllChem, Draw
import csv
from natsort import natsorted, ns, index_natsorted
import subprocess
import glob
import pandas as pd
### Section 1 ###
def process_products_1(products, heteroaryl, reaction_type, all_mod_mols_data, start_position_id):
    for product_tuple in products:
        mod_mol = Chem.RemoveHs(product_tuple[0])
        mod_smiles = Chem.MolToSmiles(mod_mol, canonical=True)
        if mod_smiles not in all_mod_mols_data:
            mol_id = f"{heteroaryl['id']}_{start_position_id}_{reaction_type}"
            all_mod_mols_data[mod_smiles] = {
                'ID': mol_id,
                'SMILES': mod_smiles,
                'parent_name': heteroaryl['name']
            }
            start_position_id += 1
    return start_position_id

def process_heteroaryls(input_csv, output_csv):
    heteroaryls = []
    with open(input_csv, mode='r') as csvfile:
        reader = csv.DictReader(csvfile, delimiter=',')
        headers = reader.fieldnames
        for row in reader:
            heteroaryls.append({'id': row['ID'], 'name': row['parent_name'], 'smiles': row['SMILES']})

    all_mod_mols_data = {}

    formylation_reaction = AllChem.ReactionFromSmarts('[cH:1]>>[c:1]C(=O)O')
    carboxylation_reaction = AllChem.ReactionFromSmarts('[cH:1]>>[c:1]C(=O)[O-]')

    for heteroaryl in heteroaryls:
        parent_mol = Chem.MolFromSmiles(heteroaryl['smiles'])
        parent_smiles = heteroaryl['smiles']
        parent_id = f"{heteroaryl['id']}_1_a"
        all_mod_mols_data[parent_smiles] = {
            'ID': parent_id,
            'SMILES': parent_smiles,
            'parent_name': heteroaryl['name']
        }

        position_id_a = 1
        position_id_b = 1

        formylated_products = formylation_reaction.RunReactants((parent_mol,))
        process_products_1(formylated_products, heteroaryl, 'b', all_mod_mols_data, position_id_a)
        carboxylated_products = carboxylation_reaction.RunReactants((parent_mol,))
        process_products_1(carboxylated_products, heteroaryl, 'c', all_mod_mols_data, position_id_b)

    with open(output_csv, 'w', newline='') as csvfile:
        fieldnames = ['ID', 'SMILES', 'parent_name']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for data in natsorted(all_mod_mols_data.values(), key=lambda x: x['ID']):
            writer.writerow(data)

input_csv = "batch1.csv" 
output_csv = "step_1.csv"
process_heteroaryls(input_csv, output_csv)
number_of_entries = sum(1 for row in csv.reader(open(input_csv, 'r'))) - 1
print(f"Number of SMILES in '{input_csv}': {number_of_entries}")

### Section 2 ###
def process_products(products, heteroaryl, reaction_name, all_mod_mols_data, fg_position_map, seen_smiles):
    position_id = fg_position_map.get(reaction_name, 1)
    removed_duplicates = 0
    id_prefix, id_suffix = heteroaryl['id'].split('_')[0:2], heteroaryl['id'].split('_')[-1]

    for product_set in products:
        for product in product_set:
            mol_smiles = Chem.MolToSmiles(product)
            # Skip this product if it's a duplicate
            if mol_smiles in seen_smiles:
                removed_duplicates += 1
                continue
            
            seen_smiles.add(mol_smiles)
            # mod_id = f"{heteroaryl['id'].split('_')[0:2]}_{position_id}_{reaction_name}_{heteroaryl['id'].split('_')[2]}"
            mod_id = f"{'_'.join(id_prefix)}_{position_id}_{reaction_name}_{id_suffix}"
            all_mod_mols_data[mol_smiles] = {'ID': mod_id, 'SMILES': mol_smiles, 'parent_name': heteroaryl['name']}
            position_id += 1

    fg_position_map[reaction_name] = position_id
    return removed_duplicates

# Define a list of SMARTS patterns for the reactions
reactions_smarts = [
    ('1', '[cH:1]>>[c:1]N(C)C'),
    ('2', '[cH:1]>>[c:1]N'),
    ('3', '[cH:1]>>[c:1]O'),
    ('4', '[cH:1]>>[c:1]OC'),
    ('5', '[cH:1]>>[c:1]C'),
    ('6', '[cH:1]>>[c:1][Si](C)(C)C'),
    ('7', '[cH:1]>>[c:1]F'),
    ('8', '[cH:1]>>[c:1]Cl'),
    ('9', '[cH:1]>>[c:1]Br'),
    ('10', '[cH:1]>>[c:1]C(=O)C'),
    ('11', '[cH:1]>>[c:1]C#N'),
    ('12', '[cH:1]>>[c:1][N+](=O)[O-]'),
    ('13', '[cH:1]>>[cH:1]'),
]

heteroaryls = []
with open('step_1.csv', mode='r') as csvfile:
    reader = csv.DictReader(csvfile, delimiter=',')
    for row in reader:
        heteroaryls.append({'id': row['ID'], 'name': row['parent_name'], 'smiles': row['SMILES']})

all_mod_mols_data = {} 
seen_smiles = set() 
fg_position_map = {} 
total_removed_duplicates = 0 

for heteroaryl in heteroaryls:
    parent_mol = Chem.MolFromSmiles(heteroaryl['smiles'])
    fg_position_map = {}  # Reset for each new heteroaryl parent
    for reaction_name, smarts in reactions_smarts:
        reaction = AllChem.ReactionFromSmarts(smarts)
        products = reaction.RunReactants((parent_mol,))
        removed_duplicates = process_products(products, heteroaryl, reaction_name, all_mod_mols_data, fg_position_map, seen_smiles)
        total_removed_duplicates += removed_duplicates

csv_file_path = 'step_2.csv'
with open(csv_file_path, 'w', newline='') as csvfile:
    fieldnames = ['ID', 'SMILES', 'parent_name']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    
    writer.writeheader()
    # sort naturally based on ID
    sorted_data = natsorted(all_mod_mols_data.values(), key=lambda x: x['ID'], alg=ns.IGNORECASE)
    for data in sorted_data:
        writer.writerow(data)

### Section 3 ###

def step_2_cleanup(filename):
  """
  This function performs three sed operations on the specified file:

  1. Deletes lines containing "_b".
  2. Deletes lines containing "_c".
  3. Replaces all occurrences of "_a" with an empty string.

  Args:
      filename (str): The name of the file to modify.
  """

  # Define commands
  command1 = ["sed", "-i", "", "/_b/d", filename]  # Delete lines with "_b"
  command2 = ["sed", "-i", "", "/_c/d", filename]  # Delete lines with "_c"
  command3 = ["sed", "-i", "", "s/_a//g", filename]  # Replace all "_a" with ""

  subprocess.run(command1, check=True)
  subprocess.run(command2, check=True)
  subprocess.run(command3, check=True)

step_2_cleanup("step_2.csv") 

### section 4 ###
input_csv = "step_2.csv" 
output_csv = "step_3.csv"
process_heteroaryls(input_csv, output_csv) # Process the heteroaryls again

### section 5 ###
def process_csv(input_file, output_file):
  """
  Processes a CSV file by modifying the first column ID and writes the result to a new CSV file.

  Args:
      input_file (str): Path to the input CSV file.
      output_file (str): Path to the output CSV file.
  """
  with open(input_file, 'r', newline='') as infile, open(output_file, 'w', newline='') as outfile:
    reader = csv.DictReader(infile, delimiter=',')
    writer = csv.DictWriter(outfile, fieldnames=reader.fieldnames)
    writer.writeheader()
    for row in reader:
      id_parts = row['ID'].split('_')
      new_id = f"{id_parts[0]}_{id_parts[2]}_{id_parts[3]}_{id_parts[4]}_{id_parts[5]}"
      row['ID'] = new_id
      writer.writerow(row)

input_csv = "step_3.csv"
output_csv = "final_temp.csv"
process_csv(input_csv, output_csv)

files_to_delete = glob.glob("step*.csv")
command4 = ["rm"] + files_to_delete
subprocess.run(command4, check=True)

print(f"CSV processing completed. IDs and SMILES are written to {output_csv}.")
with open('final_temp.csv', 'r') as csvfile:
  reader = csv.DictReader(csvfile, delimiter=',')
  count = sum(1 for row in reader if row['ID'].endswith('_a'))
  print(f"Total number of unique substituted heteroaryls: {count}")
number_of_entries = sum(1 for row in csv.reader(open('final_temp.csv', 'r'))) - 1
print(f"Total number of calculations with CO2H and CO2-: {number_of_entries}")
print(f"Total number of duplicate SMILES removed: {total_removed_duplicates}")

def replace_and_sort(file_path):
    df = pd.read_csv(file_path)
    df['ID'] = df['ID'].str.replace('_13_', '_0_')
    df = df.iloc[index_natsorted(df['ID'])]
    output_file_path = 'final.csv'
    df.to_csv(output_file_path, index=False)
    print(f"File saved to {output_file_path}")
file_path = 'final_temp.csv'
replace_and_sort(file_path)

files_to_delete = glob.glob("final_te*.csv")
command4 = ["rm"] + files_to_delete
subprocess.run(command4, check=True)

def execute_gen_opt():
    # if you want to run gen-gallery.py
    user_request = input("run gen-opt.py to generate xyz and gaussian input files? (yes/no): ").lower()
    if user_request == 'yes':
        subprocess.run(['python', 'gen-opt.py'])
        print("xyz and com files have been generated")
    elif user_request == 'no':
        print("no gallery is generated")
    else:
        print("Invalid input. gen-opt.py execution was not requested.")
execute_gen_opt()

