import pandas as pd
from rdkit import Chem
from rdkit.Chem import HybridizationType, Descriptors, rdMolDescriptors
import csv

# charge
def has_charge(smiles):
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return False
        return sum(atom.GetFormalCharge() for atom in mol.GetAtoms()) != 0
    except:
        return False

# radicals
def has_radicals(smiles):
    try:
        mol = Chem.MolFromSmiles(str(smiles))
        if mol is None:
            return False
        return any(atom.GetNumRadicalElectrons() > 0 for atom in mol.GetAtoms())
    except:
        return False

# isotopes and others (checked manually)
def has_deteriorated_structure(smiles):
    isotopes = ['2H', '15N', '13C', '18O','3H','14C','I','Ge','nan','P','Si','B','[NH]1','[NH]C2']
    return any(isotope in smiles for isotope in isotopes)

data = pd.read_excel('reaxys_search.xlsx')
data['SMILES'] = data['SMILES'].astype(str)
charged_mask = data['SMILES'].apply(has_charge)
radical_mask = data['SMILES'].apply(has_radicals)
deteriorated_structure_mask = data['SMILES'].apply(has_deteriorated_structure)
print("Number of charged entries excluded:", charged_mask.sum())
print("Number of radicals entries excluded:", radical_mask.sum())
print("Number of deteriorated structure entries excluded:", deteriorated_structure_mask.sum())

filtered_data = data[~(charged_mask | radical_mask | deteriorated_structure_mask)]
smiles_data = filtered_data.to_dict(orient='records')

# sort by number of rings, number of heavy atoms, then molecular weight
def get_molecule_properties(mol):
    if not mol:
        return (0, 0, 0, 0.0)
    num_rings = rdMolDescriptors.CalcNumRings(mol)
    heavy_atoms = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() > 1)
    mol_weight = Descriptors.MolWt(mol)
    return (num_rings, heavy_atoms, mol_weight)

for data in smiles_data:
    mol = Chem.MolFromSmiles(data['SMILES'])
    data['num_rings'], data['heavy_atoms'], data['mol_weight'] = get_molecule_properties(mol)
sorted_smiles_data = sorted(smiles_data, key=lambda x: (x['num_rings'],  x['heavy_atoms'], x['mol_weight']))

# cleanup
del_col = {'mol_weight', 'Structure: Image'}
fieldnames = ['SMILES', 'num_rings', 'heavy_atoms']
additional_fieldnames = [key for key in smiles_data[0].keys() if key not in fieldnames and key not in del_col]
fieldnames.extend(additional_fieldnames)
with open('reaxys_processed.csv', 'w', newline='') as csvfile:
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for data in sorted_smiles_data:
        data_to_write = {key: value for key, value in data.items() if key not in del_col}
        writer.writerow(data_to_write)

# condensed csv 
with open('batch1.csv', 'w', newline='') as csvfile:
    fieldnames = ['ID', 'SMILES', 'parent_name']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
    writer.writeheader()
    for i, data in enumerate(sorted_smiles_data, start=1):
        parent_name = data.get('Chemical Name', '').split(';')[0]  # split by ';' and take the first name
        writer.writerow({'ID': i, 'SMILES': data['SMILES'], 'parent_name': parent_name})

print("generated reaxys_processed.csv and condnsed version batch1.csv")
