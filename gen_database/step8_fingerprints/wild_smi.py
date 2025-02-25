import pandas as pd
import numpy as np
import natsort
import pandas as pd
from rdkit import Chem
from rdkit.Chem import RWMol

combined_df = pd.read_csv('final.csv')
combined_df['BaseFilename'] = combined_df['ID'].str[:-2]
base_filenames = combined_df['BaseFilename'].unique()
new_columns = []
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
output_path = 'wild_smi.csv'
combined_df_final.to_csv(output_path, index=False)

def replace_carboxylic_acid_with_wildcard(smiles):
    if pd.isna(smiles):
        return None 
    try:
        carboxylic_acid_smarts = "C(=O)O"
        mol = Chem.MolFromSmiles(smiles)
        query = Chem.MolFromSmarts(carboxylic_acid_smarts)
        match = mol.GetSubstructMatch(query)
        
        if match:
            emol = RWMol(mol)
            atom_idx = match[0]
            emol.ReplaceAtom(atom_idx, Chem.Atom('*'))
            sorted_indices = sorted(match[1:], reverse=True)
            for idx in sorted_indices:
                emol.RemoveAtom(idx)
            new_mol = emol.GetMol()
            modified_smiles = Chem.MolToSmiles(new_mol)
            return modified_smiles
        else:
            return smiles 
    except Exception as e:
        return None 


df = pd.read_csv("wild_smi.csv")
df['Modified_SMILES_b'] = df['SMILES_b'].apply(replace_carboxylic_acid_with_wildcard)
missing_ids = df[df['Modified_SMILES_b'].isnull()]['ID']  
converted_count = df['Modified_SMILES_b'].notnull().sum()
skipped_count = df['Modified_SMILES_b'].isnull().sum()
df.to_csv("wild_smi.csv", index=False)

# Print summary
print(f"Number of converted SMILES: {converted_count}")
print(f"Number of skipped (missing) SMILES: {skipped_count}")
if skipped_count > 0:
    print("IDs with missing SMILES:")
    print(missing_ids.to_string(index=False))

print('wild_smi.csv')