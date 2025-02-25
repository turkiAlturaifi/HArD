
import natsort
import csv
import pandas as pd
import subprocess
import numpy as np

def extract_blocks(file_path, output_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    blocks = []
    current_block = []
    capture = False
    for i, line in enumerate(lines):
        if '*****' in line:
            if capture:
                blocks.append(current_block)
                current_block = []
            capture = not capture
            if capture:
                current_block.append(lines[i-1].strip())  #the line before '*****'
        if capture:
            current_block.append(line.strip())

    if current_block:
        blocks.append(current_block)

    with open(output_path, 'w') as output_file:
        for block in blocks:
            for line in block:
                output_file.write(line + '\n')
            output_file.write('\n')

if __name__ == "__main__":
    input_file_path = 'gvibes_all.out' 
    output_file_path = 'out.out'
    extract_blocks(input_file_path, output_file_path)

#  remove 'o' and '*****'
def process_blocks_to_csv(blocks, header):
    processed_data = [header]  

    for block in blocks:
        lines = block.split('\n')
        for line in lines:
            if line.startswith('o  '):
                cleaned_line = line[2:].strip().split()
                processed_data.append(' '.join(cleaned_line))  
            elif not line.startswith('*') and line != header:
                cleaned_line = line.strip().split()
                processed_data.append(' '.join(cleaned_line)) 
    
    return '\n'.join(processed_data)

def convert_to_csv_format(processed_data):
    lines = processed_data.split('\n')
    csv_data = []
    for line in lines:
        csv_data.append(line.split())
    return csv_data

def sort_csv_data_naturally(csv_data):
    header = csv_data[0]
    rows = csv_data[1:]
    sorted_rows = natsort.natsorted(rows, key=lambda x: x[0])
    return [header] + sorted_rows

file_path = 'out.out'
with open(file_path, 'r') as file:
    content = file.read()

blocks = content.strip().split('\n\n')
headers_type1 = "Structure                                           E        ZPE             H        T.S     T.qh-S          G(T)       qh-G(T)"
headers_type2 = "Structure                                       E_SPC             E        ZPE         H_SPC        T.S     T.qh-S      G(T)_SPC   qh-G(T)_SPC"

blocks_type1 = []
blocks_type2 = []

for block in blocks:
    if headers_type1 in block:
        blocks_type1.append(block)
    elif headers_type2 in block:
        blocks_type2.append(block)

processed_csv_content_type1 = process_blocks_to_csv(blocks_type1, headers_type1)
processed_csv_content_type2 = process_blocks_to_csv(blocks_type2, headers_type2)
csv_data_type1 = convert_to_csv_format(processed_csv_content_type1)
csv_data_type2 = convert_to_csv_format(processed_csv_content_type2)
sorted_csv_data_type1 = sort_csv_data_naturally(csv_data_type1)
sorted_csv_data_type2 = sort_csv_data_naturally(csv_data_type2)
output_sorted_csv_path_type1 = 'gas.csv'
output_sorted_csv_path_type2 = 'sp.csv'
with open(output_sorted_csv_path_type1, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(sorted_csv_data_type1)
with open(output_sorted_csv_path_type2, 'w', newline='') as file:
    writer = csv.writer(file)
    writer.writerows(sorted_csv_data_type2)


gas_df = pd.read_csv('gas.csv')
sp_df = pd.read_csv('sp.csv')
gas_columns = ['Structure', 'E', 'H', 'G(T)', 'qh-G(T)']
sp_columns = ['Structure', 'E_SPC', 'H_SPC', 'G(T)_SPC', 'qh-G(T)_SPC']
gas_selected = gas_df[gas_columns]
sp_selected = sp_df[sp_columns]


merged_df = pd.merge(gas_selected, sp_selected, on='Structure')
merged_df_sorted = merged_df.sort_values(by='Structure', key=natsort.natsort_keygen())
output_path = 'energies_master.csv'
merged_df_sorted.to_csv(output_path, index=False)

def modify_csv():
    input_file_path = 'energies_master.csv'
    output_file_path = 'energies.csv'
    df = pd.read_csv(input_file_path)
    df.drop(columns=['E', 'H', "G(T)", "qh-G(T)", 'E_SPC', 'H_SPC', 'qh-G(T)_SPC'], inplace=True)
    df.to_csv(output_file_path, index=False)

if __name__ == "__main__":
    modify_csv()


# combine _b and _c with _a
combined_df = pd.read_csv('energies.csv')
combined_df['BaseFilename'] = combined_df['Structure'].str[:-2]
base_filenames = combined_df['BaseFilename'].unique()
new_columns = []

# Triple the amount of headers, adding _a, _b, _c to each column except 'Structure'
for column in combined_df.columns:
    if column not in ['Structure', 'BaseFilename']:
        new_columns.extend([f'{column}_a', f'{column}_b', f'{column}_c'])

new_columns = ['Structure'] + new_columns
combined_rows = []
for base_filename in base_filenames:
    rows = combined_df[combined_df['BaseFilename'] == base_filename]
    combined_row = {'Structure': base_filename}
    for suffix in ['a', 'b', 'c']:
        row = rows[rows['Structure'].str.endswith(suffix)]
        if not row.empty:
            row = row.iloc[0]
            for column in combined_df.columns:
                if column not in ['Structure', 'BaseFilename']:
                    combined_row[f'{column}_{suffix}'] = row[column]
        else:
            for column in combined_df.columns:
                if column not in ['Structure', 'BaseFilename']:
                    combined_row[f'{column}_{suffix}'] = np.nan

    combined_rows.append(combined_row)

combined_df_final = pd.DataFrame(combined_rows, columns=new_columns)
combined_df_final = combined_df_final.sort_values(by='Structure', key=natsort.natsort_keygen())
output_path = 'energies.csv'
combined_df_final.to_csv(output_path, index=False)

files_to_delete = ['sp.csv','gas.csv']
command4 = ["rm"] + files_to_delete
subprocess.run(command4, check=True)

print(f"energies.csv, run sigma-het.py")