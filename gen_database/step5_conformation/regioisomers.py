import os
import pandas as pd
import shutil

def check_missing_log_files_and_duplicates(csv_file, directory):
    df = pd.read_csv(csv_file)
    ids = df['ID'].unique()
    expected_files = {}
    for id in ids:
        for suffix in ['a_SP.log', 'b_SP.log', 'c_SP.log']:
            expected_files[f"{id}_{suffix}"] = suffix.split('_')[0]
    existing_files = [file for file in os.listdir(directory) if file.endswith('.log')]
    missing_files = [file for file in expected_files if file not in existing_files]
    reports = [] 
    if missing_files:
        print("Missing log files:")
        for file in missing_files:
            # print(file)
            id = file.split('_')[0] + '_' + file.split('_')[1] + '_' + file.split('_')[2] + '_' + file.split('_')[3]
            smile_type = expected_files[file]  # e.g., SMILES_a, SMILES_b, or SMILES_c
            smile_column = f"SMILES_{smile_type}"
            missing_smile = df.loc[df['ID'] == id, smile_column].values[0]
            
            duplicates = df[df[smile_column] == missing_smile]['ID'].unique()
            if len(duplicates) > 1:
                duplicate_files = [f"{dup_id}_{smile_type}_SP.log" for dup_id in duplicates if f"{dup_id}_{smile_type}_SP.log" in existing_files]
                if duplicate_files:
                    for original_file in duplicate_files:
                        reports.append({'Missing': file, 'Original': original_file})
                else:
                    print(file)
                    print("No existing log files for found duplicates.")
            else:
                print(file)
                print(f"No duplicate SMILES for {smile_column}.")
    else:
        print("No log files are missing.")
    
    if reports:
        report_df = pd.DataFrame(reports)
        report_df.to_csv('regioisomers.csv', index=False)
        print("Report saved to regioisomers.csv.")
    else:
        print("No report to save, no missing files found.")

check_missing_log_files_and_duplicates('ipso_a.csv', '.')

def copy_original_to_missing(csv_file):
    if not os.path.exists(csv_file):
        print("CSV file does not exist.")
        return

    df = pd.read_csv(csv_file)
    if df.empty:
        print("No data found in the CSV file.")
        return

    regioisomers_dir = os.path.join('.', 'regioisomers') # creat a 'regioisomers' folder
    if not os.path.exists(regioisomers_dir):
        os.makedirs(regioisomers_dir)
        print(f"Directory '{regioisomers_dir}' created.")

    for index, row in df.iterrows():
        original_file = row['Original']
        missing_file = row['Missing']
        source_path = os.path.join('.', original_file)
        target_path = os.path.join(regioisomers_dir, missing_file)
        if not os.path.exists(source_path):
            print(f"Original file {original_file} does not exist.")
            continue

        shutil.copy(source_path, target_path)
        print(f"Copied {original_file} to {missing_file} in '{regioisomers_dir}' directory.")

csv_file_path = 'regioisomers.csv'
copy_original_to_missing(csv_file_path)
