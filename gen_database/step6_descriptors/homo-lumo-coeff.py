import os
import re
import csv
import pandas as pd
def extract_coefficients(log_file, metal_index):
    homo_coeff = 0.0
    lumo_coeff = 0.0
    homo_orbital_type = None
    lumo_orbital_type = None
    start_reading = False
    last_alpha_occ_line = None
    
    try:
        with open(log_file, 'r') as file:
            content = file.readlines()
    except FileNotFoundError:
        print(f"File not found: {log_file}")
        return homo_coeff, homo_orbital_type, lumo_coeff, lumo_orbital_type
    
    # First pass: find the last "Alpha occ" and the first "Alpha vir"
    for line in content:
        if "Atomic contributions to Alpha molecular orbitals:" in line:
            start_reading = True
        if start_reading:
            if "Alpha occ" in line:
                last_alpha_occ_line = line
            if "Alpha vir" in line:
                matches = re.findall(rf'C{metal_index}-([sp])=(-?[0-9\.]+)', line)
                if matches:
                    lumo_orbital_type, lumo_coeff = matches[0]  # First occurrence
                break
    
    # Process the last "Alpha occ" line if available
    if last_alpha_occ_line:
        matches = re.findall(rf'C{metal_index}-([sp])=(-?[0-9\.]+)', last_alpha_occ_line)
        if matches:
            homo_orbital_type, homo_coeff = matches[-1]  # Last occurrence

    return homo_coeff, homo_orbital_type, lumo_coeff, lumo_orbital_type

def main():
    try:
        df = pd.read_csv('ring_atoms_final.csv')
    except FileNotFoundError:
        print("CSV file not found")
        return
    results = []
    for index, row in df.iterrows():
        file_name = row['File']
        if 'SP' not in file_name:
            continue
        base_name = '_'.join(file_name.split('_')[:6])
        log_file_name = f"{base_name.split('.')[0]}_homo-lumo.log"
        if not os.path.exists(log_file_name):
            print(f"File not found: {log_file_name}")
            results.append({
                'File': log_file_name,
                'homo_coeff': '',
                'homo_orbital_type': '',
                'lumo_coeff': '',
                'lumo_orbital_type': '',
                'atom_number': row['metal_index'] if not pd.isna(row['metal_index']) else ''
            })
            continue
        metal_index = row['metal_index']
        if pd.isna(metal_index):
            print(f"Warning: metal_index is NaN for {log_file_name}")
            results.append({
                'File': log_file_name,
                'homo_coeff': '',
                'homo_orbital_type': '',
                'lumo_coeff': '',
                'lumo_orbital_type': '',
                'atom_number': ''
            })
            continue
        homo_coeff, homo_orbital_type, lumo_coeff, lumo_orbital_type = extract_coefficients(log_file_name, int(metal_index))
        if homo_coeff == 0.0 and lumo_coeff == 0.0:
            print(f"Warning: Missing values in {log_file_name} - HOMO: {homo_coeff}, LUMO: {lumo_coeff}")
        results.append({
            'File': log_file_name,
            'homo_coeff': homo_coeff if homo_coeff != 0.0 else '',
            'homo_orbital_type': homo_orbital_type if homo_orbital_type is not None else '',
            'lumo_coeff': lumo_coeff if lumo_coeff != 0.0 else '',
            'lumo_orbital_type': lumo_orbital_type if lumo_orbital_type is not None else '',
            'atom_number': int(metal_index) if metal_index is not None else ''
        })
    with open('homo-lumo-coeff1.csv', 'w', newline='') as csvfile:
        fieldnames = ['File', 'homo_coeff', 'homo_orbital_type', 'lumo_coeff', 'lumo_orbital_type', 'atom_number']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        for row in results:
            writer.writerow(row)
    
    print(f"Results written to homo-lumo-coeff.csv. Number of entries: {len(results)}")

if __name__ == '__main__':
    main()
