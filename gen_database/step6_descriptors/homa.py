import pandas as pd
import numpy as np
import os
from cclib.io import ccread

# Define bond parameters
bond_parameters = {
    ('C', 'C'): (1.3880, 257.7000),
    ('C', 'N'): (1.3340, 93.5200),
    ('C', 'O'): (1.2650, 157.3800),
    ('C', 'S'): (1.6770, 94.0900),
    ('N', 'N'): (1.3090, 130.3300),
    ('N', 'O'): (1.2480, 57.2100),
    ('N', 'S'): (1.61, 71.875),
}

def calculate_homa(file, ring_atoms):
    print(f"Processing file: {file}")
    log_file = file.replace('.xyz', '.log')
    data = ccread(log_file)
    # data = ccread(file)
    coordinates = data.atomcoords[-1]
    atom_types = data.atomnos
    element_symbols = {1: 'H', 6: 'C', 7: 'N', 8: 'O', 16: 'S'}
    
    bond_lengths = []
    valid_bond_lengths = []
    seen_bonds = set()
    n = len(ring_atoms)
    expected_bonds = n - 1 if n > 6 else n

    for i in range(n):
        for j in range(i + 1, n):
            if ring_atoms[i] != ring_atoms[j]:
                sorted_indices = tuple(sorted((ring_atoms[i], ring_atoms[j])))
                if sorted_indices not in seen_bonds:
                    dist = np.linalg.norm(coordinates[ring_atoms[i]-1] - coordinates[ring_atoms[j]-1])
                    bond_lengths.append(dist)
                    if dist < 1.85:
                        valid_bond_lengths.append((dist, sorted_indices[0], sorted_indices[1]))
                    seen_bonds.add(sorted_indices)
                    
    warning = ""
    if len(valid_bond_lengths) < expected_bonds:
        warning = "Fewer valid bonds than expected"
    elif len(valid_bond_lengths) > expected_bonds:
        warning = "More valid bonds than expected"

    if not valid_bond_lengths:
        warning += "; No valid bond lengths under 1.85 found. Check the structure."
    
    homa_values = []
    N = len(valid_bond_lengths)  # Number of valid bonds
    for length, atom1, atom2 in valid_bond_lengths:
        atom1_symbol = element_symbols[atom_types[atom1 - 1]]
        atom2_symbol = element_symbols[atom_types[atom2 - 1]]
        bond_type = tuple(sorted((atom1_symbol, atom2_symbol)))
        # print(f"Bond type: {bond_type}, Length: {length:.3f}")
        if bond_type in bond_parameters:
            R_opt, alpha = bond_parameters[bond_type]
            squared_deviation = (R_opt - length) ** 2
            homa_value = (alpha / N) * squared_deviation
            homa_values.append(homa_value)
        else:
            warning += f"; Unknown bond type {bond_type}."

    if homa_values:
        sum_of_homa_values = sum(homa_values)
        homa = 1 - sum_of_homa_values
    else:
        homa = None

    return homa, ring_atoms, expected_bonds, len(valid_bond_lengths), valid_bond_lengths, warning

def main():
    df = pd.read_csv('ring_atoms_final.csv')
    results = []

    for index, row in df.iterrows():
        file_name = row['File']
        if 'log' in file_name or not os.path.exists(file_name):
            continue

        ring_atoms = eval(row['ring_atoms'])
        
        try:
            homa, ring_atoms_str, expected_bonds, used_bonds, valid_bonds, warning = calculate_homa(file_name, ring_atoms)
            result = {
                'File': file_name,
                'HOMA': homa,
                'ring_atoms': ring_atoms_str,
                'Number of Expected Bonds': expected_bonds,
                'Number of Used Bonds': used_bonds,
                'Used Bond Lengths': [f"{bond[1]}-{bond[2]}: {bond[0]:.3f}" for bond in valid_bonds],
                'Warning': warning
            }
            results.append(result)

        except Exception as e:
            print(f"Failed to process {file_name}: {e}")

    results_df = pd.DataFrame(results)
    results_df.to_csv('homa.csv', index=False)
    print('HOMA calculations completed and saved to homa.csv.')

if __name__ == "__main__":
    main()
