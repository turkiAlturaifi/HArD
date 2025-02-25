import pandas as pd
from morfeus import Dispersion, read_xyz, SASA

# Read the specific CSV file for ring_atoms7_b
csv_file_path = 'combined_ring.csv'
csv_data = pd.read_csv(csv_file_path)
# csv_data = csv_data[~csv_data['File'].str.contains('SP')]

# Function to calculate dispersion and SASA properties from XYZ file
def calculate_dispersion_properties(xyz_file_path, row):
    try:
        elements, coordinates = read_xyz(xyz_file_path)
        disp = Dispersion(elements, coordinates)
        sasa = SASA(elements, coordinates)

        metal_index = int(row['metal_index']) if not pd.isna(row['metal_index']) else None
        properties = {
            'File': row['File'],
            'dispersion_ipso': disp.atom_p_int[metal_index] if metal_index is not None else None,
            'dispersion': disp.p_int,
            'disp_surface_area': disp.area,
            'disp_volume': disp.volume,
            'sasa': sasa.area,
            'volume_inside_sas': sasa.volume,
            'sasa_ipso': sasa.atom_areas[metal_index] if metal_index is not None else None
        }
        return properties
    except Exception as e:
        print(f"Error calculating dispersion for {xyz_file_path}: {e}")
        return None

# Main function to process files
def main():
    output_csv_file = 'dispersion-sasa.csv'
    all_results = []

    for _, row in csv_data.iterrows():
        log_file = row['File']
        xyz_file_path = log_file.replace('.log', '.xyz')

        properties = calculate_dispersion_properties(xyz_file_path, row)
        if properties:
            all_results.append(properties)

    # Convert all results to a DataFrame and save to CSV
    results_df = pd.DataFrame(all_results)
    results_df.to_csv(output_csv_file, index=False)
    print("created dispersion-sasa.csv")

if __name__ == '__main__':
    main()


import os
import csv
from morfeus import BuriedVolume, read_xyz

print('-----section 2---- generate vbur')

def process_file(file_name, metal_index, z_axis_atom, xz_plane_atom, excluded_atoms):
    try:
        elements, coordinates = read_xyz(file_name)
        bv = BuriedVolume(
            elements, 
            coordinates, 
            metal_index=metal_index, 
            excluded_atoms=excluded_atoms,
            z_axis_atoms=[z_axis_atom],
            xz_plane_atoms=[xz_plane_atom],
            radius=3.5
        )
        bv.compute_distal_volume(method="buried_volume")
        data_for_csv = {
            'vbur': {
                'File': file_name,
                'Program': 'morfeus',
                'Metal index': metal_index,
                'Z-axis atom': z_axis_atom,
                'XZ-plane atom': xz_plane_atom,
                'Buried volume': float(bv.buried_volume),
                'distal_volume': bv.distal_volume,
                'fraction_buried_volume': bv.fraction_buried_volume
            },
        }
        return data_for_csv
    except Exception as e:
        print(f"Error processing file {file_name}: {e}")
        return None

def append_dict_to_csv(file_name, data, filename, key):
    file_exists = os.path.isfile(filename)
    with open(filename, 'a', newline='') as csvfile:
        writer = csv.writer(csvfile)
        headers = {
            'vbur': ['File', 'Program', 'Metal index', 'Z-axis atom', 'XZ-plane atom', 'Buried volume', 'distal_volume', 'fraction_buried_volume'],
        }
        if not file_exists:
            writer.writerow(headers[key])
        row = [file_name] + [data.get(h, 0) for h in headers[key][1:]]
        writer.writerow(row)

def main_vbur():
    input_csv_file = 'combined_ring.csv'
    log_files_info = []

    with open(input_csv_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            log_file = row['File']
            log_files_info.append(row)
    
    for file_info in log_files_info:
        xyz_file = file_info['File'].replace('.log', '.xyz')
        try:
            metal_index = int(float(file_info['metal_index']))
            xz_atom = int(float(file_info['xz_atom']))
            z_axis_atom = int(float(file_info['z_axis_atom']))
            excluded_atoms = [int(float(atom)) for atom in file_info['excluded_atoms'].strip('[]').split(',') if atom.strip()]
        except Exception as e:
            print(f'Skipped {file_info["File"]} due to error: {e}')
            continue

        data_for_csv = process_file(xyz_file, metal_index, z_axis_atom, xz_atom, excluded_atoms)
        if data_for_csv:
            for key, value in data_for_csv.items():
                filename = f"{key}.csv"
                append_dict_to_csv(xyz_file, value, filename, key)

if __name__ == '__main__':
    main_vbur()


print('-----section 3---- generate sterimol')

import os
import csv
from morfeus import Sterimol, read_xyz

def process_sterimol(file_name, dummy_index, attached_index, excluded_atoms):
    try:
        elements, coordinates = read_xyz(file_name)
        
        # Set up Sterimol instance and perform analysis
        sterimol = Sterimol(
            elements, 
            coordinates, 
            dummy_index=dummy_index, 
            attached_index=attached_index, 
            excluded_atoms=excluded_atoms
        )

        # Extract Sterimol parameters
        sterimol_params = {
            'File': file_name,
            'L_value': sterimol.L_value,
            'min_width_sterimol_b1': sterimol.B_1_value,
            'max_width_sterimol_b5': sterimol.B_5_value,
            'bond_length': sterimol.bond_length
        }

        return sterimol_params
    except Exception as e:
        print(f"Error processing Sterimol for {file_name}: {e}")
        return None

def main_sterimol():
    # Read from ring_atoms7_updated.csv
    input_csv_file = 'combined_ring.csv'
    log_files_info = []

    with open(input_csv_file, 'r') as csvfile:
        reader = csv.DictReader(csvfile)
        for row in reader:
            log_files_info.append(row)

    sterimol_data = []

    for file_info in log_files_info:
        try:
            log_file = file_info['File']
            xyz_file = log_file.replace('.log', '.xyz')

            metal_index = int(float(file_info['metal_index']))
            xz_atom = int(float(file_info['xz_atom']))
            z_axis_atom = int(float(file_info['z_axis_atom']))
            excluded_atoms = [int(float(atom)) for atom in file_info['excluded_atoms'].strip('[]').split(',') if atom.strip()]

            sterimol_params = process_sterimol(xyz_file, z_axis_atom, metal_index, excluded_atoms)
            if sterimol_params:
                sterimol_data.append(sterimol_params)
        except Exception as e:
            print(f'Skipped {log_file} due to error: {e}')
            continue

    # Write Sterimol parameters to CSV
    output_csv_file = 'sterimol_params.csv'
    with open(output_csv_file, 'w', newline='') as csvfile:
        fieldnames = ['File', 'L_value', 'min_width_sterimol_b1', 'max_width_sterimol_b5', 'bond_length']
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

        writer.writeheader()
        for data in sterimol_data:
            writer.writerow(data)

if __name__ == '__main__':
    main_sterimol()
