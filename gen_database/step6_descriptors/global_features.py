import os
import csv
import cclib
import numpy as np

output_csv = 'electronic.csv'
log_files = [f for f in os.listdir('.') if f.endswith('.log')]
with open(output_csv, mode='w', newline='') as csv_file:
    csv_writer = csv.writer(csv_file)
    # Write the header
    csv_writer.writerow([
        'File', 'homo', 'lumo', 'homo-lumo_gap', 'chemical_potential',
            'global_electrophilicity', 'global_nucleophilicity', 'total_dipole_moment',
            'dipole_moment_max', 'dipole_moment_middle', 'dipole_moment_min',
            'quadrupole_moment_max', 'quadrupole_moment_middle', 'quadrupole_moment_min',
            'quadrupole_moment_amplitude'
        ])

    for log_file in log_files:
        try:
            data = cclib.io.ccread(log_file)
            homo_index = data.homos[0]
            lumo_index = homo_index + 1
            homo_energy = data.moenergies[0][homo_index]
            lumo_energy = data.moenergies[0][lumo_index]

            gap = lumo_energy - homo_energy
            chemical_potential = (lumo_energy + homo_energy) / 2
            electrophilicity = chemical_potential**2 / (2 * gap)
            nucleophilicity = 1 / electrophilicity

            dipole_values = [abs(x) for x in data.moments[1]]
            sorted_dipole_values = sorted(dipole_values, reverse=True)
            dipole_x, dipole_y, dipole_z = sorted_dipole_values[0], sorted_dipole_values[1], sorted_dipole_values[2]
            total_dipole_moment = np.sqrt(dipole_x**2 + dipole_y**2 + dipole_z**2)

            # quadrupole moment
            if hasattr(data, 'moments') and len(data.moments) > 2:
                quadrupole_values = [abs(x) for x in data.moments[2]]
                sorted_quadrupole_values = sorted(quadrupole_values, reverse=True)
                quadrupole_xx, quadrupole_yy, quadrupole_zz = sorted_quadrupole_values[0], sorted_quadrupole_values[1], sorted_quadrupole_values[2]
                quadrupole_amplitude = np.sqrt(quadrupole_xx**2 + quadrupole_yy**2 + quadrupole_zz**2)
            else:
                quadrupole_xx = quadrupole_yy = quadrupole_zz = quadrupole_amplitude = np.nan

            csv_writer.writerow([
                log_file, homo_energy, lumo_energy, gap, chemical_potential,
                electrophilicity, nucleophilicity, total_dipole_moment,
                dipole_x, dipole_y, dipole_z,
                quadrupole_xx, quadrupole_yy, quadrupole_zz, quadrupole_amplitude
            ])

        except Exception as e:
            print(f"Error processing {log_file}: {e}")

print(f"Properties have been written to {output_csv}")