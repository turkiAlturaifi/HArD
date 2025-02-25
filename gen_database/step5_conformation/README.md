
## Overview

For conformers, the dihedral angle is adjusted to consider both coplanar and perpendicular geometries for carboxylic acid and carboxylate-containing heteroaryls. For conformers, the dihedral angle is adjusted to consider both coplanar and perpendicular geometries for carboxylic acid and carboxylate-containing heteroaryls. Additionally, for the parent heteroaryl, duplicate log files with different names are generated to easily calculate the descriptors for different regioisomers.

## Requirements

- `rdkit 2022.09.3`

## Usage

1. Ensure the file `ring_atoms_final.csv` is in the same directory. This script assesses whether the ring and carboxylic acid are coplanar. If they are, it generates the other rotamer. Otherwise, it generates two rotamers by rotating about the C(ipso)–C(carbonyl) axis by 180 degrees.

   ```bash
   python calc_dihedral_b.py
   ```

2. This script adjusts coplanar conformers to be 90 degrees and generates a new `.xyz` file.

   ```bash
   python calc_dihedral_c.py
   ```

3. generating Gaussian input files and placing them into two directories: `confs_b` and `confs_c`.

   ```bash
   python gen-com.py
   ```
follow the instructions in step3_hpc_calculations

4. Generate regioisomers for the parent heteroaryls for descriptor calculations
   ```bash
   python regioisomers.py
   ```
This will create regioisomers folder with log files that have the correct names. Add these files to the master log files folder. The regioisomers are just the duplicate sturctures that were removed from step 1 and now they will be used calculate descriptors in the next step with the other files