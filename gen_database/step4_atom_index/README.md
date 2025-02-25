## Overview

This section handles atom indices and generates two CSV files that are used to extract descriptors in later steps. The script identifies the following information:

- **Ring Atoms**: Atom indices in the ring.
- **Atom with FG1**: The ring atom if it has a substitution.
- **Atom with FG2**: Another ring atom with a substitution, if any.
- **Attached to Carboxyl Group**: Boolean indicating if the ring atom is attached to a carboxylic group.
- **metal_index**: The ipso carbon, which will be used for various descriptors such as HOMO/LUMO coefficients and Sterimol.
- **xz_atom**: The atom on the ring in the xz plane used to compute Sterimol and buried volume. It is identified as the atom directly bonded to the ipso carbon in the ring.
- **z_axis_atom**: Used for vbur and Sterimol calculations. It is the hydrogen atom for the parent or the carbon for carboxylic acid and carboxylate anion, which is bonded to the ipso carbon.
- **excluded_atoms**: CO2 atoms to be excluded from vbur and Sterimol calculations, along with the ipso carbon.
- **CO2H_atoms**: The CO2H atoms, necessary for conformations.

The second CSV file is needed to compute descriptors for regioisomers of the parent molecule.

## Requirements

- `rdkit 2022.09.3`

## Usage

1. Convert all Gaussian outputs to XYZ format:
   
   ```bash
   obabel *.log -O .xyz -m
   ```

2. Make sure to have the file final.csv (from step 2) in the current working dirctory. Run the following scripts, which will generate the ipso carbon atom number in the parent that corresponds to the formylation site. This is required to calculate atom-specific descriptors:
   
   ```bash
   python substitution_site.py
   python ring_atoms.py
   ```
The script will generate "ring_atoms_final.csv" which will be used for most descriptor calculations and "ipso_a.csv" which will be used in next step