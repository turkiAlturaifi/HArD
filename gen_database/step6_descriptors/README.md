
## Overview

This repository contains scripts designed to generate various chemical descriptors from Gaussian log files and other molecular data.

### `global_features.py`
This script processes Gaussian log files to extract electronic properties such as:
- HOMO and LUMO energies
- HOMO-LUMO gap
- Chemical potential
- Global electrophilicity and nucleophilicity
- Total and directional dipole moments (maximum, median, minimum)
- Quadrupole moments (maximum, median, minimum, amplitude)

All dipole and quadrupole moments are taken as absolute values. The extracted properties are saved to `electronic.csv`.

### `homa.py`
Calculates the Harmonic Oscillator Model of Aromaticity (HOMA) descriptor as described in the manuscript. It requires `ring_atoms_final.csv` and outputs the results to `homa.csv`.

### `morfeus_descriptors.py`
Generates molecular descriptors including:
- Dispersion of the molecule and the ipso carbon
- Solvent accessible surface area (SASA) of the molecule and the ipso carbon
- Volume inside the SASA
- Fraction of buried volume and distal volume
- Sterimol parameters: L_value, min_width_sterimol_b1, max_width_sterimol_b5
- Bond lengths: C–H for unsubstituted heteroaryls and C(ipso)–C(carbonyl) for carboxylic acids and carboxylates

Three CSV files will be generated: `dispersion-sasa.csv`, `vbur.csv`, and `sterimol_params.csv`.

### `homo-lumo-coeff.py`
Calculates the HOMO and LUMO coefficients at the ipso carbon using `ring_atoms_final.csv` to locate the relevant atom number in log files. If no coefficient is found, the value is left empty.

### `charges_pt1.py`
Utilizes `ring_atoms_final.csv` to identify ring atoms and the ipso carbon, generating descriptors:
- `npa_ipso`: Natural Population Analysis charge on the ipso carbon
- `npa_ring_average`: Average NPA charge across the ring
- `mulliken_ipso`: Mulliken charge on the ipso carbon
- `mulliken_ring_average`: Average Mulliken charge across the ring
- `npa_hetaryl`: NPA charge on heteroaryl rings
- `mulliken_hetaryl`: Mulliken charge on heteroaryl rings
We didn't include the mulliken in the final set because many values didn't make sense but were exctarted correctly 

### `charges_pt2.py`
Utilizes `ring_atoms_final.csv` to identify ring atoms and the ipso carbon, generating descriptors using Hirshfeld and Charge Model 5 analyses:
- `hirshfeld_ipso`: Hirshfeld charge on the ipso carbon
- `hirshfeld_ring_average`: Average Hirshfeld charge across the ring
- `cm5_ipso`: Charge Model 5 on the ipso carbon
- `cm5_ring_average`: Average Charge Model 5 across the ring
- `hirshfeld_hetaryl`: Hirshfeld charge on heteroaryl rings
- `cm5_hetaryl`: Charge Model 5 on heteroaryl rings

## Requirements
- Python 3.x
- Morfeus
- pandas
- RDKit
- numpy
- cclib

For Morfeus, follow the installation instructions provided at [Morfeus's official documentation](https://github.com/digital-chemistry-laboratory/morfeus).

## Usage
To run any script, make sure to have single point log files and xyz files in the same dirctory then:
```
python script.py
```
