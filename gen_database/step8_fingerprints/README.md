## Overview

This repository contains scripts to generate fingerprint descriptors and substitution SMILES for heteroarenes, ensuring entries for structures differing only by H/CO2H/CO2- are consolidated. The `fingerprint.py` script outputs the following descriptors:

- Number_of_rings
- Number_of_heteroatoms
- Number_of_hydrogen_bond_donors
- Number_of_hydrogen_bond_acceptors
- Topological_polar_surface_area_tpsa
- Molecular_weight
- LogP
- Rotatable_bonds_count
- Number_of_heavy_atoms
- Number_of_carbon_atoms
- Number_of_nitrogen_atoms
- Number_of_oxygen_atoms
- Number_of_sulfur_atoms
- Fraction_of_aromatic_bonds
- Fraction_of_sp3_hybridized_carbons
- Fraction_of_sp2_hybridized_carbons

## Requirements

- `rdkit 2022.09.3`

## Usage

1. Run `fingerprint.py` in the same folder as `final.csv` from "step1_reaxys", as well as `wild_smi.py`:

```bash
python fingerprint.py  # Generates fingerprint.csv
python wild_smi.py     # Generates wild_smi.csv
```

2. Consolidate the results for the descriptors using the script below. 
```bash
python combine.py vbur.csv # Generates vbur_sorted.csv
python combine.py charges-npa.csv # Generates charges-npa_sorted.csv

Where file.csv could be "vbur.csv", "charges-npa.csv", or any other computed descriptor file.

### Final note
The script processes descriptors for all files. In the final version, descriptors with suffixes _b and _c are excluded, except for charges and bond lengths, which are retained due to their importance in models like pKa prediction.

