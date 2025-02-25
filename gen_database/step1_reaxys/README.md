## Overview

`reaxys_filter.py` cleans up Reaxys search results by removing:

- Charged molecules
- Radicals
- Isotopes
- Molecules with atoms other than C, H, N, O, S
- Invalid valence SMILES: `N1N=NC=[NH]1` and `N1C=[NH]C2=NC=CC=C12`

## Requirements

- Python 3.x , `pandas`, `rdkit`

## Usage

Run the script in the same directory as `reaxys_search.xlsx`. It outputs `batch1.csv` and `reaxys_processed.csv`. Use `batch1.csv` as input for the next script.

   ```bash
   python reaxys_filter.py
   