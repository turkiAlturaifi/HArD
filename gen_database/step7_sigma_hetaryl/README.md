## Overview

This repository contains workflow and scripts to generate Hammett-type substituent constants. The first step is to extract the energies using goodvibes then compute the pKa and the sigma_het. The sigma_het is computed by using the energies of the carboxylic acid and carboxylate anion containing molecules, determining their pKa and using the pKa of benzoic acid as reference to yield the descriptor.

## Requirements
- aqme (or goodvibes)

## Usage

1. After running `gen-gvibes.sh` from "step3_hpc_calculations", you should have a number of `goodvibes_i-j.out` files where `i` and `j` are the range of the folders. Combine these files into a master file named "gvibes_all.out". Then run:
   ```bash
   python gvibes_energies.py
   ```
   This will create two CSV files: `energies_master.csv` with the energies in gas and single point, and `energies.csv` which has a summary of the energies to be used for the next step.
   
2. Generate sigma_hetaryl:
   ```bash
   python sigma-het.py
   ```
   This will generate `sigma_het.csv` containing the pKa and Hammett-type substituent constants.