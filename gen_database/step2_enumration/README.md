## Overview

This script processes heteroaryl molecules by adding CO2H and CO2- groups and then performing functional group substitutions using RDKit's `ReactionFromSmarts`. The final output includes processed derivatives with clear labels based on the substitution patterns.

The script `gen-smi.py` processes heteroaryl SMILES from batch1.csv by adding CO2H and CO2- groups to all possible ring atoms in the parent heteroaryl and performing functional group substitutions using RDKit's `ReactionFromSmarts`, while removing duplicates. 

Then the script `gen-opt.py` generates Gaussian input files with the suffix `_a` for the parent (H), `_b` for carboxylic acid, and `_c` for carboxylate derivatives.

## Functional Groups

The script adds the following functional groups to the heteroaryls, including the parent with H, and derivatives with CO2H and CO2-:

0.  no substitution (H) 
1.  NMe<sub>2</sub>  
2.  NH<sub>2</sub>  
3.  OH  
4.  OMe  
5.  Me  
6.  SiMe<sub>3</sub>  
7.  F  
8.  Cl  
9.  Br  
10. COCH<sub>3</sub> (acetyl)
11. CN  
12. NO<sub>2</sub>

## Unique Identifier System

Each molecule is assigned a unique identifier using the format: `a_b_c_d`

- `a`: Parent heteroaryl molecule
- `b`: Position of the substituent on the heteroaryl ring
- `c`: Identity of the substituent added (as described above)
- `d`: Regioisomer of the unsubstituted heteroaryl

## Requirements

- Python 3.x, `pandas`, `rdkit`, `natsort`

## Usage

Run `gen-smi.py` in the same dirctory as `batch1.csv` and `gen-opt.py`. it will generate `final.csv`, which has the IDs, SMILES, and parent names for ubstituted heteroaryl. Additionally two folder are generated with gaussian inputs and xyz files

   ```bash
   python gen-smi.py