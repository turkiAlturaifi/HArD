### Structure of `gen_database`

1. **step1_reaxys**
   - Scripts and results for filtering and processing raw data from Reaxys exports.
   - Files: Batch CSVs, Python scripts for processing, and documentation in README.

2. **step2_enumration**
   - Scripts for generating SMILES strings and other preliminary data.
   - Files: Enumeration scripts and final CSV outputs.

3. **step3_hpc_calculations**
   - High-performance computing scripts for executing computational chemistry tasks like geometry optimization and property calculation. You have to change this a little bit based on which HPC you're using like the cluster names etc.
   - Files: job submission scripts, and utility scripts for managing HPC tasks. 

4. **step4_atom_index**
   - to find key atom numbers (ring atoms, ipso carbon etc)
   - Files: Scripts for determining ring atoms and substitution sites.

5. **step5_conformation**
   - Scripts for calculating conformations and dihedral angles.
   - Files: python scripts and csv files

6. **step6_descriptors**
   - Scripts for generating various molecular descriptors.
   - Files: Python scripts for charges, HOMA, HOMO-LUMO coefficients, and more.

7. **step7_sigma_hetaryl**
   - Tools for calculating sigma constants specific to heteroarenes.
   - Files: Scripts for deriving and managing vibrational energies and sigma constants.

8. **step8_fingerprints**
   - Scripts for creating fingerprint descriptors and handling SMILES strings for database entries.
   - Files: Combination and fingerprint generation scripts.
