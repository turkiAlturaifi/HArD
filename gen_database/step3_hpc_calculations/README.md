## Overview

This workflow demonstrates how to perform parallel high-performance computing (HPC) calculations using SLURM job scheduling. Gaussian calculations are submitted in batches of 200 files per job, checked for errors, and resubmitted up to two times if necessary. After successful convergence, different single-point calculations are performed for further analysis.

## Requirements

- `aqme`
– `open babel`

## Usage

1. Ensure the `com` directory from the previous step, along with `separate.sh` and `g16_h2p.cmd`, are in the current working directory.
   
   Run the `hpc-folders.py` script to organize files and prepare batch job submission:
   ```bash
   python hpc-folders.py
   ```

2. Submit the jobs to the HPC cluster:
   ```bash
   chmod u+x submit_jobs.sh && ./submit_jobs.sh
   ```

3. Once the jobs are completed, run the `gen-run2-cmd.sh` script to generate SLURM jobs for resubmitting failed calculations:
   ```bash
   ./gen-run2-cmd.sh
   # Submit the resubmission jobs
   ```

4. After the second round of jobs is completed, run `run3.cmd` to finalize all calculations:
   ```bash
   ./run3.cmd
   ```

5. Mark any remaining failed calculations using the `mark_failed_opt.sh` script:
   ```bash
   ./mark_failed_opt.sh
   ```
6. to validate the connectivity, the number of rings before and after optimization is compared with the following script:
   ```bash
   python ring_count_check.py # first check
   python connectivity.py # second check
   ```
which will generate a csv file "ring_counts_with_warnings.csv" that labels any geometries that determine to be 'bad'

## Single-Point Calculations

The following SLURM job scripts are provided for performing single-point calculations:
- `g16_sp.cmd`: For single-point energy corrections and Natural Population Analysis (NPA) charges.
- `g16_homo-lumo-cof.cmd`: For calculating HOMO-LUMO coefficients.
- `g16_hirsh.cmd`: For Hirshfeld charge analysis.
- `gen-gvibes.sh`: to generate slurm files to compute the energies
- `g16_sp_ssas.cmd`: to run single point calculations with scaled solvent accessible surface area for sigma het


