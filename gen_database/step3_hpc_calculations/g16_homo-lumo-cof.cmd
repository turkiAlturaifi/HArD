#!/usr/bin/env bash
#SBATCH --job-name=g16
#SBATCH --output=g16_sp_submit_homo-lumo.out
#SBATCH --ntasks-per-node=8
#SBATCH --mem=14GB
#SBATCH --time=23:00:00
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --nodes=1

source /ihome/pliu/tma53/miniconda3/etc/profile.d/conda.sh
conda activate aqme

for parent_folder in {1..250}; do
    if [ -d "$parent_folder" ]; then
        pushd "$parent_folder"
    else
        echo "No errors in $parent_folder"
        echo "moving on to next one...."
        continue
    fi
    
    for folder in */; do
        pushd "$folder"
        cd opt_outputs
        python -m aqme --qprep --files "*.log" --qm_input "m062x/6-31+G(d) scrf=(smd,solvent=water) pop=(orbitals=10,ThreshOrbitals=5)" --mem 14GB --nproc 8 --suffix SP_homo-lumo --program gaussian
        mv QCALC sp_homo-lumo
        mv sp_homo-lumo ../
        cd ../sp_homo-lumo
        cp ../g16_h2p.cmd .
        sed -i 's/r-node=16/r-node=8/g' *.cmd
        sed -i 's/ime=95/ime=23/g' *.cmd
        sbatch g16_h2p.cmd
        cd ..
        popd
    done
    popd
done
