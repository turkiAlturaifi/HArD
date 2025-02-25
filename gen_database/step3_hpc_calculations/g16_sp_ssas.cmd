#!/usr/bin/env bash
#SBATCH --job-name=g16
#SBATCH --output=g16_sp_ssas.out
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8GB
#SBATCH --time=23:00:00
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --nodes=1

source {path_to_miniconda3}/miniconda3/etc/profile.d/conda.sh
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
        python -m aqme --qprep --files "*.log" --qm_input "m062x/6-31+G(d) scrf=(smd,solvent=water,read)" --mem 16GB --nproc 8 --suffix SP_ssas --program gaussian --qm_end "Surface=SAS"$'\n'"Alpha=0.485"
        mv QCALC sp_ssas
        mv sp_ssas ../
        cd ../sp_ssas
        cp ../g16_h2p.cmd .
        sed -i 's/r-node=16/r-node=8/g' *.cmd
        sed -i 's/ime=95/ime=23/g' *.cmd
        sbatch g16_h2p.cmd
        cd ..
        popd
    done
    popd
done
