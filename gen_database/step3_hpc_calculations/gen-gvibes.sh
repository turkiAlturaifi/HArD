#!/usr/bin/env bash

# Loop from 1 to 240 with steps of 5
for i in $(seq 1 5 240); do
    # Calculate the end of the current range
    j=$((i + 4))

    # Create the file name based on the current range
    filename="gvibes_${i}-${j}.cmd"

    # Generate the .cmd file
    cat << EOF > "$filename"
#!/usr/bin/env bash
#SBATCH --job-name=g16
#SBATCH --output=goodvibes_${i}-${j}.out
#SBATCH --ntasks-per-node=4
#SBATCH --mem=8GB
#SBATCH --time=23:00:00
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --nodes=1

source {path_to_miniconda3}/miniconda3/etc/profile.d/conda.sh
conda activate aqme

# Loop through the parent folders ${i} through ${j}
for parent_folder in {$i..$j}; do
    if [ -d "\$parent_folder" ]; then
        pushd "\$parent_folder"
    else
        echo "Directory \$parent_folder does not exist."
        continue
    fi

    for folder in */; do
         ## change it check if the folder has failed_run1
        pushd "\$folder"
        cd opt_outputs
        python -m goodvibes *.log --csv && mv Goodvibes_output.csv ../gas.csv
        cp ../sp_sas/*.log .
        python -m goodvibes *.log --csv -c 1 -t 298.15 --spc SP && mv Goodvibes_output.csv ../m062x_sp_ssas.csv
        mv *SP* ../sp_ssas/
        cd ..
        popd
    done
    # Go back to the base directory 'opt'
    popd
done
EOF
done
