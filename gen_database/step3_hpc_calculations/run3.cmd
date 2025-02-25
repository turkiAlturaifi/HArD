#!/usr/bin/env bash
#SBATCH --job-name=g16
#SBATCH --output=g16_run3.out
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8GB
#SBATCH --time=23:00:00
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --nodes=1

# example for H2P supercomputer
source {path_to_miniconda3}/miniconda3/etc/profile.d/conda.sh
conda activate aqme

for parent_folder in {1..239}; do
    if [ -d "$parent_folder" ]; then
        pushd "$parent_folder"
    else
        echo "No errors in $parent_folder"
        echo "moving on to next one...."
        continue
    fi

    for folder in */; do
        if [ -d "$folder/run2_inputs" ]; then
            pushd "$folder"
        else
            echo "Directory $folder does not exist."
            continue
        fi
        cd run2_inputs

        echo "path: $folder"
        echo "number of log files: $(ls *.log 2>/dev/null | wc -l)"
        python -m aqme --qcorr --file "*.log"
        rm -f *.csv *.dat

        # Handle operations within 'failed' directory
        if [ -d "failed/run_1/extra_imag_freq" ]; then
            mv "failed/run_1/extra_imag_freq" failed_run1
        else
            echo "Directory failed/run_1/extra_imag_freq does not exist."
        fi

        if [ -d "failed/run_1/fixed_QM_inputs" ]; then
            mv "failed/run_1/fixed_QM_inputs" run3_inputs
        else
            echo "Directory failed/run_1/fixed_QM_inputs does not exist."
        fi

        if [ -d "failed/run_1/error" ]; then
            mv "failed/run_1/error" failed_run1/
        else
            echo "Directory failed/run_1/error does not exist."
        fi

        # Check if the 'failed' directory has any remaining log files
        if find failed -type f -name '*.log' | read; then
            echo "Log files found in the failed directory."
        else
            [ -d "failed" ] && rm -r failed
        fi

        # Handle other directory operations
        [ -d "inputs" ] && mv inputs opt_inputs || echo "Directory inputs does not exist."
        [ -d "success" ] && mv success opt_outputs || echo "Directory success does not exist."
        [ -d "opt_outputs/json_files" ] && rm -r "opt_outputs/json_files" || echo "Directory opt_outputs/json_files does not exist."

        # Copy g16_h2p.cmd if exists
        if [ -f "../g16_h2p.cmd" ]; then
            cp ../g16_h2p.cmd run3_inputs/
            if [ -d "run3_inputs" ]; then
                pushd "run3_inputs"
                sed -i 's/ime=95/ime=23/g' *.cmd
                sbatch g16_h2p.cmd 
                popd
            else
                echo "Directory run3_inputs does not exist."
            fi
        else
            echo "File g16_h2p.cmd does not exist in the parent folder."
        fi
        popd
    done
    popd
done
