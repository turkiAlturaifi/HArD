#!/bin/bash

for parent_folder in {1..250}; do
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
        mv run2_inputs run3
        cd run3
        [ -d "opt_inputs" ] && mv opt_inputs ../opt_inputs_run2 || echo "Directory opt_inputs does not exist."
        [ -d "opt_outputs" ] && mv opt_outputs ../opt_outputs_run2 || echo "Directory opt_outputs does not exist."
        [ -d "failed_run1" ] && mv failed_run1 ../failed_run2 || echo "Directory failed_run2 does not exist."
        [ -d "run3_inputs" ] && mv run3_inputs/* . && rm -r run3_inputs || echo "Directory opt_inputs does not exist."
        popd
    done
    popd
done
