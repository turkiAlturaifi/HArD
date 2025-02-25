#!/bin/bash
for folder in {1..100}; do
    if [ -d "$folder" ] && [ -f "$folder/seperate.sh" ]; then
        echo "Processing folder $folder"
        (cd "$folder" && chmod u+x seperate.sh && ./seperate.sh)
    else
        echo "Skipping folder $folder: seperate.sh not found"
    fi
done

