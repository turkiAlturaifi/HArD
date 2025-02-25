#!/bin/bash

if [ ! -f "g16_h2p.cmd" ]; then
    echo "g16_h2p.cmd file not found in the current directory."
    exit 1
fi
readarray -t files < <(printf '%s\n' *.com | sort -V)

batchSize=200

numFiles=${#files[@]}
numBatches=$((numFiles / batchSize))
if ((numFiles % batchSize != 0)); then
    ((numBatches++))
fi

for ((batch=1; batch<=numBatches; batch++)); do
    startIndex=$(((batch - 1) * batchSize))
    endIndex=$((startIndex + batchSize - 1))
    if [ $endIndex -ge $numFiles ]; then
        endIndex=$((numFiles - 1))
    fi

    folderName="${startIndex}-${endIndex}"
    mkdir -p "$folderName"
    cp "g16_h2p.cmd" "$folderName"
    for ((i=startIndex; i<=endIndex; i++)); do
        mv "${files[$i]}" "$folderName"
    done
done
echo "Operation completed successfully."
for folder in */; do cd "$folder" && sbatch g16_h2p.cmd && cd ..; done