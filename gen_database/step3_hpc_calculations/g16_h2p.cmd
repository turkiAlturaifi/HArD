#!/usr/bin/env bash
#SBATCH --job-name=g16
#SBATCH --output=g16.out
#SBATCH --ntasks-per-node=16
#SBATCH --mem=32GB
#SBATCH --time=95:00:00
#SBATCH --cluster=smp
#SBATCH --partition=smp
#SBATCH --nodes=1

module purge
module load gaussian/16-C.01

export GAUSS_SCRDIR=$SLURM_SCRATCH
ulimit -s unlimited
export LC_COLLATE=C

for filename in $SLURM_SUBMIT_DIR/*.com; do
    [ -e "$filename" ] || continue
    basename=$(basename "$filename" .com)
    echo "running $filename"
    cp "$filename" $SLURM_SCRATCH
    cd $SLURM_SCRATCH
    /usr/bin/time g16 < "$basename".com > "$SLURM_SUBMIT_DIR"/"$basename".log
    cd $SLURM_SUBMIT_DIR
done
