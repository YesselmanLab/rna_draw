#!/bin/bash
#SBATCH -o /dev/null # Dump std out to null device
#SBATCH -e /dev/null # Dump std err to null device
#SBATCH -t 24:00:00
#SBATCH -n 1
#SBATCH -N 1
module load anaconda
conda activate py3
cd {JOB_DIR}

for dbn_file in *.dbn; do
    python $WORK/rna_draw/rna_draw/run_rna_draw.py "$dbn_file"
done