#!/bin/bash
#SBATCH -o /dev/null # Dump std out to null device
#SBATCH -e /dev/null # Dump std err to null device
#SBATCH -t 24:00:00
#SBATCH -n 1
#SBATCH -N 1
module load anaconda
conda activate py3
cd {JOB_DIR}

# Loop through all DBN files in the current directory
for dbn_file in *.dbn; do
    # Call the Python script for each DBN file
    python $WORK/rna_draw/rna_draw/run_rna_draw.py "$dbn_file"
done


#rna_draw 
#rna-struct-design helix-rand --csv-file input.csv --param-file {WORK_DIR}/resources/helix-rand.yml --output output.csv --num-seqs 3