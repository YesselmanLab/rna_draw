#!/bin/bash
#SBATCH -o /dev/null # Dump std out to null device
#SBATCH -e /dev/null # Dump std err to null device
#SBATCH -t 24:00:00
#SBATCH -n 1
#SBATCH -N 1
module load anaconda
conda activate py3
cd {JOB_DIR}

failure_count=0
success_count=0

for dbn_file in *.dbn; do
    python $WORK/rna_draw/rna_draw/run_rna_draw.py "$dbn_file"

    if [ $? -eq 0 ]; then
        ((success_count++))
    else
        ((failure_count++))
    fi

    total=$((success_count + failure_count))

    if [ $total -ne 0 ]; then
        percent=$((100 * success_count / total))
    else
        percent=0
    fi

    echo "Success percentage: $percent%"
    echo "Total: $total"
    echo "Success: $success_count"
done