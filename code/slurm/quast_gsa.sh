#!/bin/bash

#SBATCH --job-name=gsaQUAST_benchmark     # name that will show up in the queue
#SBATCH --array=1-11%5
#SBATCH --output=slurm-gsa-%j.out             # filename of the output; the %j is equal to jobID
#SBATCH --error=slurm-gsa-%j.err              #
#SBATCH --partition=cpuqueue              #
#SBATCH --ntasks=1                        # number of tasks (analyses) to run
#SBATCH --cpus-per-task=6              # the number of threads allocated to each task
#SBATCH --mem-per-cpu=8G                  # memory per cpu-core
#SBATCH --time=04:15:00                   # time for analysis (day-hour:min:sec)
#SBATCH --mail-type=ALL                   # send all type of email
#SBATCH --mail-user=andre.bourbonnais@sund.ku.dk


# I. Define directory names [DO NOT CHANGE]
# =========================================

# get the directories
submitdir=${SLURM_SUBMIT_DIR}
workdir=${TMPDIR}
jobid=${SLURM_ARRAY_TASK_ID}

# Information
echo "$(date)    Submitted from ${submitdir}"
echo "$(date)    Accessed ${workdir}"
echo "$(date)    ArrayID: ${jobid}"


# 1. Lock and load module and data
# ============================================
module load quast

gsa=${submitdir}/bsc_thesis/gold_standards/gsa/assembly
k21=${submitdir}/bsc_thesis/processed/assembly/k21
meta_sens=${submitdir}/bsc_thesis/processed/assembly/meta_sens


# 2. Execute [MODIFY COMPLETELY TO YOUR NEEDS]
# ============================================

# Different host-contamination level
hc_level=("00" "01" "02" "03" "04" "05" "06" "07" "08" "09" "095")

# Define prefix based on array id
hc_prefix=${hc_level[$jobid-1]}


quast -o quast_${hc_prefix} \
        -r ${gsa}/${hc_prefix}_* \
        --threads 6 \
        --no-icarus \
        ${k21}/${hc_prefix}_* \
        ${meta_sens}/${hc_prefix}_*

echo "$(date)    HC: ${hc_prefix}"
