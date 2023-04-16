#!/bin/bash

#SBATCH --job-name=k21_megahit     # name that will show up in the queue
#SBATCH --array=1-11%4
#SBATCH --output=slurm-%j.out             # filename of the output; the %j is equal to jobID
#SBATCH --error=slurm-%j.err              #
#SBATCH --partition=cpuqueue              #
#SBATCH --ntasks=1                        # number of tasks (analyses) to run
#SBATCH --cpus-per-task=6               # the number of threads allocated to each task
#SBATCH --mem-per-cpu=8G                  # memory per cpu-core
#SBATCH --time=01:30:00                   # time for analysis (day-hour:min:sec)
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
module load megahit
# Get the trimmed reads
reads=${submitdir}/bsc_thesis/data/subsample/reads/01_trimmed/



# 2. Execute [MODIFY COMPLETELY TO YOUR NEEDS]
# ============================================

# Different host-contamination level
hc_level=("00" "01" "02" "03" "04" "05" "06" "07" "08" "09" "095")

# Define prefix based on array id
hc_prefix=${hc_level[$jobid-1]}


megahit -t 6 --k-min 21 --k-max 91 \
    -1 ${reads}/${hc_prefix}_trim_R1.fq.gz \
    -2 ${reads}/${hc_prefix}_trim_R2.fq.gz \
    -o k21/${hc_prefix} \
    --out-prefix "${hc_prefix}_k21"
