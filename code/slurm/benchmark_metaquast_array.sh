#!/bin/bash

#SBATCH --job-name=metaQUAST_benchmark     # name that will show up in the queue
#SBATCH --array=1-3
#SBATCH --output=slurm-benchmark-%j.out             # filename of the output; the %j is equal to jobID
#SBATCH --error=slurm-benchmark-%j.err              #
#SBATCH --partition=cpuqueue              #
#SBATCH --ntasks=1                        # number of tasks (analyses) to run
#SBATCH --cpus-per-task=12              # the number of threads allocated to each task
#SBATCH --mem-per-cpu=16G                  # memory per cpu-core
#SBATCH --time=08:15:00                   # time for analysis (day-hour:min:sec)
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

ref_genomes=${submitdir}/bsc_thesis/data/references/reference_genomes
gsa=${submitdir}/bsc_thesis/gold_standards/gsa/assembly
k21=${submitdir}/bsc_thesis/processed/assembly/k21
meta_sens=${submitdir}/bsc_thesis/processed/assembly/meta_sens


# 2. Execute [MODIFY COMPLETELY TO YOUR NEEDS]
# ============================================

# Decides which assemblies to evaluate
if [ $jobid == 1 ]; then
  assemblies=$gsa
elif [ $jobid == 2 ]; then
  assemblies=$k21
else
  assemblies=$meta_sens
fi

# Different host-contamination level
suffixes=("gsa" "k21" "meta_sens")
# Define variable based on array id
suffix=${suffixes[$jobid-1]}


metaquast -o metaquast_${suffix} \
          -r ${ref_genomes}\
          --threads 12 \
          --no-icarus \
          --unique-mapping \
          ${assemblies}/*
