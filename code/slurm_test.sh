#!/bin/bash

#SBATCH --job-name=06_quast_test    # name that will show up in the queue
#SBATCH --output=slurm-%j.out       # filename of the output; the %j is equal to jobID
#SBATCH --error=slurm-%j.err        #
#SBATCH --partition=cpuqueue        #
#SBATCH --ntasks=1                  # number of tasks (analyses) to run
#SBATCH --cpus-per-task=12           # the number of threads allocated to each task
#SBATCH --mem-per-cpu=8G           # memory per cpu-core
#SBATCH --time=02:30:00             # time for analysis (day-hour:min:sec)
#SBATCH --mail-type=ALL             # send all type of email
#SBATCH --mail-user=andre.bourbonnais@sund.ku.dk


# I. Define directory names [DO NOT CHANGE]
# =========================================

# get name of the temporary directory working directory, physically on the compute-node
workdir="${TMPDIR}"
echo "Accessed ${workdir} $(date)"

# get submit directory
# (every file/folder below this directory is copied to the compute node)
submitdir="${SLURM_SUBMIT_DIR}"

# 1. Transfer to node [DO NOT CHANGE]
# ===================================

# create/empty the temporary directory on the compute node
if [ ! -d "${workdir}" ]; then
  echo "workdir does not exist"
else
  echo "workdir exists, error?"
  ls "${workdir}"
fi

# change current directory to the location of the sbatch command
# ("submitdir" is somewhere in the home directory on the head node)
cd "${submitdir}"
# copy only the required files/folders in "submitdir" to "workdir"
# ("workdir" == temporary directory on the compute node)
cp -prf gsa megahit_results "${workdir}"
# change directory to the temporary directory on the compute-node
cd ${workdir}

# 3. Function to transfer back to the head node [DO NOT CHANGE]
# =============================================================

# define clean-up function
function clean_up {
  # - delete temporary files from the compute-node, before copying
  # rm -r ...
  # - change directory to the location of the sbatch command (on the head node)
  cd "${workdir}"
  # - copy only the required files/folders from the temporary directory on the compute-node
  cp -prf quast_results "${submitdir}"
  # - erase the temporary directory from the compute-node
  rm -rf gsa megahit_results quast_results
  # - exit the script
  exit
}

# call "clean_up" function when this script exits, it is run even if SLURM cancels the job
trap 'clean_up' EXIT

# 2. Execute [MODIFY COMPLETELY TO YOUR NEEDS]
# ============================================

# Load latest version of quast
module load quast

quast -o quast_results \
          -r 00_gsa.fasta \
          --threads 12 \
          --unique-mapping\
          --no-icarus \
          --min-alignment 200\
          megahit_results/platanthera_mock_assembly/00/00.contigs.fa