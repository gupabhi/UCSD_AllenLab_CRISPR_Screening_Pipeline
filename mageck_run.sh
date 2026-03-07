#!/bin/bash
#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A sio141
#SBATCH -N 1 -n 1 -c 32
#SBATCH --mem=128gb
#SBATCH -J CRISPR_Pipeline_Run
#SBATCH -t 7-00:00:00
#SBATCH -o logs/slurm-%j.out
#SBATCH -e logs/slurm-%j.err
#SBATCH --mail-type END,FAIL

# ---------------------------------------------------------
# 1. Environment Preparation
# ---------------------------------------------------------
echo "Initializing environment..."
module purge
module load cpu slurm gcc

# Initialize Conda for your shell
source ~/.bashrc

# Activate your specific environment
conda activate mageckenv

# Add lab-specific bin folder to path
export PATH=$PATH:/tscc/projects/ps-allenlab/projdata/common/bin

# ---------------------------------------------------------
# 2. Directory Management
# ---------------------------------------------------------
# Move to the directory from which the script was submitted
cd $SLURM_SUBMIT_DIR

# Ensure necessary directories exist
mkdir -p logs

# ---------------------------------------------------------
# 3. Execution Pipeline
# ---------------------------------------------------------
echo "Pipeline started at $(date)"
echo "Running on node: $(hostname)"

# ---------------------------------------------------------
# ---------------------------------------------------------
# ---------------------------------------------------------

echo "------------------------------------------------"
echo "Running CRISPR MAGeCK Pipeline Steps"
python mageck_pipeline.py
if [ $? -ne 0 ]; then 
    echo "CRITICAL: Running CRISPR MAGeCK Pipeline FAILED. Exiting."
    exit 1
fi

# ---------------------------------------------------------

# ---------------------------------------------------------
#END. Completion
echo "------------------------------------------------"
echo "PIPELINE SUCCESSFUL"
echo "Finished at $(date)"


## Command Notes:
# To submit the job which runs the .sh script: "sbtach run_crispr.sh"
# To check the status of the job: "squeue -u a8gupta"
# To cancel the job: "scancel <job_id>"
# To view memory usage: "ssh tscc-9-16" and then "top -u a8gupta" (press 'M' to sort by memory)
