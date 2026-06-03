#!/bin/bash
# =============================================================================
# AllenLab CRISPR MAGeCK Pipeline — SLURM submission script
# TSCC | Partition: condo | Allocation: sio141
#
# Submit:  sbatch mageck_run.sh
# Monitor: squeue -u a8gupta
# Cancel:  scancel <job_id>
# Logs:    tail -f logs/slurm-<job_id>.out
# =============================================================================

#SBATCH -p condo
#SBATCH -q condo
#SBATCH -A sio141
#SBATCH -N 1 -n 1 -c 32
#SBATCH --mem=128gb
#SBATCH -J CRISPR_Pipeline
#SBATCH -t 7-00:00:00
#SBATCH -o logs/slurm-%j_%x.out
#SBATCH -e logs/slurm-%j_%x.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=a8gupta@ucsd.edu
#SBATCH --requeue

# =============================================================================
# 1. Environment Setup
# =============================================================================
echo "============================================================"
echo "  Job ID   : $SLURM_JOB_ID"
echo "  Node     : $(hostname)"
echo "  Start    : $(date)"
echo "  Dir      : $SLURM_SUBMIT_DIR"
echo "============================================================"

module purge
module load cpu slurm gcc

# Initialise conda for this shell
source ~/.bashrc
conda activate mageckenv

# Lab-specific binaries (MAGeCK etc.)
export PATH=$PATH:/tscc/projects/ps-allenlab/projdata/common/bin

# =============================================================================
# 2. Directory setup
# =============================================================================
cd $SLURM_SUBMIT_DIR
mkdir -p logs

# =============================================================================
# 3. Run pipeline
# =============================================================================
echo ""
echo "--- Pipeline started at $(date) ---"
echo ""

python mageck_pipeline.py

EXIT_CODE=$?

echo ""
echo "============================================================"
if [ $EXIT_CODE -eq 0 ]; then
    echo "  ✅ PIPELINE COMPLETE  — $(date)"
else
    echo "  ❌ PIPELINE FAILED (exit $EXIT_CODE) — $(date)"
    echo "  Check logs/slurm-${SLURM_JOB_ID}_CRISPR_Pipeline.err for details."
fi
echo "============================================================"
exit $EXIT_CODE


# =============================================================================
# Quick reference
# =============================================================================
# Submit job:
#   sbatch mageck_run.sh
#
# Check job status:
#   squeue -u a8gupta
#
# Watch live log output:
#   tail -f logs/slurm-<job_id>_CRISPR_Pipeline.out
#
# Cancel job:
#   scancel <job_id>
#
# Check memory usage (SSH to compute node first):
#   ssh <node_name>   # e.g. tscc-9-16
#   top -u a8gupta    # press M to sort by memory, q to quit
#
# Run interactively (for testing / debugging):
#   srun -p condo -q condo -A sio141 -N 1 -n 1 -c 4 --mem=16gb -t 2:00:00 --pty bash
#   conda activate mageckenv
#   python mageck_pipeline.py
