#!/bin/bash
#SBATCH --job-name=plantomicsgwas
#SBATCH --partition=compute
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=32
#SBATCH --mem=128G
#SBATCH --time=48:00:00


set -e

echo "======================================"
echo " PlantOmicsGWAS HPC Job"
echo " Scheduler: SLURM"
echo "======================================"

cd "/home/ahmed/Documents/Projects/PlantOmicsGWAS/examples"

source ~/.bashrc

conda activate plantomicsgwas-hpc



plantomicsgwas-compute run --config "/home/ahmed/Documents/Projects/PlantOmicsGWAS/examples/hpc_config_template.yaml"

echo "Finished."
