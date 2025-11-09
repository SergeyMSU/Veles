#! /bin/bash
#SBATCH --job-name="Veles2"
#SBATCH --gpus=1
#SBATCH --nodes=1
#SBATCH --gpus-per-node=1
#SBATCH --ntasks=1
#SBATCH --time=0-5:0
#SBATCH --mail-user=sergey.korolkov@cosmos.ru
#SBATCH --mail-type=ALL
#SBATCH --constraint="type_e"
#SBATCH --partition=normal

# Executable
srun ./Veles2
