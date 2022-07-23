#!/bin/bash

#SBATCH --job-name=merge
#SBATCH --mail-user=marcase@umich.edu
#SBATCH --mail-type=END,FAIL
#SBATCH --cpus-per-task=1
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-cpu=128000m
#SBATCH --time=16:00:00
#SBATCH --account=gthurber0
#SBATCH --partition=standard
#SBATCH --output=/home/%u/%x-%j-merge.log

# The application(s) to execute along with its input arguments and options:
module load python3.8-anaconda/2021.05
python NGS_merge.py
