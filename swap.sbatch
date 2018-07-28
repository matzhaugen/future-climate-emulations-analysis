#!/bin/bash
#SBATCH --nodes=1
#SBATCH --partition=bigmem2
#SBATCH --ntasks-per-node=4
#SBATCH --output=main.out

# --ntasks-per-node will be used in doParallel.R to specify the number of 
# cores to use on the machine. Using 16 will allow us to use all cores
# on a sandyb node

module load R

Rscript swapNcdfDim.R
