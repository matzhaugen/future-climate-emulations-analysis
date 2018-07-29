#!/bin/bash
#SBATCH --partition=bigmem2
#SBATCH --ntasks-per-node=24
#SBATCH --mem-per-cpu=6000
#SBATCH --output=multi_main_jack_bigmem.out
#SBATCH --time=36:00:00
module purge
module load R/3.4.3

Rscript multiMain_jack.R
