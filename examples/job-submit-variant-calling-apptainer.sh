#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake

bash $HOME/ecoli-analysis-container/run-variant-calling-apptainer.sh 
