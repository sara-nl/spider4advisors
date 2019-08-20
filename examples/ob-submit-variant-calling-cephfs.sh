#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake

export PATH="/project/surfadvisors/Software/ecoli-analysis-software/miniconda2/bin:$PATH"

bash $HOME/ecoli-analysis-cephfs/run-variant-calling-cephfs.sh 
