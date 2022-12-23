#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake

mkdir "$TMPDIR"/var-calling
cd "$TMPDIR"/var-calling

cp $HOME/ecoli-analysis-tmpdir/run-variant-calling-tmpdir.sh .

export PATH="/project/surfadvisors/Software/ecoli-analysis-software/miniconda3/bin:$PATH"

time bash run-variant-calling-tmpdir.sh 
