#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake

mkdir "$TMPDIR"/var-calling
cd "$TMPDIR"/var-calling

cp $HOME/ecoli-analysis-ada/run-variant-calling-ada.sh .

export PATH="/project/surfadvisors/Software/ecoli-analysis-software/miniconda2/bin:$PATH"

bash run-variant-calling-ada.sh 
