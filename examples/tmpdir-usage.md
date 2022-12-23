## High throughput data processing model 

We will be using a genomics pipeline example to test some of the high-throughput functionalities of Spider.

The data we are going 
to use is part of a long-term evolution [experiment](https://en.wikipedia.org/wiki/E._coli_long-term_evolution_experiment) 
led by Richard Lenski to assess adaptation in E. coli. A population was propagated for more than 50,000 
generations in a glucose-limited minimal medium. We will be working with three sample events from the Ara-3 strain of this 
experiment, one from 5,000 generations, one from 15,000 generations, and one from 50,000 generations to study how the 
population changed. Generally, the quality of raw data is assessed and data is 'trimmed'. In this example, you will download 
a small set of data that has already been trimmed and will run the variant calling workflow.

Let us download the scripts that will run the job for you:

```sh
cd $HOME
mkdir ecoli-analysis-tmpdir
cd ecoli-analysis-tmpdir
wget https://raw.githubusercontent.com/sara-nl/spider4advisors/master/examples/job-submit-variant-calling-tmpdir.sh

wget https://raw.githubusercontent.com/sara-nl/spider4advisors/master/examples/run-variant-calling-tmpdir.sh
```
We copy the files and scripts to the local scratch space on the worker node where your job lands. Let us inspect these scripts.

```sh
cat job-submit-variant-calling-tmpdir.sh

#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake

mkdir "$TMPDIR"/var-calling
cd "$TMPDIR"/var-calling

cp $HOME/ecoli-analysis-tmpdir/run-variant-calling-tmpdir.sh .

export PATH="/project/surfadvisors/Software/ecoli-analysis-software/miniconda3/bin:$PATH"

time bash run-variant-calling-tmpdir.sh 

```
Here we first created a directory with the help of a globally defined variable $TMPDIR. This directory will be 
created at the start of the job on the local scratch space and removed when the job is done. We copy the variant calling 
script to this directory and run it. To compare the performance with jobs that ran with data located on the project spaces,
we will 'time' the job - this will tell us how long it took for the full job to finish.

Let us submit the job first and then inspect the steps while the job runs

```sh
sbatch --job-name=var-call -J 'var-call' --output=%x-%j.out job-submit-variant-calling-tmpdir.sh
squeue -u $USER

cat run-variant-calling-tmpdir.sh

#!/bin/bash
set -e
set -x
ecolipath=$PWD

mkdir -p data/ref_genome
cp /project/surfadvisors/Data/ecoli-analysis/data/ref_genome/ecoli_rel606.fasta data/ref_genome/
ls data/ref_genome

mkdir data/trimmed_fastq_small
cp /project/surfadvisors/Data/ecoli-analysis/data/trimmed_fastq_small/*fastq data/trimmed_fastq_small/
ls data/trimmed_fastq_small

mkdir results
cd $ecolipath/results

genome=$ecolipath/data/ref_genome/ecoli_rel606.fasta

bwa index $genome

mkdir -p sam bam bcf vcf

for fq1 in $ecolipath/data/trimmed_fastq_small/*_1.trim.sub.fastq
    do
    echo "working with file $fq1"

    base=$(basename $fq1 _1.trim.sub.fastq)
    echo "base name is $base"

    fq1=$ecolipath/data/trimmed_fastq_small/${base}_1.trim.sub.fastq
    fq2=$ecolipath/data/trimmed_fastq_small/${base}_2.trim.sub.fastq
    sam=$ecolipath/results/sam/${base}.aligned.sam
    bam=$ecolipath/results/bam/${base}.aligned.bam
    sorted_bam=$ecolipath/results/bam/${base}.aligned.sorted.bam
    raw_bcf=$ecolipath/results/bcf/${base}_raw.bcf
    variants=$ecolipath/results/bcf/${base}_variants.vcf
    final_variants=$ecolipath/results/vcf/${base}_final_variants.vcf 

    bwa mem $genome $fq1 $fq2 > $sam
    samtools view -S -b $sam > $bam
    samtools sort -o $sorted_bam $bam 
    samtools index $sorted_bam
    bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam
    bcftools call --ploidy 1 -m -v -o $variants $raw_bcf 
    vcfutils.pl varFilter $variants > $final_variants
   
    done

cp -r $TMPDIR/var-calling/results $HOME/ecoli-analysis-tmpdir/
```
Here we copy the input data to the $TMPDIR. In the end we copy the output to our $HOME directory as the $TMPDIR is removed after the job finishes and we will lose our results. You can compare if the performance was better/worse/equivalent to the performance with the jobs when the data is in project spaces.

> **_Food for brain:_**
>
> * What does the time command do? How do you interpret the output?
> * You need to rerun the previous example with data in the project space by adding the 'time' command.
> * Does the $TMPDIR example have better performance? When is it advantageous to use it?

To test the status of your job you can run the command

```sh
squeue -u $USER
```

If your job has finished you can check the output log for errors

```sh
cat var-call-tmpdir-jobID.out #replace the jobID

#Another check would be the output of the following command
grep -v "#" $HOME/ecoli-analysis-tmpdir/results/vcf/SRR2589044_final_variants.vcf | wc -l

#The answer should be 10 (the number if expected variants detected in this population)
```

In this example you used the SSD available locally on the worker nodes instead of having your data on the shared 
project spaces. For large datasets and heavy processing (particularly heavy I/O), the overall gain in using the scratch
space can be higher.
