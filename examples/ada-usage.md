### The Advanced dCache API (ADA)


We will reuse the genomics pipeline from the [dCache example](#integration-with-scalable-external-storage) to demonstrate the capabilities of ADA on Spider.

The data we are going to use is part of a long-term evolution [experiment](https://en.wikipedia.org/wiki/E._coli_long-term_evolution_experiment)
led by Richard Lenski to assess adaptation in E. coli. A population was propagated for more than 50,000
generations in a glucose-limited minimal medium. We will be working with three sample events from the Ara-3 strain of this
experiment, one from 5,000 generations, one from 15,000 generations, and one from 50,000 generations to study how the
population changed. Generally, the quality of raw data is assessed and data is 'trimmed'. In this example, you will download
a small set of data that has already been trimmed and will run the variant calling workflow.  

Let us download the scripts that will run the job for you:

```sh
cd $HOME
mkdir ecoli-analysis-ada
cd ecoli-analysis-ada
wget https://raw.githubusercontent.com/sara-nl/spider4advisors/master/examples/job-submit-variant-calling-ada.sh

wget https://raw.githubusercontent.com/sara-nl/spider4advisors/master/examples/run-variant-calling-ada.sh
```

Before we run the script we also need a macaroon that will allow us to approach dCache through the ADA interface. 
Luckily for us, the macaroon is already available at `/project/surfadvisors/Share/maca_spider4advisors_ada.conf` and can be read by
SURF advisors following this course. This macaroon can only read, (un)stage and download files, to ensure the data is not accidentally
deleted, among other security features.

First, let's see if we can list the file that we want to download in the script:

```sh
ada --tokenfile /project/surfadvisors/Share/maca_spider4advisors_ada.conf --longlist ./trimmed_fastq_small.tar

./trimmed_fastq_small.tar  347176960  2022-12-08 14:15 UTC  tape  NEARLINE
```

With `--longlist` we see more information than with the default `--list`. What we see is that this file is currently 
stored and is in the _nearline_ status. This means that the file is directly available for reading, just that a disk
needs to be spinned up, which will only take a few seconds. 

The script we are about to run actually _stages_ the file for us, meaning that it is taken to the _online_ state, which means it is
immediately available for reading. After the file is downloaded to local storage it is also _unstaged_, putting it back in the nearline (disk)
or offline (tape) state. Unstaging happens automatically after some time, but in this example the script explicitly does it to show the user it can be done by hand.

Let us inspect the script that submits the job:

```sh
cat job-submit-variant-calling-ada.sh

#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake

mkdir "$TMPDIR"/var-calling
cd "$TMPDIR"/var-calling

cp $HOME/ecoli-analysis-ada/run-variant-calling-ada.sh .

export PATH="/project/surfadvisors/Software/ecoli-analysis-software/miniconda2/bin:$PATH"

bash run-variant-calling-ada.sh
```

Here we first created a directory with the help of a globally defined variable $TMPDIR. This directory will be created at
the start of the job on the local scratch space and removed when the job is done. We copy the variant calling script to this
directory and run it. Let us first run the job and while it runs we can inspect how the data transfer happens within the job.

```sh
sbatch --job-name=var-call-ada -J 'var-call-ada' --output=%x-%j.out job-submit-variant-calling-ada.sh
```

Now lets inspect the variant calling script:

```sh
cat run-variant-calling-ada.sh

#!/bin/bash
set -e
set -x
ecolipath=$PWD

mkdir -p data/ref_genome
cp /project/surfadvisors/Data/ecoli-analysis/data/ref_genome/ecoli_rel606.fasta data/ref_genome/

cd data/

ada --tokenfile /project/surfadvisors/Share/maca_spider4advisors_ada.conf --stage ./trimmed_fastq_small.tar
ada --tokenfile /project/surfadvisors/Share/maca_spider4advisors_ada.conf --longlist ./trimmed_fastq_small.tar
rclone --config=/project/surfadvisors/Share/maca_spider4advisors_ada.conf copy maca_spider4advisors_ada:./trimmed_fastq_small.tar . -P
ada --tokenfile /project/surfadvisors/Share/maca_spider4advisors_ada.conf --unstage ./trimmed_fastq_small.tar
ada --tokenfile /project/surfadvisors/Share/maca_spider4advisors_ada.conf --longlist ./trimmed_fastq_small.tar

tar xvf trimmed_fastq_small.tar

mkdir $ecolipath/results
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

cp -r $TMPDIR/var-calling/results $HOME/ecoli-analysis-ada/
```

As you can see `ada` takes care of the file staging and unstaging, while `rclone` is used to move the file to local storage.

> **_Food for brain:_**
>
> * The ada commands explicitly show the file that is downloaded. Can you find other files accessible with this macaroon?
> * Is this data freely available to anyone? Try running the ada list command outside of the SURF network (not on Spider) and see what happens.

You can check the status of the job and inspect the output log file (even if the job is not completed).

```sh
squeue -u $USER
cat var-call-ada-jobid.out #replace the jobid with your jobid

#Another check would be the output of the following command
grep -v "#" $HOME/ecoli-analysis-ada/results/vcf/SRR2589044_final_variants.vcf | wc -l

#The answer should be 10 (the number if expected variants detected in this population)
```
