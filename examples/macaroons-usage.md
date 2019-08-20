### Integration with scalable external storage 


We will be using a genomics pipeline example to demonstrate the integration capabilities of Spider with highly scalable external storage systems.

The data we are going to use is part of a long-term evolution [experiment](https://en.wikipedia.org/wiki/E._coli_long-term_evolution_experiment) 
led by Richard Lenski to assess adaptation in E. coli. A population was propagated for more than 50,000 
generations in a glucose-limited minimal medium. We will be working with three sample events from the Ara-3 strain of this 
experiment, one from 5,000 generations, one from 15,000 generations, and one from 50,000 generations to study how the 
population changed. Generally, the quality of raw data is assessed and data is 'trimmed'. In this example, you will download 
a small set of data that has already been trimmed and will run the variant calling workflow.  

Let us download the scripts that will run the job for you:

```sh
cd $HOME
mkdir ecoli-analysis-dcache
cd ecoli-analysis-dcache
wget https://raw.githubusercontent.com/sara-nl/spider4advisors/master/examples/job-submit-variant-calling-dcache.sh

wget https://raw.githubusercontent.com/sara-nl/spider4advisors/master/examples/run-variant-calling-dcache.sh
```

Let us inspect the script that submits the job

```sh
cat job-submit-variant-calling-dcache.sh

#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake

mkdir "$TMPDIR"/var-calling
cd "$TMPDIR"/var-calling

cp $HOME/ecoli-analysis-dcache/run-variant-calling-dcache.sh .

export PATH="/project/surfadvisors/Software/ecoli-analysis-software/miniconda2/bin:$PATH"

bash run-variant-calling-dcache.sh 
```

Here we first created a directory with the help of a globally defined variable $TMPDIR. This directory will be created at 
the start of the job on the local scratch space and removed when the job is done. We copy the variant calling script to this
directory and run it. Let us first run the job and while it runs we can inspect how the data transfer happens within the job.

```sh
sbatch --job-name=var-call-dcache -J 'var-call-dcache' --output=%x-%j.out job-submit-variant-calling-dcache.sh
```

Now lets inspect the variant calling script
```sh
cat run-variant-calling-dcache.sh

#!/bin/bash
set -e
set -x
ecolipath=$PWD

mkdir -p data/ref_genome
cp /project/surfadvisors/Data/ecoli-analysis/data/ref_genome/ecoli_rel606.fasta data/ref_genome/

cd data/
curl https://webdav.grid.surfsara.nl:2880/?authz=MDAxY2xvY2F0aW9uIE9wdGlvbmFsLmVtcHR5CjAwMThpZGVudGlmaWVyIDVMdFI5S29QCjAwMzJjaWQgaWQ6NDM2MzI7NDEzODUsNDQ0MzYsNDI1MjksMzAwMTM7bWFpdGhpbGsKMDAyOGNpZCBiZWZvcmU6MjAxOS0wOS0xMlQxMDoxMzoyNy42NzVaCjAwNWFjaWQgcm9vdDovcG5mcy9ncmlkLnNhcmEubmwvZGF0YS9sc2dyaWQvU1VSRnNhcmEvc3BpZGVyY291cnNlL3RyaW1tZWRfZmFzdHFfc21hbGwudGFyCjAwMWZjaWQgYWN0aXZpdHk6RE9XTkxPQUQsTElTVAowMDJmc2lnbmF0dXJlIGL5MfchTf7sH1Ela025OBtIiYmsB3LAbutPyTbgW73yCg --output trimmed_fastq_small.tar
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

cp -r $TMPDIR/var-calling/results $HOME/ecoli-analysis-dcache/
```

> **_Food for brain:_**
>
> * The https link from where you download the data looks rather funny. Is this a normal URL? If not, do you know what it is?
> * Is this data freely available to anyone? Try copying the link in a browser on your laptop and see what happens.
> * Dosen't look like it requires any authentication. What if your data cannot be publicly made available?

The Ecoli probably do not mind their data being public but we are all very aware of data privacy - and this is not the case
only for genomic data
but research data in most domains. So how was the authentication performed for your input data? 

The data was shared with you with Macaroons - these are bearer tokens that you can use to authorize someone to
dwonload/upload/delete data stored on dCache. These macaroons can be used with clients that can support bearer tokens
(e.g., curl, Rclone). For this exercise a macaroon was created with certain restrictions (called as caveats) on the 
lifetime of the macaroon, the IP address you can use the macaroon from, the file that you can access, etc. Depending 
on who you want to share the data with, for how long and from which systems, these caveats can be adjusted. 

You can check the status of the job and inspect the output log file (even if the job is not completed).

```sh
squeue -u $USER
cat var-call-dcache-jobid.out #replace the jobid with your jobid 

#Another check would be the output of the following command
grep -v "#" $HOME/ecoli-analysis-dcache/results/vcf/SRR2589044_final_variants.vcf | wc -l

#The answer should be 10 (the number if expected variants detected in this population)
```

You can download the input data on the fly (if you have good network connectivity to the storage system from Spider) 
within each job without the hassle of downloading all the data to Spider. This is particularly handy if you want to 
automate large scale data analysis. In this example we saved the results locally on Spider, but you can also push the
output to dCache or another external storage system. dCache also supports username/password authentication and certificate 
based proxy authentication.   
