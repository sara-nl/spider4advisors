##  Collaboration within your project

We will be using a genomics pipeline example to test some of the collaboration functionalities of Spider. 

The data we are 
going to use is part of a long-term evolution experiment led by [Richard Lenski](https://en.wikipedia.org/wiki/E._coli_long-term_evolution_experiment)
to assess adaptation in E. coli. A population was propagated for more than 50,000 
generations in a glucose-limited minimal medium. We will be working with three sample events from the Ara-3 strain of this 
experiment, one from 5,000 generations, one from 15,000 generations, and one from 50,000 generations to study how the 
population changed. Generally, the quality of raw data is assessed and data is 'trimmed'. In this example, you will download 
a small set of data that has already been trimmed and will run the variant calling workflow.

Let us download the scripts that will run the job for you - one script for setting the job environment that calls the script that runs the analysis.

```sh
cd $HOME
mkdir ecoli-analysis-cephfs
cd ecoli-analysis-cephfs
wget https://raw.githubusercontent.com/sara-nl/spidercourse/master/scripts/job-submit-variant-calling-cephfs-adv.sh

wget https://raw.githubusercontent.com/sara-nl/spidercourse/master/scripts/run-variant-calling-cephfs-adv.sh
```
Let us inspect what theese scripts will do:

```sh
cat job-submit-variant-calling-cephfs-adv.sh

#!/bin/bash
#SBATCH -c 1
#SBATCH --constraint=skylake

export PATH="/project/surfadvisors/Software/ecoli-analysis-software/miniconda2/bin:$PATH"

bash $HOME/ecoli-analysis-cephfs/run-variant-calling-cephfs-adv.sh 
```

The #SBATCH flags that you see in the script have the following function:

-c: 1 core requested
--constraint: here we request nodes with skylake processors 

We then set up the software environment for your analysis and then run the script for data analysis.

> **_Food for brain:_**
>
> * Nice to have the software ready to use, but who installed it? And do I have the permission to use this software?
> * The software was installed by the 'Software manager' of the surfadvisors project. Do you know who this angel is ;)? Hint: getent group surfadvisors-sw
> * Do you know if you are a Software manager? If not, do you know what 'role' you have in the project?

Let us submit the job first and then inspect the steps while the job runs

```sh
sbatch --job-name=var-call -J 'var-call' --output=%x-%j.out job-submit-variant-calling-cephfs-adv.sh
squeue -u $USER

cat run-variant-calling-cephfs-adv.sh

#!/bin/bash
set -e
ecolipath=$HOME/ecoli-analysis-cephfs

mkdir -p data/ref_genome
cp /project/surfadvisors/Data/ecoli-analysis/data/ref_genome/ecoli_rel606.fasta data/ref_genome/
ls data/ref_genome

mkdir -p data/trimmed_fastq_small
cp /project/surfadvisors/Data/ecoli-analysis/data/trimmed_fastq_small/*fastq data/trimmed_fastq_small/
ls data/trimmed_fastq_small

mkdir -p $ecolipath/results
cd $ecolipath/results

genome=$ecolipath/data/ref_genome/ecoli_rel606.fasta

#index the reference genome
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
    
    #Align reads to reference genome
    bwa mem $genome $fq1 $fq2 > $sam
    
    #Convert from sam to bam format
    samtools view -S -b $sam > $bam

    #Sort the bam files    
    samtools sort -o $sorted_bam $bam 
    
    #Calculate the read coverage of positions in the genome
    bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam
    
    #Detect the single nucleotide polymorphisms (SNPs)
    bcftools call --ploidy 1 -m -v -o $variants $raw_bcf 
    
    #Filter the SNPs for the final output in VCF format
    vcfutils.pl varFilter $variants > $final_variants
   
    done
```
Here we copy the input data files first to your working directory from the Data space of your project. The output is directly saved in your home directory. Naturally you can also read the input data directly from the Data project space. 

> **_Food for brain:_**
>
> * Nice to have the raw data and reference genome ready to use, but who downloaded it? 
> * The data was downloaded by the 'Data manager' of the surfadvisors project. Do you know who this angel is ;)? Hint: getent group surfadvisors-data
> * How would you go about sharing your results with your collagues in the project? The results are located in your home directory that no one else has access to. Hint: ls /project

Let us see if the job is ready for you to inspect the output

```sh
squeue -u $USER
```

If your job has finished you can check the output log for errors

```sh
cat var-call-tmpdir-jobID.out #replace the jobID

#Another check would be the output of the following command
grep -v "#" $HOME/ecoli-analysis-cephfs/results/vcf/SRR2589044_final_variants.vcf | wc -l

#The answer should be 10 (the number if expected variants detected in this population)
```

In this example you ran analysis on data that was downloaded by Data manager in the Data project space, using Software that was installed by a Software manager and looked into the possibility of sharing the results with your colleagues in your project. Pretty cool collaboration features right? 

> **_Food for brain:_**
>
> * What if you want to share results with collaborators who are not a part of the project at all? Hint: ls /project. 
> * How about opening a browser on your laptop and going to this site -  https://public.spider.surfsara.nl/project/spidercourse/ecoli-analysis/.  Whose project space is this? This is the 'untrimmed' raw data that you used for your analysis.
> * Can you accidentally remove the raw data from your project space Data folder?
> * Can you accidentally remove the Software from your project space Software folder?
> * Can you accidentally delete files from your project space Share/Public folder? Warning: Be nice now - collaborate with someone before going all terminator here!




