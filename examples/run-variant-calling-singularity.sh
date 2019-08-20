#!/bin/bash
set -x
set -e
ecolipath=$HOME/ecoli-analysis-container

mkdir -p data/ref_genome
cp /project/surfadvisors/Data/ecoli-analysis/data/ref_genome/ecoli_rel606.fasta data/ref_genome/
ls data/ref_genome

mkdir data/trimmed_fastq_small
cp /project/surfadvisors/Data/ecoli-analysis/data/trimmed_fastq_small/*fastq data/trimmed_fastq_small/
ls data/trimmed_fastq_small

mkdir results
cd $ecolipath/results

genome=$ecolipath/data/ref_genome/ecoli_rel606.fasta

singularity exec $ecolipath/elixir-singularity.sif bwa index $genome

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
    singularity exec $ecolipath/elixir-singularity.sif bwa mem $genome $fq1 $fq2 > $sam
    
    #Convert from sam to bam format
    singularity exec $ecolipath/elixir-singularity.sif samtools view -S -b $sam > $bam

    #Sort the bam files    
    singularity exec $ecolipath/elixir-singularity.sif samtools sort -o $sorted_bam $bam 
    
    #Calculate the read coverage of positions in the genome
    singularity exec $ecolipath/elixir-singularity.sif bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam
    
    #Detect the single nucleotide polymorphisms (SNPs)
    singularity exec $ecolipath/elixir-singularity.sif bcftools call --ploidy 1 -m -v -o $variants $raw_bcf 
    
    #Filter the SNPs for the final output in VCF format
    singularity exec $ecolipath/elixir-singularity.sif vcfutils.pl varFilter $variants > $final_variants
   
done
