#!/bin/bash
set -e
set -x

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
