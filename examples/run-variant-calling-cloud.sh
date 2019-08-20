#!/bin/bash
set -x
set -e
ecolipath=$HOME/ecoli-analysis-cloud
singularitypath=/cvmfs/softdrive.nl/anatolid/spidercourse/

cd $ecolipath
mkdir -p data/ref_genome
curl https://webdav.grid.surfsara.nl:2880/?authz=MDAxY2xvY2F0aW9uIE9wdGlvbmFsLmVtcHR5CjAwMThpZGVudGlmaWVyIDVMdFI5S29QCjAwMzJjaWQgaWQ6NDM2MzI7NDEzODUsNDQ0MzYsNDI1MjksMzAwMTM7bWFpdGhpbGsKMDAyOGNpZCBiZWZvcmU6MjAxOS0wOS0xNFQxNTowMDoxNS44MzNaCjAwNTVjaWQgcm9vdDovcG5mcy9ncmlkLnNhcmEubmwvZGF0YS9sc2dyaWQvU1VSRnNhcmEvc3BpZGVyY291cnNlL2Vjb2xpX3JlbDYwNi5mYXN0YQowMDFmY2lkIGFjdGl2aXR5OkRPV05MT0FELExJU1QKMDAyZnNpZ25hdHVyZSDzlhQyMCATd3ItdfCQxu4fIEwwq8iRBWBc_3NxPyHbbQo -o data/ref_genome/ecoli_rel606.fasta
ls data/ref_genome

cd $ecolipath/data/
curl https://webdav.grid.surfsara.nl:2880/?authz=MDAxY2xvY2F0aW9uIE9wdGlvbmFsLmVtcHR5CjAwMThpZGVudGlmaWVyIDVMdFI5S29QCjAwMzJjaWQgaWQ6NDM2MzI7NDEzODUsNDQ0MzYsNDI1MjksMzAwMTM7bWFpdGhpbGsKMDAyOGNpZCBiZWZvcmU6MjAxOS0wOS0xMlQxMDoxMzoyNy42NzVaCjAwNWFjaWQgcm9vdDovcG5mcy9ncmlkLnNhcmEubmwvZGF0YS9sc2dyaWQvU1VSRnNhcmEvc3BpZGVyY291cnNlL3RyaW1tZWRfZmFzdHFfc21hbGwudGFyCjAwMWZjaWQgYWN0aXZpdHk6RE9XTkxPQUQsTElTVAowMDJmc2lnbmF0dXJlIGL5MfchTf7sH1Ela025OBtIiYmsB3LAbutPyTbgW73yCg --output trimmed_fastq_small.tar
tar xvf trimmed_fastq_small.tar

mkdir $ecolipath/results
cd $ecolipath/results

genome=$ecolipath/data/ref_genome/ecoli_rel606.fasta

singularity exec $singularitypath/elixir-singularity.sif bwa index $genome

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
    singularity exec $singularitypath/elixir-singularity.sif bwa mem $genome $fq1 $fq2 > $sam

    #Convert from sam to bam format
    singularity exec $singularitypath/elixir-singularity.sif samtools view -S -b $sam > $bam

    #Sort the bam files
    singularity exec $singularitypath/elixir-singularity.sif samtools sort -o $sorted_bam $bam

    #Calculate the read coverage of positions in the genome
    singularity exec $singularitypath/elixir-singularity.sif bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam

    #Detect the single nucleotide polymorphisms (SNPs)
    singularity exec $singularitypath/elixir-singularity.sif bcftools call --ploidy 1 -m -v -o $variants $raw_bcf

    #Filter the SNPs for the final output in VCF format
    singularity exec $singularitypath/elixir-singularity.sif vcfutils.pl varFilter $variants > $final_variants

done
