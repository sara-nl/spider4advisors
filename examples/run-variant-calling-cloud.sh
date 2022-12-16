#!/bin/bash
set -x
set -e
ecolipath=$HOME/ecoli-analysis-cloud
apptainerpath=/cvmfs/softdrive.nl/anatolid/spidercourse/

cd $ecolipath
mkdir -p data/ref_genome
curl https://webdav.grid.surfsara.nl:2880/?authz=MDAxY2xvY2F0aW9uIE9wdGlvbmFsLmVtcHR5CjAwMThpZGVudGlmaWVyIDNoN1hqZjcvCjAwMTVjaWQgaWlkOkdhYmdHTWFRCjAwMWRjaWQgaWQ6MzY0OTg7MzE4Nzk7cHZpZXIKMDAyYmNpZCBiZWZvcmU6MjAyMy0xMi0xNlQxNToyNjowMC44MzA5OTFaCjAwNTBjaWQgcm9vdDovcG5mcy9ncmlkLnNhcmEubmwvZGF0YS9wdmllci9TVVJGL3NwaWRlcmNvdXJzZS9lY29saV9yZWw2MDYuZmFzdGEKMDAxZmNpZCBhY3Rpdml0eTpET1dOTE9BRCxMSVNUCjAwYmNjaWQgaXA6MTQ1LjM4LjAuMC8xNiwxNDUuMTAwLjUuMC8yNywxNDUuMTAwLjUuMjEwLzI2LDE0NS4xMDAuMzIuMC8yMiwxNDUuMTAwLjQ4LjAvMjMsMTQ1LjEwMC41MC4wLzIzLDE0NS4xMDAuMjAwLjAvMjEsMTQ1LjEwMC45LjY0LzI5LDE0NS4xMDEuMzIuMC8yMSwxNDUuMTAwLjU2LjAvMjIsMjAwMTo2MTA6MTA4OjovNDgKMDAyZnNpZ25hdHVyZSCUYMBQEHcDbOIedn-szW1XziGGwj1qmpY2c_zYdduzRgo -o data/ref_genome/ecoli_rel606.fasta
ls data/ref_genome

cd $ecolipath/data/
curl https://webdav.grid.surfsara.nl:2880/?authz=MDAxY2xvY2F0aW9uIE9wdGlvbmFsLmVtcHR5CjAwMThpZGVudGlmaWVyIDNoN1hqZjcvCjAwMTVjaWQgaWlkOmFxUFRha1dqCjAwMWRjaWQgaWQ6MzY0OTg7MzE4Nzk7cHZpZXIKMDAyYmNpZCBiZWZvcmU6MjAyMy0xMi0wOVQxNjowNDo1NC45NTUxMjFaCjAwNTVjaWQgcm9vdDovcG5mcy9ncmlkLnNhcmEubmwvZGF0YS9wdmllci9TVVJGL3NwaWRlcmNvdXJzZS90cmltbWVkX2Zhc3RxX3NtYWxsLnRhcgowMDFmY2lkIGFjdGl2aXR5OkRPV05MT0FELExJU1QKMDBiY2NpZCBpcDoxNDUuMzguMC4wLzE2LDE0NS4xMDAuNS4wLzI3LDE0NS4xMDAuNS4yMTAvMjYsMTQ1LjEwMC4zMi4wLzIyLDE0NS4xMDAuNDguMC8yMywxNDUuMTAwLjUwLjAvMjMsMTQ1LjEwMC4yMDAuMC8yMSwxNDUuMTAwLjkuNjQvMjksMTQ1LjEwMS4zMi4wLzIxLDE0NS4xMDAuNTYuMC8yMiwyMDAxOjYxMDoxMDg6Oi80OAowMDJmc2lnbmF0dXJlIIrtFFFlJs6jjhPcHMd2LgodarueaalJJsxcoS5O89zYCg --output trimmed_fastq_small.tar
tar xvf trimmed_fastq_small.tar

mkdir $ecolipath/results
cd $ecolipath/results

genome=$ecolipath/data/ref_genome/ecoli_rel606.fasta

apptainer exec $apptainerpath/elixir-singularity.sif bwa index $genome

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
    apptainer exec $apptainerpath/elixir-singularity.sif bwa mem $genome $fq1 $fq2 > $sam

    #Convert from sam to bam format
    apptainer exec $apptainerpath/elixir-singularity.sif samtools view -S -b $sam > $bam

    #Sort the bam files
    apptainer exec $apptainerpath/elixir-singularity.sif samtools sort -o $sorted_bam $bam

    #Calculate the read coverage of positions in the genome
    apptainer exec $apptainerpath/elixir-singularity.sif bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam

    #Detect the single nucleotide polymorphisms (SNPs)
    apptainer exec $apptainerpath/elixir-singularity.sif bcftools call --ploidy 1 -m -v -o $variants $raw_bcf

    #Filter the SNPs for the final output in VCF format
    apptainer exec $apptainerpath/elixir-singularity.sif vcfutils.pl varFilter $variants > $final_variants

done
