#!/bin/bash
#SBATCH --mem=8g
#SBATCH --nodes=1
#SBATCH --job-name="headstar"
#SBATCH --ntasks-per-node=1
#SBATCH --time=02:59:00

##pre_1 split the vcf and omit the sample info to speed up the annotation process

module load biobuilds

cd /home/gc1a20/hgig_me/ukbb/wes_200k/run_202206/original
#tabix -f -p vcf $1
bcftools view -G $1 -Ov -o $2


module load ensembl-vep


vep -i $2 \
    --assembly GRCh38 \
    --vcf \
    --fork 1 \
    --cache --force_overwrite \
    --pick_allele \
    --plugin CADD,/home/gc1a20/pirate/software/CADD-scripts/data/prescored/GRCh38_v1.6/no_anno/gnomad.genomes.r3.0.indel.tsv.gz,/home/gc1a20/pirate/software/CADD-scripts/data/prescored/GRCh38_v1.6/no_anno/ukbb23156_complex_SNP_INDELS.tsv.gz,/home/gc1a20/pirate/software/CADD-scripts/data/prescored/GRCh38_v1.6/no_anno/whole_genome_SNVs.tsv.gz \
    --custom /mainfs/hgig/private/software/gnomAD/GRCh38/gnomad.genomes.v3.1.1.RF_flag.vcf.gz,gnomadRF,vcf,exact,,RF_flag \
    --fields "Allele,Consequence,SYMBOL,Gene,gnomadE_AF,CADD_RAW,gnomadRF_RF_flag" \
    -o $3 \
    --offline



