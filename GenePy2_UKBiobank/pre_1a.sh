#!/bin/bash
#SBATCH --nodes=1
#SBATCH --mem=16g
#SBATCH --job-name="carrot"
#SBATCH --ntasks-per-node=1
#SBATCH --time=15:59:00

##pre_2 extract the indel and complex variants which are often not prescored to add the CADD score to the database


cd /home/gc1a20/hgig_me/ukbb/wes_200k/run_202206

##variants extraction from the output vcf of pre_1

awk -F"\t" '$1 ~/#/ || length($4)>1||length($5)>1' $1 | sed '3383,$s/chr//g' $1 > $2


##CADD annoation of the variants


~/pirate/software/CADD-scripts/CADD.sh \
    -o $3 $2

gunzip $3
rm $2


##concate the scores and save to the CADD pre-score database
# head -2 original/ukb23156_chr14_31113681_35469448_v1_p1_noCADD.tsv >  ~/pirate/software/CADD-scripts/data/prescored/GRCh38_v1.6/no_anno/ukbb23156_complex_SNP_INDELS.tsv
# cat original/ukb23156*noCADD.tsv | grep -v '#' >> ~/pirate/software/CADD-scripts/data/prescored/GRCh38_v1.6/no_anno/ukbb23156_complex_SNP_INDELS.tsv
