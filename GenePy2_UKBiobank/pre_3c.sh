#!/bin/bash
#SBATCH --mem=16g
#SBATCH --nodes=1
#SBATCH --job-name="Yrabbit"
#SBATCH --ntasks-per-node=1
#SBATCH --time=23:59:00

##pre_3c is the filtration step of pre_3*.sh
##modified on July 12th 2023 to include all variants rather than exonic-only variants; random-forest based filtration is also defaulted in this modification due to the lack of evidence in qc (flag for bench marking)


cd /home/gc1a20/hgig_me/ukbb/wes_200k/run_202206/original

gunzip -c $1 | grep -v '##' |cut -f 9-> $2

cp ~/scratch/header $3
grep -v '##' $4 >$5

##the change
#paste $5 $2 | grep -v 'FAIL' | fgrep -w -f /home/gc1a20/runukbb/script/GenePy-2/exonic_variants >> $3
paste $5 $2 >> $3

rm $2 $5
### QC-based filtration modified with addtion of missingness based filtration to the original DP/AB-based filtration
##module load biobuilds
#bgzip $3
#tabix -p vcf $3.gz
#module load picard
#java -Xmx12g -jar /local/software/picard/2.18.14/picard.jar \
#    FilterVcf \
#    I=$3.gz O=$6 MIN_DP=8 MIN_AB=0.15

bcftools +fill-tags $3 -- -t 'FORMAT/AB:1=float((AD[:1]) / (DP))' | bgzip -c > $3.gz
rm $3
tabix -p vcf $3.gz
bcftools filter -S . --include 'FORMAT/DP>=8 & FORMAT/AB>=0.15 |FORMAT/GT="0/0"' -Oz -o $6.gz $3.gz
tabix -p vcf $6.gz
bcftools +fill-tags $6.gz -- -S ~/ref/ethnicity.txt -t 'HWE,F_MISSING' | bcftools view -e '(CHROM=="chrY" & INFO/F_MISSING>=0.56 & INFO/HWE_1>(0.05/15922704))'| bcftools view -i 'INFO/F_MISSING<0.12 & INFO/HWE_1>(0.05/15922704)' -Oz -o $7.gz

tabix -p vcf $7.gz

bcftools view $7.gz -R ~/ref/target_ref/xgen_plus_spikein.GRCh38.bed -Oz -o $8
rm $3.gz* $6.gz* $7.gz*

