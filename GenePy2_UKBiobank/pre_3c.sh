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
module load biobuilds

bgzip $3
tabix -p vcf $3.gz

module load picard

java -Xmx12g -jar /local/software/picard/2.18.14/picard.jar \
    FilterVcf \
    I=$3.gz O=$6 MIN_DP=8 MIN_AB=0.15

bgzip $6
tabix -p vcf $6.gz

#bcftools view $6.gz -R ~/ref/target_ref/xgen_plus_spikein.GRCh38.bed -Oz -o $7
bcftools view $6.gz -R ~/ref/target_ref/xgen_plus_spikein.GRCh38.bed -Ov -o filtered.vcf
uniq filtered.vcf | bgzip -c > $7

rm $3.gz $6.gz filtered.vcf


