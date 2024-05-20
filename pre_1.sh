#!/bin/bash

module load biobuilds

cd ${wkdir}

bcftools view -G $1 -Ov -o p1.vcf
##
#awk -F"\t" '$1 ~/#/ || length($4)>1||length($5)>1' p1.vcf | sed '3408,$s/chr//g' > p11.vcf
##
##
##CADD annoation of the variants
#
#module unload biobuilds
conda activate cadd
#
~/pirate/software/CADD-scripts/CADD.sh \
    -c 8 \
    -o wes_2023010_patch.tsv.gz p11.vcf
#
#module load biobuilds
tabix -p vcf wes_2023010_patch.tsv.gz
#
mv wes_2023010_patch* ~/pirate/software/CADD-scripts/data/prescored/GRCh38_v1.6/no_anno/
#
#
module load ensembl-vep
#
#
vep -i p1.vcf \
    --offline \
    --assembly GRCh38 \
    --vcf \
    --fork 10 \
    --cache --force_overwrite \
    --pick_allele \
    --plugin CADD,/home/gc1a20/pirate/software/CADD-scripts/data/prescored/GRCh38_v1.6/no_anno/wes_2023010_patch.tsv.gz,/home/gc1a20/pirate/software/CADD-scripts/data/prescored/GRCh38_v1.6/no_anno/whole_genome_SNVs.tsv.gz,/home/gc1a20/pirate/software/CADD-scripts/data/prescored/GRCh38_v1.6/no_anno/gnomad.genomes.r3.0.indel.tsv.gz \
    --custom /mainfs/hgig/private/software/gnomAD/GRCh38/gnomad.genomes.v3.1.1.RF_flag.vcf.gz,gnomadRF,vcf,exact,,RF_flag \
    --af_gnomade --af_gnomadg \
    --fields "Allele,Consequence,SYMBOL,Gene,gnomADg_AF,gnomADg_NFE_AF,gnomADe_AF,gnomADe_NFE_AF,CADD_RAW,gnomadRF_RF_flag" \
    -o p1.vep.vcf

gunzip -c $1 | grep -v '##' |cut -f 9-> p2
grep -v '##' p1.vep.vcf >p1
grep '##' p1.vep.vcf >p12.vcf
paste p1 p2 >> p12.vcf
rm p1 p2


##filtrtion for genepy 20230420
module load biobuilds

awk -F"\t" '$7~/PASS/ || $1~/#/' p12.vcf >f1.vcf
vcftools --vcf f1.vcf --minDP 8 --minGQ 20 --recode-INFO-all --recode --out f2.vcf

~/bin/meanGQ_filter.sh f2.vcf.recode.vcf f3.vcf

vcftools --vcf f3.vcf --max-missing 0.88 --recode-INFO-all --recode --out f4
vcftools --vcf f3.vcf --chr chrY --max-missing 0.44 --recode-INFO-all --recode --out f4_y
grep -v '#' f4_y.recode.vcf >> f4.recode.vcf
bgzip f4.recode.vcf
tabix -p vcf f4.recode.vcf.gz

bcftools view f4.recode.vcf.gz -R /home/gc1a20/ref/target_ref/SSV5_6_pad100_intersect_uniq.bed -Ov -o f5.vcf


awk -F"\t" '{OFS=FS}{gsub(/,.*$/, "", $5)}1' f5.vcf \
    | awk '$5 !~/\*/' | \
    awk -F"\t" '{OFS=FS}{for (i=10;i<=NF;i++) if ($i ~/\/2|\/3/) $i="./."; else $i=substr($i,1,3)}1' | \
    grep -v 'FAIL' | \
    fgrep -w -f /home/gc1a20/runukbb/script/GenePy-2/exonic_variants >> f6.vcf

module load picard

java -Xmx240G -jar /local/software/picard/2.18.14/picard.jar \
    SortVcf I=f5.vcf O=x.vcf SD=/mainfs/hgig/public/HUMAN_REFS/HG38/REF_HLA/GRCh38_full_analysis_set_plus_decoy_hla.dict
bcftools norm -d both -Ov x.vcf -o x

grep '#' p12.vcf >f6.vcf
awk '$1 !~/#/' f5.vcf | \
    awk -F"\t" '{OFS=FS}{for (i=10;i<=NF;i++) if ($i ~/\.\//) $i="./."; else $i=substr($i,1,3)}1' | \
    awk -F"\t" '{OFS=FS}{for (i=10;i<=NF;i++) gsub(/\|/,"/",$i)}1' >>f6.vcf
#    grep -v 'FAIL | \
#    fgrep -w -f /home/gc1a20/runukbb/script/GenePy-2/exonic_variants >> f6.vcf
