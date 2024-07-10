#!/bin/bash
#SBATCH --mem=2g
#SBATCH --nodes=1
#SBATCH --job-name="genarch"
#SBATCH --ntasks-per-node=1
#SBATCH --time=11:59:00

##pre_4 transform the final vcf to the genepy input format, retaining only necessary information of the variants



##kggseq succinct the genotype field?


###transform the variant information
module load biobuilds

cd /home/gc1a20/hgig_me/ukbb/wes_200k/run_202206

DIR=$(echo $1 | sed 's/_v1_p12_filterNew_onTarget\.vcf\.gz//g')

mkdir -p ${DIR}
cd ${DIR}

cp ../header.meta meta.txt
cp ../header.meta meta_CADD15.txt
cp ../header.meta meta_CADD20.txt

bcftools view -G ../original/$1 -Ov -o p1.vcf


sed '/#/d' p1.vcf >p1


##variant info
cut -f 1-2,4-5 p1 > c1
cut -f 1-5 p1 >c1a

#align the order of alt allele as appears in c1

#cut -f 4 c1 >alt
cut -f 8 p1 >c
cut -f 3 -d';' c | sed 's/CSQ\=//g' >csq

paste p1 csq | awk '$5 !~/,/' | cut -f 1-7,9 >p1_s
awk '$5 ~/,/' p1 >p1_m


paste p1 csq | awk '$5 ~/,/' |while read i; do echo $i | cut -f 5 -d' ' |sed 's/\,/\n/g'>j; echo $i | cut -f 9 -d' ' | sed 's/\,/\n/g' >k ; cat j |while read l; do grep -w "$l" k; done |paste -sd',';done >order

paste p1_m order |awk '$9 !~/,/' |cut -f 5 > alt_re
paste p1_m order |awk '$9 !~/,/' |cut -f 1-8 > p1_re
paste p1_m order |awk '$9 !~/,/' |cut -f 8 | cut -f 3 -d';' | cut -f 1 > csq_re
paste p1_m order |awk '$9 ~/,/' |cut -f 1-7,9 > p1_1

##repeat
paste alt_re csq_re |sed 's/CSQ\=//g' |while read i; do echo $i | cut -f 1 -d' ' |sed 's/\,/\n/g' | awk '{if (length($1)==1) print"--"; else print$i}'>j; echo $i | cut -f 2 -d' ' | sed 's/\,/\n/g' >k ; cat j |while read l; do m=${l:1};grep -w "$m" k; done |paste -sd',';done > order_re


paste p1_re order_re |awk '$9 ~/,/' |cut -f 1-7,9 > p1_2

cat p1_s p1_1 p1_2 >p1_order
cat c1a | while read i; do grep -w "$i" p1_order;done >p1_u
#cat p1_s p1_1 p1_2 |\
#        sort -k1,1 -k2,2n |\
#        awk -F"\t" '{print$1"_"$2"_"$4"_"$5,$6,$7,$8}' > p1_order
#join -j 1 c1a p1_order |sed 's/ /\t/g'> p1_u


#awk 'NR==FNR{x++} END{ if(x!=FNR){print"mismatch ERROR on ma ordering"} }' p1 p1_u
##tba: n of alleles? suppose n.alt=10 atm

cut -f 8 p1_u >c_u


##allele funtional consequence
cut -f 2 -d'|' c_u  >c2

##gene with ensemblID; Note: there are 806 x-genes crossing chunks

#cut -f 3-4 -d'|' c_u|sed 's/|/_/g' >c3
#sort -u c3 > gene.lst
bedtools intersect \
        -wao \
        -a p1.vcf \
        -b ~/ref/gencode.v43.basic.p50.bed ~/ref/IBD_gwas_locus_0927.bed |\
        cut -f 1-5,13 >p1.bed
~/bin/datamash-1.8/datamash -g 1,2,3,4,5 collapse 6 <p1.bed |\
    cut -f 6 >c3

#cut -f 3-4 -d'|' c_u|sed 's/|/_/g' >c3
#perl -ne 'print join("\n", split(/\,/,$_));print("\n")' c3 |sort -u |grep -E 'locus|ENSG'>gene.lst
perl -ne 'print join("\n", split(/\,/,$_));print("\n")' c3 |sort -u |grep -E 'locus'>gene.lst

##AF
cut -f 1 -d';' c | sed 's/AF\=//g' >c4a
awk -F, -v OFS=, 'NR==FNR{if(max<10)max=10;next};
                           {NF=10}1' c4a{,} | sed 's/\,/\t/g' >c4


##raw_score_all
cut -f 3 -d';' c_u |awk -F"|" '{OFS="\t"}{print$6,$12,$18,$24,$30,$36,$42,$48,$54,$60}' >c5

##phred_score >=15, which set smaller scores as 0
awk -F"\t" '{OFS=FS}{for(i=1;i<=NF;i++)if($i<1.387112){$i="";}}1' c5 >c5a

##phred_score >=20
awk -F"\t" '{OFS=FS}{for(i=1;i<=NF;i++)if($i<2.097252){$i="";}}1' c5 >c5b

##genotype
zgrep -v '#' ../original/$1 | cut -f 10- | awk -F"\t" '{OFS=FS}{for(i=1;i<=NF;i++) $i=substr($i,1,3)}1' >c6



##merge;
paste c1 c2 c3 c4 c5 c6 >> meta.txt
paste c1 c2 c3 c4 c5a c6 >> meta_CADD15.txt
paste c1 c2 c3 c4 c5b c6 >> meta_CADD20.txt

rm c6
rm c* p1* alt_re order*
rm k j


##CADD_all
while read gene;
    do
        cp /home/gc1a20/runukbb/header.meta ${gene}.meta
        grep -w "$gene" meta.txt >> ${gene}.meta
        bgzip ${gene}.meta
    done < gene.lst

##CADD_15
while read gene;
    do
        cp /home/gc1a20/runukbb/header.meta ${gene}_CADD15.meta
        grep -w "$gene" meta_CADD15.txt >> ${gene}_CADD15.meta
        bgzip ${gene}_CADD15.meta
    done < gene.lst

##CADD_20
while read gene;
    do
        cp /home/gc1a20/runukbb/header.meta ${gene}_CADD20.meta
        grep -w "$gene" meta_CADD20.txt >> ${gene}_CADD20.meta
        bgzip ${gene}_CADD20.meta
    done < gene.lst

#rm meta*.txt

####x-genes
##..
##mkdir -p x-genes
##
##for i in x-gene.lst;
##    do cp header $gene.meta
##        ls -1 ukb23156*/*.meta.gz |grep "$i" |xargs zcat | grep -v '#'  >> x-genes/$gene.meta
##        bgzip x-genes/$gene.meta
##    done
##
##
##ls -1 ukb23156*/*.meta.gz |fgrep -f x-gene.lst - |xargs rm
##

