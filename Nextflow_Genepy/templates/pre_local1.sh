#cut -f 4 -d' ' p1_u|awk -F"," '{OFS=FS}{for (i=1;i<=NF;i++) if ($i ~/*/) $i="*|*|*|||||||"}1' > c_u
cat p1_u |head
cut -f 4 -d' ' p1_u|awk -F"," '{OFS=FS}{for (i=1;i<=NF;i++) if ($i ~/\*/) $i="*|*|*|||||||"}1' > c_u
##allele funtional consequence
cut -f 2 -d'|' c_u  >c2 #Not used for GenePy

##gene with ensemblID; Note: there are 806 x-genes crossing chunks
cut -f 1-8 f6.vcf >> f61.vcf


bedtools intersect \
    -wao \
    -a f61.vcf \
    -b /mainfs/hgig/public/GENE_DATABASES/gencode.v45.annotation_2.bed |\
    cut -f 1-5,12 >f61.bed

/mainfs/hgig/private/software/datamash-1.8/datamash -g 1,2,3,4,5 collapse 6 <f61.bed |\
    cut -f 6 >c3

perl -ne 'print join("\n", split(/\,/,$_));print("\n")' c3 |sort -u >gene.lst


##AF
cut -f 3 -d';' c_u |awk -F"|" '{OFS="\t"}{if ($5>0) print$6,$15,$24,$33,$42,$51,$60,$69,$78,$87; else print$8,$17,$26,$5,$44,$53,$62,$71,$80,$89}' >c4


##raw_score_all
cut -f 3 -d';' c_u |awk -F"|" '{OFS="\t"}{print$9,$18,$27,$36,$45,$54,$63,$72,$81,$90}' >c5

##phred_score >=15, which set smaller scores as 0
awk -F"\t" '{OFS=FS}{for(i=1;i<=NF;i++)if($i<1.387112){$i="";}}1' c5 >c5a

##phred_score >=20
awk -F"\t" '{OFS=FS}{for(i=1;i<=NF;i++)if($i<2.097252){$i="";}}1' c5 >c5b

#genotype
awk -F"," '{if ($0 !~/\*/) print"0"; else for (i=1;i<=NF;i++) if ($i ~/\*/) print i}' alt >x1
cut -f 10- f6 >x2
paste x1 x2 |while read i; do pos=$(echo $i |cut -f 1 -d' '); echo $i | cut -f 2- -d' '| sed "s/${pos}/0/g" | sed 's/ /\t/g';done >c6


##merge;
paste c1 c2 c3 c4 c5 c6 > CADDALL.txt
paste c1 c2 c3 c4 c5a c6 > CADD15.txt
paste c1 c2 c3 c4 c5b c6 > CADD20.txt
cat header_CADDALL.txt CADDALL.txt > meta_CADDALL.txt
cat header_CADD15.txt CADD15.txt > meta_CADD15.txt
cat header_CADD20.txt CADD20.txt > meta_CADD20.txt
######rm c6