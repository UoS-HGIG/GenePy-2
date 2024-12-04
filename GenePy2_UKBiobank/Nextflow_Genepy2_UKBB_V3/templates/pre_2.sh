#zcat f5.vcf.gz | awk 'BEGIN { OFS="\t" } !/^#/ { for (i = 10; i <= NF; i++) { $i = substr($i, 1, 3) } print }' > c6
#awk 'BEGIN { OFS="\t" } !/^#/ { for (i = 10; i <= NF; i++) { $i = substr($i, 1, 3) } print }' f5.vcf > c6
> meta_CADDALL.txt
> meta_CADD15.txt
> meta_CADD20.txt
region=$1
echo "$region"
cp header.meta meta_CADDALL.txt
cp header.meta meta_CADD15.txt
cp header.meta meta_CADD20.txt

zgrep -v '#' f5.vcf.gz | cut -f 10- | awk -F"\t" '{OFS=FS}{for(i=1;i<=NF;i++) $i=substr($i,1,3)}1' >c6


##merge;
paste c1 c2 c3 c4 c5 c6 >> meta_CADDALL.txt
paste c1 c2 c3 c4 c5a c6 >> meta_CADD15.txt
paste c1 c2 c3 c4 c5b c6 >> meta_CADD20.txt




##CADD_all
while read gene;
    do
        cp header.meta "metafilesALL_${region}/${gene}_CADDALL.meta"
        grep -a -w "$gene" meta_CADDALL.txt >> "metafilesALL_${region}/${gene}_CADDALL.meta"
        ###bgzip ${gene}_CADDALL.meta
    done < gene.lst

##CADD_15
while read gene;
    do
        cp header.meta "metafiles15_${region}/${gene}_CADD15.meta"
        grep -a -w "$gene" meta_CADD15.txt >> "metafiles15_${region}/${gene}_CADD15.meta"
        ###bgzip ${gene}_CADD15.meta
    done < gene.lst

##CADD_20
while read gene;
    do
        cp header.meta "metafiles20_${region}/${gene}_CADD20.meta"
        grep -a -w "$gene" meta_CADD20.txt >> "metafiles20_${region}/${gene}_CADD20.meta"
        ###bgzip ${gene}_CADD20.meta
    done < gene.lst
