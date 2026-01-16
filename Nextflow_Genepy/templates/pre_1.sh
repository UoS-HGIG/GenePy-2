grep '#' p12.vcf >f6_unsorted.vcf

awk '$1 !~/#/' f5.vcf |awk  -F"\t" '{OFS=FS}{for (i=10;i<=NF;i++) if ($i ~/\.\//) $i="./."; else $i=substr($i,1,3)}1' |  awk -F"\t" '{OFS=FS}{for (i=10;i<=NF;i++) gsub(/\|/,"/",$i)}1' >> f6_unsorted.vcf
LC_ALL=C sort -u f6_unsorted.vcf > f6.vcf
