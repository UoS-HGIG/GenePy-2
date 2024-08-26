####### ./jobscript.sh /vcf_file_path/
SCRIPT=$(readlink -f .)
vcf_directory=$1   ## input vcf path

#chromosomes=$(ls $vcf_directory/*.vcf.gz | grep -oP 'chr[0-9XYM]+' | sort -u)
chromosomes=$(ls $vcf_directory/*.vcf.gz | grep -oP 'chr6' | sort -u)
echo "all chromosomes : $chromosomes"
for chrom in $chromosomes; do
    sbatch $SCRIPT/split_chr.sh $vcf_directory $chrom
    echo "jobscript send for $chrom"
done




