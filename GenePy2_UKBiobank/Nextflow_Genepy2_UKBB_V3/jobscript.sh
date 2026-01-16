####### ./jobscript.sh /iridisfs/hgig/private/in1f24/singularity_genepy/Nextflow/Ukbb200k/vcf_files/ukb23156_chr6_160136439_166464695_v1
SCRIPT=$(readlink -f .)
vcf_directory=$1   ## input vcf path

#chromosomes=$(ls $vcf_directory/*.vcf.gz | grep -oP 'chr[0-9M]+' | sort -u)
#chromosomes=$(ls $vcf_directory/*.vcf.gz | grep -oP 'chr2(?![0-9])|chr21' | sort -u)
#chromosomes=$(ls $vcf_directory/*.vcf.gz | grep -oP 'chr7|chr21' | sort -u)
chromosomes=$(ls $vcf_directory/*.vcf.gz | grep -oP $2 | sort -u)
echo "all chromosomes : $chromosomes"
for chrom in $chromosomes; do
    sbatch $SCRIPT/split_chr.sh $vcf_directory $chrom
    echo "jobscript send for $chrom"
done

