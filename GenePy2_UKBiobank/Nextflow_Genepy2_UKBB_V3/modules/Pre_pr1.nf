process Pre_processing_1 {
  //publishDir "${params.output}/${vcf_n}", mode: "copy", overwrite: true
  //maxForks 10

  input:
  tuple path(x), val(vcf_n), file(vcfFile), val(chrx) 
  
  output:
  tuple path("f5.vcf.gz"), val(vcf_n), val(chrx) 
  
  
  shell:
    """
    gunzip -c ${vcfFile} | grep -v '##'|cut -f 9-> p2
    grep -v '##' ${x} > p1
    grep '##' ${x} > f3.vcf
    paste p1 p2 >> f3.vcf
    rm -r p1 p2
    bcftools view -h  ${vcfFile} | grep '^##FORMAT=' > format.txt
    sed -i '1 r format.txt' f3.vcf
    #####
    bcftools +fill-tags f3.vcf -- -t 'FORMAT/AB:1=float((AD[:1]) / (DP))' | bgzip -c > f3.vcf.gz
    rm f3.vcf
    tabix -p vcf f3.vcf.gz
    bcftools filter -S . --include 'FORMAT/DP>=8 & FORMAT/AB>=0.15 |FORMAT/GT="0/0"' -Oz -o f3b.vcf.gz f3.vcf.gz
    tabix -p vcf f3b.vcf.gz

    cat ${params.ethnicity} > ethnicity.txt
    bcftools +fill-tags f3b.vcf.gz -- -S ethnicity.txt -t 'HWE,F_MISSING' | bcftools view -e '(CHROM=="chrY" & INFO/F_MISSING>=0.56 & INFO/HWE_1>(0.05/15922704))'| bcftools view -i 'INFO/F_MISSING<0.12 & INFO/HWE_1>(0.05/15922704)' -Oz -o f4.vcf.gz

    tabix -p vcf f4.vcf.gz
    
    bcftools view f4.vcf.gz -R ${params.xgen_bed} -Oz -o f5.vcf.gz
    
    
    """
}