process VEP_score {
   //publishDir "${params.output}/${vcf_n}", mode: "copy", overwrite: true
   //maxForks 10

  input:
  tuple val(chrx), path("p1.vcf"), path("wes.tsv.gz"), path("wes.tsv.gz.tbi"), val(vcf_n) , file(vcfFile)
  
  output:
   tuple path("${chrx}.p1.vep.vcf"), val(vcf_n), file(vcfFile), val(chrx) ,emit: vep_out
  
  script:
  
    """
     vep  -i "p1.vcf" --offline --assembly GRCh38 --vcf --fork 10 --cache --force_overwrite --pick_allele --plugin CADD,${params.plugin1},${params.plugin2},"wes.tsv.gz"  --fields "Allele,Consequence,SYMBOL,Gene,gnomadE_AF,CADD_RAW,gnomadRF_RF_flag" -o "${chrx}.p1.vep.vcf" --dir_cache ${params.homos_vep} --dir ${params.homos_vep} --dir_plugins ${params.vep_plugins}
    """
}