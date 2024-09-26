#!/usr/bin/env nextflow

nextflow.enable.dsl=2
// import java.nio.file.Files
// import java.nio.file.Paths

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check input path parameters to see if they exist
def checkPathParamList = [ 
    params.vcf
]


process CADD_score {
  publishDir "${params.output}", mode: "copy", overwrite: true
  maxForks 5
  input:
  tuple val(chrx),file(vcfFile)

  output:
  tuple val(chrx), path("p1.vcf"), path("wes_${chrx}.tsv.gz"), path("wes_${chrx}.tsv.gz.tbi"), emit: pre_proc_1
  
  script:
    """
    
    bcftools view -G ${vcfFile} -Ov -o p1.vcf
    
    awk -F"\t" '\$1 ~/#/ || length(\$4)>1||length(\$5)>1' p1.vcf | sed '3383,\$s/chr//g' p1.vcf > ${chrx}.p11.vcf
    CADD.sh -c 8 -o wes_${chrx}.tsv.gz ${chrx}.p11.vcf
    tabix -p vcf wes_${chrx}.tsv.gz
    
    
    """
}

process VEP_score {
   publishDir "${params.output}", mode: "copy", overwrite: true
   maxForks 5

  input:
  tuple val(chrx), path("p1.vcf"), path("wes.tsv.gz"), path("wes.tsv.gz.tbi")
  
  output:
   path("${chrx}.p1.vep.vcf"),emit: vep_out
  
  script:
  
    """
     vep  -i "p1.vcf" --offline --assembly GRCh38 --vcf --fork 10 --cache --force_overwrite --pick_allele --plugin CADD,${params.plugin1},${params.plugin2},"wes.tsv.gz"  --fields "Allele,Consequence,SYMBOL,Gene,gnomadE_AF,CADD_RAW,gnomadRF_RF_flag" -o "${chrx}.p1.vep.vcf" --dir_cache ${params.homos_vep} --dir ${params.homos_vep} --dir_plugins ${params.vep_plugins}
    """
}


process Pre_processing_1 {
  publishDir "${params.output}", mode: "copy", overwrite: true
  maxForks 5

  input:
  path(x)
  output:
  path("f5.vcf.gz")
  
  
  shell:
    """
    gunzip -c ${params.vcf} | grep -v '##'|cut -f 9-> p2
    grep -v '##' ${x} > p1
    grep '##' ${x} > f3.vcf
    paste p1 p2 >> f3.vcf
    rm -r p1 p2
    bcftools view -h  ${params.vcf} | grep '^##FORMAT=' > format.txt
    sed -i '1 r format.txt' f3.vcf
    #####
    bcftools +fill-tags f3.vcf -- -t 'FORMAT/AB:1=float((AD[:1]) / (DP))' | bgzip -c > f3.vcf.gz
    rm f3.vcf
    tabix -p vcf f3.vcf.gz
    bcftools filter -S . --include 'FORMAT/DP>=8 & FORMAT/AB>=0.15 |FORMAT/GT="0/0"' -Oz -o f3b.vcf.gz f3.vcf.gz
    tabix -p vcf f3b.vcf.gz
    bcftools +fill-tags f3b.vcf.gz -- -t 'F_MISSING' | bcftools view -e '(CHROM=="chrY" & INFO/F_MISSING>=0.56)'| bcftools view -i 'INFO/F_MISSING<0.12' -Oz -o f4.vcf.gz
    
    tabix -p vcf f4.vcf.gz
    bcftools view f4.vcf.gz -R ${params.xgen_bed} -Oz -o f5.vcf.gz
    
   
    
    """
}

process Pre_processing_2 {
  publishDir "${params.output}", mode: "copy", overwrite: true
  maxForks 5

  input:
  file("f5.vcf.gz")
  output:
  tuple path("c1"), path("c2"), path("c3"), path("c4"),path("c5"),path("c5a"),path("c5b"),path("meta_CADDALL.txt"),path("meta_CADD15.txt"),path("meta_CADD20.txt"),path("gene.lst"),path("f5.vcf.gz"),path("header.meta")
  
  shell:
    """
    cat ${params.header_meta} > meta_CADD_head
    cat ${params.IBD_gwas.bed} > IBD.bed
    cat ${params.Genecode_p50.bed} > p50.bed
    ## bgzip -c "f5.vcf" > f5.vcf.gz
    bcftools view -h f5.vcf.gz | grep -v "##" | cut -f 10- >p
    ${template("pre_1.sh")}
    """
}

process Pre_processing_3 {
  publishDir "${params.output}", mode: "copy", overwrite: true
  maxForks 5

  input:
  tuple path("c1"), path("c2"), path("c3"), path("c4"),path("c5"),path("c5a"),path("c5b"),path("meta_CADDALL.txt"),path("meta_CADD15.txt"),path("meta_CADD20.txt"),path("gene.lst"),path("f5.vcf.gz"),path("header.meta")
  output:
  // tuple file("meta_CADD${params.cadd_filter}.txt"),file("gene.lst"), emit: gene_list
  // tuple file("meta_CADDALL.txt"),file("meta_CADD15.txt"),file("meta_CADD20.txt")
  tuple file("metafiles/*CADD${params.cadd_filter}.meta"),file("gene.lst"), emit: meta_files
  
  shell:
    """
    echo "Processing 3"
    mkdir -p metafiles
     ${template("pre_2.sh")}
    """
}



// Define workflow
workflow {
    println """\
         G E N E P Y           P I P E L I N E
          ===================================
     G E N O M I C --------------- M E D I C I N E 
                         UoS
                   Sarah Ennis Lab
                     Iman Nazari
          ===================================
         Samples         : ${params.vcf}
         """.stripIndent()
     
  
       // def chromosomeList = params.chromosomes.split(',').collect { it.trim().replaceAll('"', '') }
       chromosomeList = params.chromosomes
       println "Chromosome list: $chromosomeList"
       chr = channel.from(chromosomeList).map { chrx -> tuple(chrx, file(params.vcf)) }.view()
       CADD_score(chr)
       VEP_score(CADD_score.out.pre_proc_1)
//       Merging_chunks(VEP_score.out.vep_out.flatten().collect())
      Pre_processing_1(VEP_score.out.collect().view())
      Pre_processing_2(Pre_processing_1.out.collect().view())
      Pre_processing_3(Pre_processing_2.out.collect().view())
 
}

//workflow.onError {
//  file(params.output).deleteDir()
//}

                      
