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
    bgzip -c f3.vcf > f3.vcf.gz
    tabix -p vcf f3.vcf.gz
    
  #####################################################################
   ## picard FilterVcf  -I f3.vcf.gz -O f4.vcf --MIN_DP 8 --MIN_AB 0.15 --MAX_RECORDS_IN_RAM 300000
   #################################
    LINES=3000
    vcf=f3.vcf.gz
    mkdir -p chunks
    bcftools view -h \$vcf > header.txt
    bcftools view -H \$vcf | split -l \$LINES - chunks/chunk_
    
    for chunk in chunks/chunk_*; do
        cat header.txt \$chunk | bgzip -c > \${chunk}.vcf.gz
        rm \$chunk
    done
    rm header.txt
    
    # Parallel processing
    process_chunk() {
        chunk=\$1
        output=\${chunk%.vcf.gz}.filtered.vcf.gz
        picard -Xmx20g FilterVcf \
            -I \$chunk \
            -O \$output \
            --MIN_DP 8 \
            --MIN_AB 0.15 \
            --MAX_RECORDS_IN_RAM 100000
    }
    export -f process_chunk
    find chunks -name '*.vcf.gz' | xargs -n 1 -P 2 -I {} bash -c 'process_chunk "{}"'
    
    # Merging and indexing
    filtered_chunks=\$(ls chunks/*.filtered.vcf.gz)
    bcftools concat -a -Oz -o f4.vcf.gz \$filtered_chunks
    tabix -p vcf f4.vcf.gz
    rm -r chunks
  #########################################
  ##### FilterVcf -I f3.vcf.gz -O f4.vcf -MIN_DP 8 -MIN_AB 0.15
  #####################################################################
  ##  bgzip -@ 10 -c f4.vcf > f4.vcf.gz
  ##   tabix -p vcf f4.vcf.gz
    bcftools view f4.vcf.gz -R ${params.xgen_bed} -Oz -o f5.vcf.gz
    
    rm -r f3.vcf f4.vcf.gz f3.vcf.gz
    
   
    
    """
}

process Pre_processing_2 {
  publishDir "${params.output}", mode: "copy", overwrite: true

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

                      