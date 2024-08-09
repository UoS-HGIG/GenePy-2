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
  tuple val(chrx), path("${chrx}.p1.vcf"), path("wes_${chrx}.tsv.gz"), path("wes_${chrx}.tsv.gz.tbi"),path("${chrx}.p11.vcf"), emit: pre_proc_1
  
  script:
    """
    
    bcftools view -G ${vcfFile} -Oz -o p1.vcf.gz
    tabix -p vcf p1.vcf.gz
    bcftools view -O z -o ${chrx}.p1.vcf p1.vcf.gz "chr${chrx}"
    awk -F"\t" '\$1 ~/#/ || length(\$4)>1||length(\$5)>1' ${chrx}.p1.vcf | sed 's/^chr//g' > ${chrx}.p11.vcf
    CADD.sh -c 8 -o wes_${chrx}.tsv.gz ${chrx}.p11.vcf
    tabix -p vcf wes_${chrx}.tsv.gz
    bgzip -c ${chrx}.p1.vcf > ${chrx}.p1.vcf.gz
    tabix -p vcf ${chrx}.p1.vcf.gz
    
    """
}

process VEP_score {
   publishDir "${params.output}", mode: "copy", overwrite: true

  input:
  tuple val(chrx), path("p1.vcf"), path("wes.tsv.gz"), path("wes.tsv.gz.tbi"),path("p11.vcf")
  
  output:
   path("${chrx}.p1.vep.vcf"),emit: vep_out
  
  script:
  
    """
    vep  -i "p1.vcf" --offline --assembly GRCh38 --vcf --fork 10 --cache --force_overwrite --pick_allele --plugin CADD,"wes.tsv.gz",${params.plugin1},${params.plugin2} --af_gnomade --af_gnomadg --fields "Allele,Consequence,SYMBOL,Gene,gnomADg_AF,gnomADg_NFE_AF,gnomADe_AF,gnomADe_NFE_AF,CADD_RAW,gnomadRF_RF_flag" -o "${chrx}.p1.vep.vcf" --dir_cache ${params.homos_vep} --dir_plugins ${params.vep_plugins} --dir ${params.homos_vep}
    """
}


process Merging_chunks {
  publishDir "${params.output}", mode: "copy", overwrite: true

  input:
  path(x)
  output:
  path("p1.vep.sorted.vcf")
  
  script:
    """
    ls $x > vep_vcf_merge.list
    gatk MergeVcfs -I vep_vcf_merge.list -O p1.vep.vcf
    gatk SortVcf -I p1.vep.vcf -O p1.vep.sorted.vcf
    
    """
}


process Pre_processing_1 {
   publishDir "${params.output}", mode: "copy", overwrite: true
 
  input:
  path(x)
  output:
  tuple file("f3.vcf.gz"),file("f4.vcf.gz"),file("f4.vcf.gz.tbi"),file("p12.vcf"),file("f5.vcf")
  
  
  shell:
    """
    gunzip -c ${params.vcf} | awk '\$1 !~/##|\\_/' |cut -f 9-> p2
    grep -v '##' ${x} > p1
    grep '##' ${x} > p12.vcf
    paste p1 p2 >> p12.vcf 
    rm -r p1 p2
    bcftools view -h  ${params.vcf} | grep '^##FORMAT=' > format.txt
    sed -i '1 r format.txt' p12.vcf 
    awk -F"\t" '\$7 ~/PASS/ || \$1 ~/#/' p12.vcf > f1.vcf
    bcftools filter -S . --include 'FORMAT/DP>=8 & FORMAT/GQ>=20' -o f2.vcf f1.vcf
    
    bcftools query -f '[%GQ\t]' f2.vcf >temp0
    awk '{a=0; for(i = 3; i<= NF; i++) {a=a+\$i}; {print a/NF}}' temp0 > temp_meanGQ.out
    grep "^#" f2.vcf >temp_header
    grep -v "^#" f2.vcf >temp_body
    sed -i "1s|0|99.0|" temp_meanGQ.out
    paste temp_meanGQ.out temp_body > temp1
    awk '{ if(\$1>35) {print \$0}}' temp1 >temp2
    cut -f 2- temp2 >temp3
    rm -r temp0 temp1 temp2 temp_body temp_meanGQ.out
    
    cat temp_header temp3 > f3.vcf
    rm -r temp3
    
    bgzip -c f3.vcf > f3.vcf.gz
    tabix -p vcf f3.vcf.gz
    
    
    bcftools +fill-tags f3.vcf -- -t 'F_MISSING' | bcftools view -e '(CHROM=="chrY" & INFO/F_MISSING>=0.56)'| bcftools view -i 'INFO/F_MISSING<0.12' -Oz -o f4.vcf.gz
    tabix -p vcf f4.vcf.gz
    bcftools view f4.vcf.gz -R ${params.bed_int} -Ov -o f5.vcf
    
    rm -r f3.vcf f2.vcf f1.vcf
    """
}

process Pre_processing_2 {
   publishDir "${params.output}", mode: "copy", overwrite: true

  input:
  tuple file("f3.vcf.gz"),file("f4.vcf.gz"),file("f4.vcf.gz.tbi"),file("p12.vcf"),file("f5.vcf")
  output:
  file("f6.vcf")
  
  
  shell:
    """
    ${template("pre_1.sh")}
    """
}

process Pre_processing_3 {
  
  publishDir "${params.output}", mode: "copy", overwrite: true

  input:
  tuple file("f3.vcf.gz"),file("f4.vcf.gz"),file("f4.vcf.gz.tbi"),file("p12.vcf"),file("f5.vcf"),file("f6.vcf")
  
  output:
  tuple file("meta_CADD${params.cadd_filter}.txt"),file("gene.lst"), emit: gene_list
  file("meta_CADD*.txt")
  shell:
    """
  bcftools view -h "f5.vcf" | grep -v "##" | cut -f 10- >p
  bcftools view -G ${params.vcf} | bcftools view -h > f61.vcf
  cat ${params.header_meta} > meta_CADD_head.txt
  cat ${params.gene_code_bed} > genecode.bed
  ${template("pre_local.sh")}
    """
}


process SplitByGene {
  
   publishDir "${params.output}/chunks", mode: "copy", overwrite: true

  input:
  tuple file("meta_CADD*.txt"),file("gene.lst")
  
  output:
  tuple file("chunk_*.lst"),env("tot_chunks"), emit: out1
  shell:
    """
    mapfile -t genes < <(grep '^ENSG' "gene.lst")
    
    total_genes=\${#genes[@]}
    chunk_size=${params.chunk_size}
    num_chunks=\$((total_genes / chunk_size))
    remainder=\$((total_genes % chunk_size))
    
    for ((i = 0; i < num_chunks; i++)); do
        start=\$((i * chunk_size))
        end=\$((start + chunk_size))
        chunk_file="chunk_\$((i + 1)).lst"
        printf "%s\\n" "\${genes[@]:start:chunk_size}" > "\$chunk_file"
    done
    
    # Handle the remainder by adding it to the last chunk
    if (( remainder > 0 )); then
        chunk_file="chunk_\$((num_chunks + 1)).lst"
        start=\$((num_chunks * chunk_size))
        printf "%s\\n" "\${genes[@]:start:remainder}" >> "\$chunk_file"
       tot_chunks=\$((num_chunks+1))
    else
    tot_chunks=\$((num_chunks))
    fi
    echo \$tot_chunks
    """
}

process Meta_file_extraction {
  
    // publishDir "${params.output}/Genepy_score/Meta_files", mode: "copy", overwrite: true

  input:
  tuple path("chunk"),path(metaCADD)
  
  output:
  path("Meta_files")
   
  
  script:
    """
    mkdir -p Meta_files
    cat ${metaCADD} | head -n 1 > header.meta
   
    cat "chunk" | while read i; do
    cat header.meta > "Meta_files/\${i}_.meta"
    grep "\$i" ${metaCADD}|awk -F"\t" '{OFS=FS}{for (s=7;s<=16;s++) if(length(\$s)<1 || \$s==0) \$s="3.98e-6"}1' >> "Meta_files/\${i}_.meta"

  done
    """
}


process Genepy_score {
    publishDir "${params.output}/Genepy_score/Genepy2_score/", mode: "copy", overwrite: true
    input:
    path("Meta_files")

    output:
    file("*.txt")

    script:
    """
    source activate drop
    python "${params.genepy_py}" "Meta_files" ${params.cadd_filter}
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
     
        def chromosomeList = params.chromosomes.split(',').collect { it.trim().replaceAll('"', '') }
        println "Chromosome list: $chromosomeList"
        chr = channel.from(chromosomeList).map { chrx -> tuple(chrx, file(params.vcf)) }.view()
        CADD_score(chr)
        VEP_score(CADD_score.out.pre_proc_1)
        Merging_chunks(VEP_score.out.vep_out.flatten().collect())
        Pre_processing_1(Merging_chunks.out.collect().view())
        Pre_processing_2(Pre_processing_1.out.collect().view())
        Pre_processing_3(Pre_processing_1.out.mix(Pre_processing_2.out).collect().view())
        SplitByGene(Pre_processing_3.out.gene_list.collect())
        tot_chunks1= SplitByGene.out.out1.mix(Pre_processing_3.out.gene_list).flatten().filter { it =~ /chunk_.*\.lst/ }.view()
        tot_meta= SplitByGene.out.out1.mix(Pre_processing_3.out.gene_list).flatten().filter { it =~ /meta_CADD/ }.view()
        chr = tot_chunks1
        .combine(tot_meta)
        .filter { chunk, meta -> chunk.name.startsWith("chunk_") && meta.name.contains("meta_CADD") }
        .map { chunk, meta -> tuple(chunk, meta) }
        .view()
        Meta_file_extraction(chr)
        Genepy_score(Meta_file_extraction.out)
}

//workflow.onError {
//  file(params.output).deleteDir()
//}

                      
