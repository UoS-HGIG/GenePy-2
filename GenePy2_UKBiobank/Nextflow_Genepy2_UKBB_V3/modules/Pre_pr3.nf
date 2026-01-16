process Pre_processing_3 {
  //publishDir "${params.output}/${vcf_n}/metafiles", mode: "copy", overwrite: true
  //maxForks 10

  input:
  tuple path("c1"), path("c2"), path("c3"), path("c4"),path("c5"),path("c5a"),path("c5b"),path("gene.lst"),path("f5.vcf.gz"),path("header.meta"), val(vcf_n) , val(chrx) 
  output:
  path("metafiles15_*"), emit: meta_files15
  path("metafiles20_*"), emit: meta_files20
  path("metafilesALL_*"), emit: meta_filesALL
  
  
  shell:
    """
    echo "Processing 3"
    echo "vcf_n: ${vcf_n}"
    region=\$(echo ${vcf_n} | awk -F'[_|.]' '{print \$3"_"\$4}')

    mkdir -p metafiles15_\${region}
    mkdir -p metafiles20_\${region}
    mkdir -p metafilesALL_\${region}

    ${template("pre_2.sh")} \${region}
    """
}