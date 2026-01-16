process Genepy_score {
    publishDir "${params.output}/Genepy_score/${chr}_${cadd}", mode: "copy", overwrite: true
    // maxForks 30
    
    input:
    tuple val(path),val(chr),val(cadd)
     
    output:
    path("*.txt"),optional: true

    script:
    """
    ##echo "Path: ${path[0]}, Chromosome: ${chr}, CADD Score: ${cadd}"
    python -u "${params.genepy_py}" "${path[0]}"
    """
}