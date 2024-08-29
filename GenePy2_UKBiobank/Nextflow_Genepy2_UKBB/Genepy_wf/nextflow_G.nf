#!/usr/bin/env nextflow

nextflow.enable.dsl=2
// import java.nio.file.Files
// import java.nio.file.Paths

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


process Genepy_score {
    publishDir "${params.output}/Genepy_score/Genepy2_score/", mode: "copy", overwrite: true
    maxForks 4
    
    input:
    path("Meta_files")

    output:
    file("*.txt")

    script:
    """
    source activate drop
    echo ${params.basedir}
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
         """.stripIndent()
     
 input_files = Channel.fromPath(params.input.list)
        .splitText { it.strip() }
        .map { it -> file(it) }
        .filter { folder ->
            // Only pass non-empty folders to the next step
            def folderContent = folder.listFiles()
            if (folderContent && folderContent.size() > 0) {
                return true
            } else {
                println "Skipping empty folder: ${folder.path}"
                return false
            }
        }

    // Run the Genepy_score process for non-empty folders
    println "Genepy script folder: ${params.genepy_py}  ${params.basedir}"
    Genepy_score(input_files)
}
//workflow.onError {
//  file(params.output).deleteDir()
// }

                      
