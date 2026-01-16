#!/usr/bin/env nextflow

nextflow.enable.dsl=2

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    VALIDATE & PRINT PARAMETER SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Check input path parameters to see if they exist
def checkPathParamList = [ 
    params.vcf
]

include { CADD_score } from "./modules/CADD"  
include { VEP_score } from "./modules/VEP"  
include { Pre_processing_1 } from "./modules/Pre_pr1"  
include { Pre_processing_2 } from "./modules/Pre_pr2"  
include { Pre_processing_3 } from "./modules/Pre_pr3"  
include { GeneReattach } from "./modules/Gene_reattach"
include { Genepy_score } from "./modules/Genepy"




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
       chrx = channel.fromPath("${params.vcf}/*_${params.chromosomes}_*.vcf.gz").map { file -> 
                      def filename = file.baseName  // Extracts filename without the .vcf.gz extension
                      return [chromosomeList,filename,file]       // Returns a tuple with [full path, base filename]
                      }


      CADD_score(chrx)
      VEP_score(CADD_score.out.pre_proc_1)
      Pre_processing_1(VEP_score.out)
      Pre_processing_2(Pre_processing_1.out)
      Pre_processing_3(Pre_processing_2.out)
      
      def meta15 = Pre_processing_3.out.meta_files15.collect().map { genes_list -> ["15",chromosomeList, genes_list] }
      def meta20 = Pre_processing_3.out.meta_files20.collect().map { genes_list -> ["20",chromosomeList, genes_list] }
      def metaALL = Pre_processing_3.out.meta_filesALL.collect().map { genes_list -> ["ALL",chromosomeList, genes_list] }
      x_combo= meta15.concat(meta20).concat(metaALL)
      GeneReattach(x_combo)
      def result = GeneReattach.out.path_.flatten().map{[it]}.map { path ->
            path1 = path.toString()
            def chromosome = (path1 =~ /chr([1-9]|1[0-9]|2[0-4])\b/)[0][0]
            def cadd_score = (path1.contains('metafilesALL')) ? 'ALL' :
                             (path1.contains('metafiles20')) ? '20' :
                             (path1.contains('metafiles15')) ? '15' : 'ALL'
            [path, chromosome, cadd_score]
        }
      result.view()
      Genepy_score(result)

}

workflow.onComplete {
   println ( workflow.success ? """
       Pipeline execution summary
       ---------------------------
       Completed at: ${workflow.complete}
       Duration    : ${workflow.duration}
       Success     : ${workflow.success}
       workDir     : ${workflow.workDir}
       exit status : ${workflow.exitStatus}
       """ : """
       Failed: ${workflow.errorReport}
       exit status : ${workflow.exitStatus}
       """
   )
}
workflow.onError {
  file(params.output).deleteDir()
}

                      
