// process combine_and_rename {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "${barcode}"
//   echo false
//   cpus 1
//   container 'staphb/filtlong:latest'
//
//   input:
//   set path(barcode), file(sample_key) from fastq_pass_directory
//
//   output:
//   tuple env(sample), file("${task.process}/*.fastq.gz") into combined_fastq
//   file("logs/${task.process}/*.${workflow.sessionId}.{log,err}")
//
//   shell:
//   '''
//     mkdir -p !{task.process} logs/!{task.process}
//     log_file=logs/!{task.process}/!{barcode}.!{workflow.sessionId}.log
//     err_file=logs/!{task.process}/!{barcode}.!{workflow.sessionId}.err
//     # time stamp + capturing tool versions
//     date | tee -a $log_file $err_file > /dev/null
//     if [ "!{sample_key}" != "null" ]
//     then
//       echo "The key for the barcode was found!" >> $log_file
//       sample=$(grep !{barcode} !{sample_key} | cut -f 2 -d "," | head -n 1)
//     fi
//     if [ -z "$sample" ]; then sample=!{barcode} ; fi
//     echo "The sample is $sample for barcode !{barcode}" >> $log_file
//     if $(ls !{barcode}/*fastq 1> /dev/null 2>&1)
//     then
//       cat !{barcode}/*fastq | gzip > !{task.process}/$sample.fastq.gz
//     else
//       cat !{barcode}/*fastq.gz > !{task.process}/$sample.fastq.gz
//     fi
//   '''
// }
// process rename {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "${barcode}"
//   echo false
//   cpus 1
//   container 'staphb/filtlong:latest'
//
//   input:
//   tuple file(fastq), file(sample_key) from fastq
//
//   output:
//   tuple env(sample), file("${task.process}/*.fastq.gz") into renamed_fastq
//   file("logs/${task.process}/*.${workflow.sessionId}.{log,err}")
//
//   shell:
//   '''
//     mkdir -p !{task.process} logs/!{task.process}
//     log_file=logs/!{task.process}/!{barcode}.!{workflow.sessionId}.log
//     err_file=logs/!{task.process}/!{barcode}.!{workflow.sessionId}.err
//
//     # time stamp + capturing tool versions
//     date | tee -a $log_file $err_file > /dev/null
//
//     if [ "!{sample_key}" != "null" ]
//     then
//       echo "The key for the fastq file was found!" >> $log_file
//       sample=$(grep !{fastq} !{sample_key} | cut -f 2 -d "," | head -n 1)
//     fi
//     if [ -z "$sample" ]; then sample=$(echo !{fastq} | sed 's/.fastq.*//g' ); fi
//
//     echo "The sample is $sample for fastq !{fastq}" >> $log_file
//
//     if $(ls *fastq 1> /dev/null 2>&1)
//     then
//       cat *fastq | gzip > !{task.process}/$sample.fastq.gz
//     else
//       cat *fastq.gz > !{task.process}/$sample.fastq.gz
//     fi
//   '''
// }

params.fastq_pass = true
if ( params.fastq_pass ) {
  params.fastq_pass_directory = workflow.launchDir + '/fastq_pass'
  println("Fastq pass directory from nanopore sequencing run : " + params.fastq_pass_directory)

  if ( params.sample_key.exists() ) {
    Channel
      .fromPath("${params.fastq_pass_directory}/barcode*", type:'dir')
      .ifEmpty {
        println("Could not find 'fastq_pass' directories. Set with 'params.fastq_pass_directory'")
        exit 1
      }
      .view { "Directories identified : $it" }
      .combine(sample_key)
      .set { fastq_pass_directory }
    } else {
    Channel
      .fromPath("${params.fastq_pass_directory}/barcode*", type:'dir')
      .ifEmpty {
        println("Could not find 'fastq_pass' directories. Set with 'params.fastq_pass_directory'")
        exit 1
      }
      .map{ it -> tuple(it, null )}
      .set { fastq_pass_directory }
    }
  } else {
    fastq_pass_directory = Channel.empty()
}

// params.fastq_concat = false
// if ( params.fastq_concat) {
//   params.fastq = workflow.launchDir + '/fastq'
//   println("Directory for concatenated fastq files:" + params.fastq)
//   if ( params.sample_key.exists() ) {
//     Channel
//       .fromPath("${params.fastq}/*.fastq", "${params.fastq}/*.fastq.gz")
//       .ifEmpty {
//         println("Could not find 'fastq' files. Set with 'params.fastq'")
//         exit 1
//       }
//       .combine(sample_key)
//       .set { fastq }
//     } else {
//     Channel
//       .fromPath("${params.fastq}/*.fastq", "${params.fastq}/*.fastq.gz")
//       .ifEmpty {
//         println("Could not find 'fastq' files. Set with 'params.fastq'")
//         exit 1
//       }
//       .map{ it -> tuple(it, null )}
//       .set { fastq }
//     }
//   } else {
//     fastq = Channel.empty()
// }

params.sample_key = workflow.launchDir + '/sample_key.csv'
if (params.sample_key.exists()) {
  Channel
    .fromPath(params.sample_key, type:'file')
    .view { "File for linking files to identifiers : $it" }
    .into{sample_key ; sample_key_illumina}
} else {
  sample_key = Channel.empty()
  sample_key_illumina = Channel.empty()
}
