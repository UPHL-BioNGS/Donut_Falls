#!/usr/bin/env nextflow

println("Currently using the Donut Falls workflow for use with nanopore sequencing\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v.0.20220430")
println("")

params.outdir = workflow.launchDir + '/trycycler_2'

if ( Runtime.runtime.availableProcessors() < 6 ) {
  params.maxcpus = Runtime.runtime.availableProcessors()
} else {
  params.maxcpus = 6
}

params.cluster_directory = workflow.launchDir + '/after_1'
Channel
  .fromPath("${params.cluster_directory}/*", type:'dir')
  .map { it -> tuple(it.getName(), it ) }
  .view { "Samples with clusters to reconcile : ${it[0]}" }
  .set { clusters }

params.reads = workflow.launchDir + '/trycycler_1/filtlong'
Channel
  .fromPath("${params.reads}/*.fastq", type:'file')
  .ifEmpty {
    println("Could not find filtered fastq files from subsample2cluster. Set with 'params.reads'")
    exit 1
  }
  .map { reads -> tuple(reads.simpleName.replaceAll(~/_filtered/,""), reads ) }
  .view { "Fastq file found : ${it[0]}" }
  .set { fastq }


// params.trycycler_dotplot_options = ''
// process trycycler_dotplot {
//   publishDir "${params.outdir}", mode: 'copy', pattern: '*.{log,err}'
//   publishDir "${params.outdir}/phase2", mode: 'copy', pattern: '*dotplots.png'
//   tag "${cluster}"
//   echo false
//   cpus 1
//   // errorStrategy 'finish'
//   // container 'staphb/trycycler:latest'
//
//   input:
//   set val(sample), path(cluster) from inital_clusters
//
//   output:
//   path("${cluster}/*/dotplots.png")
//   file("logs/${task.process}/${cluster}.${workflow.sessionId}.{log,err}")
//
//   shell:
//   '''
//     mkdir -p logs/!{task.process}
//     log_file=logs/!{task.process}/!{cluster}.!{workflow.sessionId}.log
//     err_file=logs/!{task.process}/!{cluster}.!{workflow.sessionId}.err
//     # time stamp + capturing tool versions
//     date | tee -a $log_file $err_file > /dev/null
//     trycycler --version >> $log_file
//     for cluster in $(ls !{cluster}/cluster* -d )
//     do
//       name=$(echo $cluster | cut -f 2 -d "/")
//       trycycler dotplot !{params.trycycler_dotplot_options} \
//         --cluster_dir $cluster \
//         2>> $err_file >> $log_file
//     done
//   '''
// }

params.trycycler_reconcile_options = ''
params.trycycler_reconcile_minimum = '2'
process trycycler_reconcile {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.maxcpus
  //errorStrategy 'finish'
  container 'staphb/trycycler:latest'

  input:
  set val(sample), path(cluster), file(fastq) from clusters.join(fastq, by: 0).view()

  output:
  path("${cluster}")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version >> $log_file

    ls -d !{cluster}/cluster*/
    for cluster in ls -d !{cluster}/cluster*/
    do
      num_fasta=$(ls $cluster/1_contigs/*.fasta | wc -l)
      if [ "$num_fasta" -ge "!{params.trycycler_reconcile_minimum}" ]
      then
        echo "Working on cluster $cluster" | tee -a $err_file $log_file
        trycycler reconcile !{params.trycycler_reconcile_options} \
          --reads !{fastq} \
          --cluster_dir $cluster \
          --threads !{task.cpus} \
          2>> $err_file >> $log_file
      else
        echo "cluster $cluster only had $num_fasta fastas" >> $log_file
      fi
    done
  '''
}
