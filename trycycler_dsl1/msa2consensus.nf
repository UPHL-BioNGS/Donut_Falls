#!/usr/bin/env nextflow

println("Currently using the Donut Falls workflow for use with nanopore sequencing\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v.0.20220430")
println("")

params.outdir = workflow.launchDir + '/trycycler_3'

if ( Runtime.runtime.availableProcessors() < 6 ) {
  params.maxcpus = Runtime.runtime.availableProcessors()
} else {
  params.maxcpus = 6
}

params.cluster_directory = workflow.launchDir + '/trycycler_2'
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

process transpose {
  tag "${trycycler}"
  cpus 1
  container 'staphb/trycycler:latest'

  input:
  set val(sample), path(trycycler) from clusters

  output:
  tuple sample, path("${trycycler}/cluster*") optional true into transposed

  shell:
  '''
  ls !{trycycler}*
  '''
}

params.trycycler_msa_options = ''
process trycycler_msa {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}_${trycycler}"
  cpus params.maxcpus
  container 'staphb/trycycler:latest'

  input:
  set val(sample), path(trycycler) from transposed.transpose()

  output:
  tuple sample, env(cluster), path(trycycler) into msa
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err
    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version >> $log_file
    cluster=!{trycycler}
    trycycler msa !{params.trycycler_msa_options} \
      --cluster_dir !{trycycler} \
      --threads !{task.cpus} \
      2>> $err_file >> $log_file
  '''
}

params.trycycler_partition_options = ''
process trycycler_partition {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}_${cluster}"
  cpus params.maxcpus
  container 'staphb/trycycler:latest'

  input:
  set val(sample), val(cluster), path(trycycler), file(fastq) from msa.combine(fastq, by: 0)

  output:
  tuple sample, cluster, path(trycycler) into partition
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err
    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version >> $log_file
    trycycler partition !{params.trycycler_partition_options} \
      --reads !{fastq} \
      --cluster_dirs !{trycycler} \
      --threads !{task.cpus} \
      2>> $err_file >> $log_file
  '''
}

params.trycycler_consensus_options = ''
process trycycler_consensus {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}_${cluster}"
  cpus params.maxcpus
  container 'staphb/trycycler:latest'

  input:
  set val(sample), val(cluster), path(trycycler) from partition

  output:
  tuple sample, cluster, path(trycycler) into consensus
  file("${trycycler}/${sample}_${cluster}_trycycler_consensus.fasta") into consensus_all
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err
    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version >> $log_file
    trycycler consensus !{params.trycycler_consensus_options} \
      --cluster_dir !{trycycler} \
      --threads !{task.cpus} \
      2>> $err_file >> $log_file
    cp !{trycycler}/7_final_consensus.fasta !{trycycler}/!{sample}_!{cluster}_trycycler_consensus.fasta
  '''
}

process trycycler_combine_consensus {
  publishDir "${params.outdir}", mode: 'copy'
  tag "combine all"
  cpus 1
  container 'staphb/trycycler:latest'

  input:
  file(consensus) from consensus_all.collect()

  output:
  file("${task.process}/*_trycycler_consensus.fasta")

  shell:
  '''
    mkdir !{task.process}
    for sample in $(ls *fasta | sed 's/_cluster.*//g' | sort | uniq)
    do
      cat $sample*fasta > !{task.process}/${sample}_trycycler_consensus.fasta
    done
  '''
}

params.medaka_options = ''
process medaka {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}_${cluster}"
  cpus params.medcpus
  container 'ontresearch/medaka:latest'

  input:
  set val(sample), val(cluster), path(trycycler) from consensus

  output:
  path(trycycler)
  file("${trycycler}/${sample}_${cluster}_medaka_consensus.fasta") into medaka_fastas_all
  tuple sample, cluster, file("${trycycler}/7_final_consensus.fasta") into medaka_fastas, medaka_fastas2
  tuple env(sample_cluster), file("${trycycler}/7_final_consensus.fasta") into medaka_fastas_pilon
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err
    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    medaka --version >> $log_file
    sample_cluster="!{sample}_!{cluster}"
    medaka_consensus !{params.medaka_options} \
      -i !{trycycler}/4_reads.fastq \
      -d !{trycycler}/7_final_consensus.fasta \
      -o !{trycycler}/medaka \
      -t !{task.cpus}
    cp !{trycycler}/7_final_consensus.fasta !{trycycler}/!{sample}_!{cluster}_medaka_consensus.fasta
  '''
}

process medaka_combine_consensus {
  publishDir "${params.outdir}", mode: 'copy'
  tag "combine all"
  cpus 1
  container 'ontresearch/medaka:latest'

  input:
  file(consensus) from medaka_fastas_all.collect()

  output:
  file("${task.process}/*_medaka_consensus.fasta")

  shell:
  '''
    mkdir !{task.process}
    for sample in $(ls *fasta | sed 's/_cluster_.*//g' | sort | uniq)
    do
      cat $sample*fasta > !{task.process}/${sample}_medaka_consensus.fasta
    done
  '''
}
