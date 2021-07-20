#!/usr/bin/env nextflow

println("Currently using the Donut Falls workflow for use with nanopore sequencing\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v.20210730")
println("")

params.outdir = workflow.launchDir + '/donut_falls'

params.maxcpus = Runtime.runtime.availableProcessors()
println("The maximum number of CPUS used in this workflow is ${params.maxcpus}")
if ( params.maxcpus < 5 ) {
  params.medcpus = params.maxcpus
} else {
  params.medcpus = 5
}

params.sample_key = workflow.launchDir + '/sample_key.csv'
if (params.sample_key.exists()) {
  Channel
    .fromPath(params.sample_key, type:'file')
    .set{sample_key}
} else {
  sample_key = Channel.empty()
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

params.fastq_pass = true
if ( params.fastq_pass ) {
  params.fastq_pass_directory = workflow.launchDir + '/fastq_pass'
  println("Directory for fastq pass files:" + params.fastq_pass_directory)

  if ( params.sample_key.exists() ) {
    Channel
      .fromPath("${params.fastq_pass_directory}/barcode*", type:'dir')
      .ifEmpty {
        println("Could not find 'fastq_pass' directories. Set with 'params.fastq_pass_directory'")
        exit 1
      }
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

params.phase2 = false
params.cluster_directory = workflow.launchDir + '/donut_falls/trycycler_cluster'
if ( params.phase2 ) {
  println("Beginning second phase of workflow.")
  Channel
    .fromPath("${params.cluster_directory}/*", type:'dir')
    .view { "Clusters to reconcile : $it" }
    .into { clusters_dotplot ; clusters_reconcile }
} else {
  Channel
    .empty()
    .into { clusters_dotplot ; clusters_reconcile }
}

params.illumina = false
params.illumina_fastq = workflow.launchDir + '/illumina_fastq'
if ( params.phase2 && params.illumina ) {
  Channel
    .fromFilePairs(["${params.illumina_fastq}/*_R{1,2}*.fastq.gz",
                    "${params.illumina_fastq}/*_{1,2}.fastq*"], size: 2 )
    .ifEmpty {
      println("Could not find paired-end Illumina reads. Set with directory with 'params.illumina_fastq'")
      exit 1
    }
    .combine(sample_key)
    .set { illumina_fastq }
} else {
  illumina_fastq = Channel.empty()
}

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

process combine_and_rename {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${barcode}"
  echo false
  cpus 1
  container 'staphb/filtlong:latest'

  input:
  tuple path(barcode), file(sample_key) from fastq_pass_directory

  output:
  tuple env(sample), file("${task.process}/*.fastq.gz") into combined_fastq
  file("logs/${task.process}/*.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{barcode}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{barcode}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null

    if [ "!{sample_key}" != "null" ]
    then
      echo "The key for the barcode was found!" >> $log_file
      sample=$(grep !{barcode} !{sample_key} | cut -f 2 -d "," | head -n 1)
    fi
    if [ -z "$sample" ]; then sample=!{barcode} ; fi

    echo "The sample is $sample for barcode !{barcode}" >> $log_file

    if $(ls !{barcode}/*fastq 1> /dev/null 2>&1)
    then
      cat !{barcode}/*fastq | gzip > !{task.process}/$sample.fastq.gz
    else
      cat !{barcode}/*fastq.gz > !{task.process}/$sample.fastq.gz
    fi
  '''
}

// concat_and_rename_fastq
//   .concat(renamed_fastq)
//   .set { combined_fastq }

params.filtlong_options = ""
params.filtlong_min_length = "1000"
params.filtlong_keep_percent = "95"
process filtlong {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1
  container 'staphb/filtlong:latest'

  input:
  tuple val(sample), file(fastq) from combined_fastq

  output:
  tuple sample, file("${task.process}/${sample}_filtered.fastq") into filtered_fastq_subsample, filtered_fastq_cluster, filtered_fastq_reconcile
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    filtlong --version >> $log_file

    filtlong !{params.filtlong_options} \
      --min_length !{params.filtlong_min_length} \
      --keep_percent !{params.filtlong_keep_percent} \
      !{fastq} 2>> $err_file \
      > !{task.process}/!{sample}_filtered.fastq
  '''
}

params.trycycler_subsample_options = ""
params.trycycler_subsample_count = "12"
if ( params.trycycler_subsample_count.toInteger() % 3 > 0 ) {
  println("For this workflow, 'params.trycycler_subsample_count' (which you set as ${params.trycycler_subsample_count}) needs to be divisible by three. Sorry!")
  exit 1
}
subsample_count = params.trycycler_subsample_count.toInteger() / 3
twice_count = subsample_count * 2
process trycycler_subsample {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.maxcpus
  //errorStrategy 'ignore' // because this will attempt to subsample even when there aren't enough reads

  input:
  tuple val(sample), file(fastq) from filtered_fastq_subsample

  output:
  tuple sample, env(genome_size), file("${task.process}/${sample}/canu/sample*.fastq") into subsampled_canu_fastq
  tuple sample, file("${task.process}/${sample}/flye/sample*fastq") into subsampled_flye_fastq
  tuple sample, file("${task.process}/${sample}/raven/sample*fastq") into subsampled_raven_fastq
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/!{sample}/{canu,flye,raven} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version 2>> $err_file >> $log_file

    trycycler subsample !{params.trycycler_subsample_options} \
      --reads !{fastq} \
      --count !{params.trycycler_subsample_count} \
      --threads !{task.cpus} \
      --out_dir !{task.process}/!{sample} \
      2>> $err_file >> $log_file

    genome_size=$(grep "Estimated genome size:" $err_file | awk '{print $4}' | sed 's/,//g')

    ls !{task.process}/!{sample}/*fastq  | head -n !{subsample_count} | xargs -I {} cp {} !{task.process}/!{sample}/canu/.
    ls !{task.process}/!{sample}/*fastq  | tail -n !{subsample_count} | xargs -I {} cp {} !{task.process}/!{sample}/flye/.
    ls !{task.process}/!{sample}/*fastq  | tail -n !{twice_count} | head -n !{subsample_count} | xargs -I {} cp {} !{task.process}/!{sample}/raven/.
  '''
}

params.flye_options = ''
process flye {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.maxcpus
  container 'staphb/flye:latest'

  input:
  tuple val(sample), file(fastq) from subsampled_flye_fastq

  output:
  tuple sample, file("${task.process}/${sample}_flye.fa") into flye_fastas
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    flye --version 2>> $err_file >> $log_file

    for fastq in !{fastq}
    do
      name=$(echo $fastq | sed 's/.fastq//g')
      flye !{params.flye_options} \
        --nano-raw $fastq \
        --threads !{task.cpus} \
        --plasmids \
        --out-dir !{task.process}/!{sample} \
        2>> $err_file >> $log_file

        mv !{task.process}/!{sample}/assembly.fasta !{task.process}/!{sample}/${name}_!{task.process}.fasta
    done

    exit 1
  '''
}

// params.miniasm_and_minipolish_options = ''
// process miniasm_and_minipolish {
//   publishDir "${params.outdir}", mode: 'copy'
//   tag "${sample}"
//   echo false
//   cpus 4
//
//   input:
//   tuple val(sample), file(fastq) from subsampled_miniasm_fastq
//
//   output:
//   tuple sample, file("${task.process}/${sample}_flye.fa") into miniasm_fastas
//   file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")
//
//   shell:
//   '''
//     mkdir -p !{task.process} logs/!{task.process}
//     log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
//     err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err
//
//     # time stamp + capturing tool versions
//     date | tee -a $log_file $err_file > /dev/null
//     miniasm_and_minipolish.sh --version 2>> $err_file >> $log_file
//
//     miniasm_and_minipolish.sh \
//       !{fastq} \
//       !{task.cpus} \
//       > assembly_02.gfa
//
//     any2fasta assembly_02.gfa > assemblies/assembly_02.fasta
//     rm assembly_02.gfa
//   '''
// }

//miniasm_and_minipolish.sh read_subsets/sample_02.fastq "$threads" > assembly_02.gfa && any2fasta assembly_02.gfa > assemblies/assembly_02.fasta && rm assembly_02.gfa

params.raven_options = ''
params.raven_polishing_rounds = '2'
process raven {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.medcpus
  //container 'docker://quay.io/biocontainers/raven-assembler:1.5.1'

  input:
  tuple val(sample), file(fastq) from subsampled_raven_fastq

  output:
  tuple sample, file("${task.process}/${sample}/*.fasta") into raven_fastas
  file("${task.process}/${sample}/*.gfa")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    raven --version 2>> $err_file >> $log_file

    for fastq in !{fastq}
    do
      name=$(echo $fastq | sed 's/.fastq//g')
      raven !{params.raven_options} \
        --polishing-rounds !{params.raven_polishing_rounds} \
        --threads !{task.cpus} \
        --graphical-fragment-assembly !{task.process}/!{sample}/$name.gfa \
        $fastq \
        2>> $err_file \
        > !{task.process}/!{sample}/${name}_!{task.process}.fasta
    done
  '''
}

params.canu_options = ''
process canu {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.maxcpus
  container 'staphb/canu:latest'

  input:
  tuple val(sample), val(genome_size), file(fastq) from subsampled_canu_fastq

  output:
  tuple sample, file("${task.process}/${sample}_canu.fa") into canu_fastas
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    canu --version 2>> $err_file >> $log_file

    echo "Genome size used : " !{genome_size} >> $log_file

    for fastq in !{fastq}
    do
      name=$(echo $fastq | sed 's/.fastq//g')
      canu !{params.canu_options} \
       -p ${name}_!{task.process} \
       -d !{task.process}/!{sample}  \
       genomeSize=!{genome_size} \
       -nanopore $fastq |
       2>> $err_file >> $log_file
    done

    exit 1
  '''
}

canu_fastas
  .join(flye_fastas, by: 0)
  .join(raven_fastas, by: 0)
  .join(filtered_fastq_cluster, by:0)
  .view()
  .set {assembled_fastas}

params.trycycler_cluster_options = ""
process trycycler_cluster {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.maxcpus

  input:
  tuple val(sample), file(fasta) from assembled_fastas

  output:
  file("contigs/newick")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version 2>> $err_file >> $log_file

    trycycler cluster !{params.trycycler_cluster_options} \
      --threads !{task.cpus} \
      --assemblies assemblies/*.fasta \
      --reads !{fastq} \
      -out_dir !{task.process} \
      2>> $err_file >> $log_file

    exit 1
  '''
}

params.trycycler_dotplot_options = ''
process trycycler_dotplot {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${cluster}"
  echo false
  cpus 1
  //  container 'staphb/trycycler:latest'

  input:
  path(cluster) from clusters_dotplot

  output:
  tuple env(cluster), file("${task.process}/*.fastq.gz") into dotplots
  file("logs/${task.process}/${cluster}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{cluster}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{cluster}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version >> $log_file

    cluster=!{cluster}

    trycycler dotplot !{params.trycycler_dotplot_options} \
      --cluster_dir !{cluster} \
      2>> $err_file >> $log_file

    exit 1
  '''
}

dotplots
  .join(clusters_reconcile, by:0)
  .join(filtered_fastq_reconcile, by: 0)
  .view()
  .set{for_reconcile}

process trycycler_reconcile {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1
  //  container 'staphb/trycycler:latest'

  input:
  tuple val(sample), file(dotplots), path(cluster), file(fastq) from for_reconcile

  output:
  tuple sample, file("${task.process}/*.fasta") into reconciled_fastas_msa
  tuple sample, file("${task.process}/*.fasta"), fastq into reconciled_fastas_partition
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version >> $log_file

    trycycler reconcile \
      --reads !{fastq} \
      --cluster_dir !{cluster} \
      --threads !{task.cpus} \
      2>> $err_file >> $log_file

    exit 1
  '''
}

params.trycycler_msa_options = ''
process trycycler_msa {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1
  //  container 'staphb/trycycler:latest'

  input:
  tuple val(sample), file(fasta) from reconciled_fastas_msa

  output:
  tuple env(sample), file("${task.process}/*.fastq.gz") into msa
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version >> $log_file

    trycycler msa !{params.trycycler_msa_options} \
      --cluster_dir !{cluster} \
      --threads !{task.cpus} \
      2>> $err_file >> $log_file

      exit 1
  '''
}

params.trycycler_partition_options = ''
process trycycler_partition {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1
  //  container 'staphb/trycycler:latest'

  input:
  tuple val(sample), file(fasta), file(fastq) from reconciled_fastas_partition

  output:
  tuple env(sample), file("${task.process}/*.fastq.gz") into partition
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version >> $log_file

    trycycler partition !{trycycler_partition_options} \
      --reads !{fastq} \
      --cluster_dirs !{cluster} \
      --threads !{task.cpus} \
      2>> $err_file >> $log_file

    exit 1
  '''
}

partition
  .join(msa, by: 0)
  .view()
  .set {for_consensus}

params.trycycler_consensus_options = ''
process trycycler_consensus {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1
  //  container 'staphb/trycycler:latest'

  input:
  tuple val(sample), file(fasta), file(fastq) from for_consensus

  output:
  tuple env(sample), file("${task.process}/*.fastq.gz") into consensus
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version >> $log_file

    trycycler consensus !{params.trycycler_consensus_options} \
      --cluster_dir !{cluster} \
      --threads !{task.cpus} \
      2>> $err_file >> $log_file

    exit 1
  '''
}

params.medaka_options = ''
process medaka {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1
  //  container 'staphb/trycycler:latest'

  input:
  tuple val(sample), file(fasta), file(fastq) from consensus

  output:
  tuple env(sample), file("${task.process}/*.fastq.gz") into medaka_fastas
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    medaka --version >> $log_file

    exit 1
  '''
}

if (params.illumina) {
  params.fastp_options = ''
  process fastp {
    publishDir "${params.outdir}", mode: 'copy'
    tag "${fastq}"
    echo false
    cpus 1
    //  container 'staphb/trycycler:latest'

    input:
    file(fastq), file(sample_key) from illumina_fastq

    output:
    tuple env(sample), file("${task.process}/*{R1,R2}.fastq.gz"), file("${task.process}/*u.fastq.gz") into clean_reads
    file("${task.process}/*.u.fastq.gz")
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
      mkdir -p !{task.process} logs/!{task.process}
      log_file=logs/!{task.process}/!{fastq}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{fastq}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null
      fastp --version >> $log_file

      sample=$(grep -e !{fastq[0]} -e !{fastq[1]} !{sample_key} | cut -f 2 -c "," | head -n 1)

      fastp !{params.fastp_options} \
        --in1 !{fastq[0]} \
        --in2 !{fastq[1]} \
        --out1 !{task.process}/${sample}_R1.fastq.gz \
        --out2 !{task.process}/${sample}_R2.fastq.gz \
        --unpaired1 !{task.process}/$sample.u.fastq.gz \
        --unpaired2 !{task.process}/$sample.u.fastq.gz

      exit 1
    '''
  }

  medaka_fastas
    .join(clean_reads, by: 0)
    .view()
    .set {for_bwa}

  process bwa {
    publishDir "${params.outdir}", mode: 'copy', pattern: "logs/bwa/*.{log,err}"
    tag "${sample}"
    echo false
    cpus params.maxcpus
    container 'staphb/bwa:latest'

    input:
    set val(sample), file(reference_genome), file(reads), file(unpaired) from for_bwa

    output:
    tuple sample, file("${task.process}/${sample}.sam"), file("${task.process}/${sample}_unpaired.sam") into bwa_sams
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
      mkdir -p !{task.process} logs/!{task.process}
      log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null
      echo "bwa $(bwa 2>&1 | grep Version )" >> $log_file
      bwa_version="bwa : "$(bwa 2>&1 | grep Version)

      # index the reference fasta file
      bwa index !{reference_genome}

      # bwa mem command
      bwa mem -t !{task.cpus} !{reference_genome} !{reads[0]} !{reads[2]} 2>> $err_file > !{task.process}/!{sample}.sam
      bwa mem -t !{task.cpus} !{reference_genome} !{unpaired} 2>> $err_file > !{task.process}/!{sample}_unpaired.sam
    '''
  }

  process sort {
    publishDir "${params.outdir}", mode: 'copy'
    tag "${sample}"
    echo false
    cpus params.maxcpus
    container 'staphb/samtools:latest'

    input:
    set val(sample), file(sam), file(unpaired_sam) from sams

    output:
    tuple sample, file("!{task.process}/${sample}.sorted{.bam,.bam.bai}"), file("!{task.process}/${sample}_unpaired.sorted{.bam,.bam.bai}") into bams
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
      mkdir -p !{task.process} logs/!{task.process}
      log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null
      samtools --version >> $log_file

      samtools sort -@ !{task.cpus} !{sam} 2>> $err_file | \
        samtools view -F 4 -o !{task.process}/!{sample}.sorted.bam --write-index 2>> $err_file >> $log_file

      samtools sort -@ !{task.cpus} !{unpaired_sam} 2>> $err_file | \
        samtools view -F 4 -o !{task.process}/!{sample}_unpaired.sorted.bam --write-index 2>> $err_file >> $log_file

    '''
  }

  bams
    .join(medaka_fastas, by:0)
    .set{ for_pilon }

  params.pilon_options = ''
  process pilon {
    publishDir "${params.outdir}", mode: 'copy'
    tag "${sample}"
    echo false
    cpus 1
    //  container 'staphb/trycycler:latest'

    input:
    tuple val(sample), file(bam), file(unpaired_bam), file(fasta) from for_pilon

    output:
    tuple env(sample), file("${task.process}/*.fastq.gz")
    file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
      mkdir -p !{task.process} logs/!{task.process}
      log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null
      pilon --version >> $log_file

      round=1

      file=!{sample}_${round}.fasta

      pilon --genome !{fasta} --frags !{bam} --unpaired !{unpaired_bam} --output $file --changes

      exit 1
    '''
  }
}

if ( params.phase2 ) {
  workflow.onComplete {
      println("Donuf Falls completed at: $workflow.complete")
      println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
  }
} else {
  workflow.onComplete {
      println("First phase of Donut Falls complete at : $workflow.complete")
      println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
  }
}
