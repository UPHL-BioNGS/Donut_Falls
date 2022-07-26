#!/usr/bin/env nextflow

println("Currently using the Donut Falls workflow for use with nanopore sequencing\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v.20220430")
println("")

params.outdir = workflow.launchDir + '/trycycler_cluster'

params.maxcpus = Runtime.runtime.availableProcessors()
println("Maximum number of CPUS used in this workflow : ${params.maxcpus}")
if ( params.maxcpus < 12 ) {
  params.medcpus = params.maxcpus
} else {
  params.medcpus = 12
}

params.cluster_directory = workflow.launchDir + '/donut_falls/trycycler_cluster'
if ( params.phase2 ) {
  println("Beginning second phase of workflow.")
  Channel
    .fromPath("${params.cluster_directory}/*", type:'dir')
    .map { it -> tuple(it.getName(), it ) }
    .view { "Samples with clusters to reconcile : ${it[0]}" }
    .set { clusters }
} else {
  clusters = Channel.empty()
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
    .combine(sample_key_illumina)
    .view {"Illumina fastq files for pilon polishing : ${it[0]}"}
    .set { illumina_fastqs }
} else {
  illumina_fastqs = Channel.empty()
}

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
  set val(sample), file(fastq) from combined_fastq

  output:
  tuple sample, file("${task.process}/${sample}_filtered.fastq") into filtered_fastq_subsample, filtered_fastq_cluster, filtered_fastq_reconcile, filtered_fastq_partician, filtered_fastq_medaka
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

  when:
  params.phase2 == false

  input:
  set val(sample), file(fastq) from filtered_fastq_subsample

  output:
  tuple sample, file("${task.process}/${sample}/miniasm/sample*.fastq") into subsampled_miniasm_fastq
  tuple sample, file("${task.process}/${sample}/flye/sample*fastq") into subsampled_flye_fastq
  tuple sample, file("${task.process}/${sample}/raven/sample*fastq") into subsampled_raven_fastq
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/!{sample}/{miniasm,flye,raven} logs/!{task.process}
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

    ls !{task.process}/!{sample}/*fastq  | head -n !{subsample_count} | xargs -I {} cp {} !{task.process}/!{sample}/miniasm/.
    ls !{task.process}/!{sample}/*fastq  | tail -n !{subsample_count} | xargs -I {} cp {} !{task.process}/!{sample}/flye/.
    ls !{task.process}/!{sample}/*fastq  | tail -n !{twice_count} | head -n !{subsample_count} | xargs -I {} cp {} !{task.process}/!{sample}/raven/.
  '''
}

params.flye_options = ''
process flye {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.medcpus
  container 'staphb/flye:latest'

  input:
  set val(sample), file(fastq) from subsampled_flye_fastq

  output:
  tuple sample, file("${task.process}/${sample}/*_flye.fasta") into flye_fastas
  file("${task.process}/${sample}/*/{assembly_graph.gfa,assembly_graph.gv,assembly_info.txt,flye.log,params.json}")
  file("${task.process}/${sample}/*/00-assembly/draft_assembly.fasta")
  file("${task.process}/${sample}/*/10-consensus/{chunks.fasta.fai,consensus.fasta}")
  file("${task.process}/${sample}/*/20-repeat/*")
  file("${task.process}/${sample}/*/22-plasmids/*")
  file("${task.process}/${sample}/*/30-contigger/*")
  file("${task.process}/${sample}/*/40-polishing/*")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err
    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    flye --version 2>> $err_file >> $log_file
    for fastq in !{fastq}
    do
      name=$(echo $fastq | sed 's/.fastq//g')
      mkdir -p !{task.process}/!{sample}/$name
      flye !{params.flye_options} \
        --nano-raw $fastq \
        --threads !{task.cpus} \
        --plasmids \
        --out-dir !{task.process}/!{sample}/$name \
        2>> $err_file >> $log_file
        cp !{task.process}/!{sample}/$name/assembly.fasta !{task.process}/!{sample}/${name}_!{task.process}.fasta
    done
  '''
}

params.miniasm_and_minipolish_options = ''
process miniasm_and_minipolish {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.medcpus
  //container 'staphb/minipolish:latest'

  input:
  set val(sample), val(genome_size), file(fastq) from subsampled_miniasm_fastq

  output:
  tuple sample, file("${task.process}/${sample}/*gfa") into miniasm_gfa
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err
    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    echo "miniasm version : $(miniasm -V)" 2>> $err_file >> $log_file
    miniasm_and_minipolish.sh --version 2>> $err_file >> $log_file
    for fastq in !{fastq}
    do
      name=$(echo $fastq | sed 's/.fastq//g')
      miniasm_and_minipolish.sh \
        $fastq \
        !{task.cpus} \
        2>> $err_file \
        > !{task.process}/!{sample}/$name.gfa
    done
  '''
}

params.any2fasta_options = ''
process any2fasta {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus 1
  //container 'staphb/any2fasta:latest'

  input:
  tuple val(sample), file(gfa) from miniasm_gfa

  output:
  tuple sample, file("${task.process}/${sample}/*_miniasm.fasta") into miniasm_fastas
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err
    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    any2fasta -v 2>> $err_file >> $log_file
    for gfa in !{gfa}
    do
      name=$(echo $gfa | sed 's/.gfa//g')
      any2fasta !{params.any2fasta_options} \
        $gfa \
        2>> $err_file \
        > !{task.process}/!{sample}/${name}_miniasm.fasta
    done
  '''
}

params.raven_options = ''
params.raven_polishing_rounds = '2'
process raven {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.medcpus
  container 'staphb/raven:latest'

  input:
  set val(sample), file(fastq) from subsampled_raven_fastq

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

miniasm_fastas
  .join(flye_fastas, by: 0)
  .join(raven_fastas, by: 0)
  .join(filtered_fastq_cluster, by:0)
  .set {assembled_fastas}

params.trycycler_cluster_options = ""
process trycycler_cluster {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  echo false
  cpus params.medcpus
//  container 'staphb/trycycler:latest'

  when:
  params.phase2 == false

  input:
  set val(sample), file(miniasm), file(flye), file(raven), file(fastq) from assembled_fastas

  output:
  file("${task.process}/${sample}/contigs.{newick,phylip}")
  file("${task.process}/${sample}/cluster*/*contigs/*.fasta")
  tuple sample, path("${task.process}/${sample}") into inital_clusters
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err
    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version 2>> $err_file >> $log_file
    trycycler cluster !{params.trycycler_cluster_options} \
      --threads !{task.cpus} \
      --assemblies *.fasta \
      --reads !{fastq} \
      --out_dir !{task.process}/!{sample} \
      2>> $err_file >> $log_file
  '''
}

// TODO : move to before intermission
params.trycycler_dotplot_options = ''
process trycycler_dotplot {
  publishDir "${params.outdir}", mode: 'copy', pattern: '*.{log,err}'
  publishDir "${params.outdir}/phase2", mode: 'copy', pattern: '*dotplots.png'
  tag "${cluster}"
  echo false
  cpus 1
  // errorStrategy 'finish'
  // container 'staphb/trycycler:latest'

  input:
  set val(sample), path(cluster) from inital_clusters

  output:
  path("${cluster}/*/dotplots.png")
  file("logs/${task.process}/${cluster}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process}
    log_file=logs/!{task.process}/!{cluster}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{cluster}.!{workflow.sessionId}.err
    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version >> $log_file
    for cluster in $(ls !{cluster}/cluster* -d )
    do
      name=$(echo $cluster | cut -f 2 -d "/")
      trycycler dotplot !{params.trycycler_dotplot_options} \
        --cluster_dir $cluster \
        2>> $err_file >> $log_file
    done
  '''
}

clusters
  .join(filtered_fastq_reconcile, by: 0)
  .set{for_reconcile}

params.trycycler_reconcile_options = ''
params.trycycler_reconcile_minimum = '2'
process trycycler_reconcile {
  publishDir "${params.outdir}", mode: 'copy', pattern: '*.{log,err}'
  //publishDir "${params.outdir}/phase2/${sample}", mode: 'copy', pattern: '*2_all_seqs.fasta'
  tag "${sample}"
  echo true
  cpus params.medcpus
  //errorStrategy 'finish'
  //container 'staphb/trycycler:latest'

  input:
  set val(sample), path(cluster), file(fastq) from for_reconcile

  output:
  tuple sample, path("${task.process}/${sample}/*") into reconciled_fastas
  file("trycycler_reconcile/${sample}/*/2_all_seqs.fasta")
  file("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/!{cluster} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err
    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    trycycler --version >> $log_file
    rsync -rh !{cluster}/ !{task.process}/!{cluster}/.
    err_flag=''
    for cluster in $(ls !{task.process}/!{cluster}/cluster* -d )
    do
      num_fasta=$(ls $cluster/1_contigs/*.fasta | wc -l)
      if [ "$num_fasta" -lt "!{params.trycycler_reconcile_minimum}" ]
      then
        rm -rf $cluster
      else
        trycycler reconcile !{params.trycycler_reconcile_options} \
          --reads !{fastq} \
          --cluster_dir $cluster \
          --threads !{task.cpus} \
          2>> $cluster.err >> $log_file
        if [ -f "$cluster/2_all_seqs.fasta" ]
        then
          cat $cluster.err >> $err_file
        else
          err_flag='1'
          echo "User input required:"
          grep Error -A 2 $cluster.err
          cat $cluster.err >> $err_file
        fi
        rm $cluster.err
      fi
    done
    if [ -n "$err_flag" ]
    then
      exit 1
    fi
  '''
}

reconciled_fastas
  .transpose()
  .filter( it -> it[1] =~ /cluster/ )
  .set { reconciled_fastas_cluster }

params.trycycler_msa_options = ''
process trycycler_msa {
  publishDir "${params.outdir}", mode: 'copy', pattern: '*.{log,err}'
  tag "${sample}_${trycycler}"
  echo false
  cpus params.medcpus
  //container 'staphb/trycycler:latest'

  input:
  set val(sample), path(trycycler) from reconciled_fastas_cluster

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

msa
  .combine(filtered_fastq_partician, by: 0)
  .set{ for_partition }

params.trycycler_partition_options = ''
process trycycler_partition {
  publishDir "${params.outdir}", mode: 'copy', pattern: '*.{log,err}'
  tag "${sample}_${cluster}"
  echo false
  cpus params.medcpus
  //container 'staphb/trycycler:latest'

  input:
  set val(sample), val(cluster), path(trycycler), file(fastq) from for_partition

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
  publishDir "${params.outdir}", mode: 'copy', pattern: '*.{log,err}'
  publishDir "${params.outdir}", mode: 'copy', pattern: '*7_final_consensus.fasta'
  publishDir "${params.outdir}", mode: 'copy', pattern: 'final_consensus.fasta'
  tag "${sample}_${cluster}"
  echo false
  cpus params.medcpus
  //container 'staphb/trycycler:latest'

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
  echo false
  cpus 1
  //container 'staphb/trycycler:latest'

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
  echo false
  cpus 1
  container 'staphb/medaka:latest'

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
  echo false
  cpus 1
  //container 'staphb/medaka:latest'

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
