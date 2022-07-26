#!/usr/bin/env nextflow

println("Currently using the Donut Falls workflow for use with nanopore sequencing\n")
println("Author: Erin Young")
println("email: eriny@utah.gov")
println("Version: v.0.0.20220525")
println("")

// TODO : fix commas in  numbers for summary file
// TODO : Add longQC?
// TODO : Add MuliQC?

params.outdir = workflow.launchDir + '/donut_falls'

params.maxcpus = Runtime.runtime.availableProcessors()
println("Maximum number of CPUS used in this workflow : ${params.maxcpus}")
if ( params.maxcpus < 12 ) {
  params.medcpus = params.maxcpus
} else {
  params.medcpus = 12
}

params.sequencing_summary = "${workflow.launchDir}/*sequencing_summary*txt"
Channel
  .fromPath(params.sequencing_summary, type:'file')
  .view { "Summary File : $it" }
  .ifEmpty{
    params.nanoplot_summary = false
    println("Could not find sequencing summary file! Set with 'params.sequencing_summary'")
  }
  .set { sequencing_summary }

params.reads = workflow.launchDir + '/reads'
Channel
  .fromPath("${params.reads}/*.{fastq,fastq.gz,fq,fq.gz}", type:'file')
  .ifEmpty {
    println("Could not find fastq files for nanopore sequencing. Set with 'params.reads'")
    exit 1
  }
  .map { reads -> tuple(reads.simpleName, reads ) }
  .view { "Fastq file found : ${it[0]}" }
  .into { fastq ; fastq_nanoplot }

//params.assembler = 'flye'
params.assembler = 'unicycler'
//params.assembler = 'raven'
//params.assembler = 'miniasm'

params.illumina = workflow.launchDir + '/illumina'
Channel
  .fromFilePairs("${params.illumina}/*_R{1,2}*.{fastq,fastq.gz}", size: 2 )
  .map { reads -> tuple(reads[0].replaceAll(~/_S[0-9]+_L[0-9]+/,""), reads[1]) }
  .view { "Illumina fastq files for for greater accuracy : ${it[0]}" }
  .into { illumina_fastqs ; illumina_fastqs_polishing }

params.nanoplot_summary = true
params.nanoplot_options = ''
process nanoplot_summary {
  publishDir "${params.outdir}", mode: 'copy'
  tag "nanoplot"
  cpus params.medcpus
  container 'staphb/nanoplot:latest'

  when:
  params.nanoplot_summary

  input:
  file(sequencing_summary) from sequencing_summary

  output:
  path("${task.process}")
  path("logs/${task.process}/${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/basic logs/!{task.process}
    log_file=logs/!{task.process}/!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    NanoPlot --version | tee -a $log_file $err_file > /dev/null

    NanoPlot \
      --summary !{sequencing_summary} \
      --threads !{task.cpus} \
      --outdir !{task.process} \
      2>> $err_file >> $log_file
  '''
}

params.nanoplot = true
params.nanoplot = ''
process nanoplot {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus params.medcpus
  container 'staphb/nanoplot:latest'

  when:
  params.nanoplot

  input:
  tuple val(sample), file(fastq) from fastq_nanoplot

  output:
  path("${task.process}/${sample}")
  tuple sample, env(mnrdlngt) into nanoplot_mean_read_length
  tuple sample, env(mnrdqual) into nanoplot_mean_read_quality
  tuple sample, env(mdrdlngt) into nanoplot_median_read_length
  tuple sample, env(mdrdqual) into nanoplot_median_read_quality
  tuple sample, env(numreads) into nanoplot_number_of_reads
  tuple sample, env(rdlenN50) into nanoplot_read_length_N
  tuple sample, env(totbases) into nanoplot_total_bases
  path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process}/!{sample} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    NanoPlot --version | tee -a $log_file $err_file > /dev/null

    NanoPlot \
      --fastq !{fastq} \
      --threads !{task.cpus} \
      --outdir !{task.process}/!{sample} \
      2>> $err_file >> $log_file

    mnrdlngt=$(grep "Mean read length:"    !{task.process}/!{sample}/NanoStats.txt | awk '{print $4 }' )
    mnrdqual=$(grep "Mean read quality:"   !{task.process}/!{sample}/NanoStats.txt | awk '{print $4 }' )
    mdrdlngt=$(grep "Median read length:"  !{task.process}/!{sample}/NanoStats.txt | awk '{print $4 }' )
    mdrdqual=$(grep "Median read quality:" !{task.process}/!{sample}/NanoStats.txt | awk '{print $4 }' )
    numreads=$(grep "Number of reads:"     !{task.process}/!{sample}/NanoStats.txt | awk '{print $4 }' )
    rdlenN50=$(grep "Read length N50:"     !{task.process}/!{sample}/NanoStats.txt | awk '{print $4 }' )
    totbases=$(grep "Total bases:"         !{task.process}/!{sample}/NanoStats.txt | awk '{print $3 }' )
  '''
}

params.fastp_options = ''
process fastp {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'bromberglab/fastp:latest'

  when:
  sample != null

  input:
  tuple val(sample), path(fastq) from illumina_fastqs

  output:
  tuple val(sample), path("${task.process}/${sample}_{R1,R2}.fastq.gz") into clean_reads
  path("${task.process}/*")
  path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    fastp --version >> $log_file 2>> $err_file

    fastp !{params.fastp_options} \
      --in1 !{fastq[0]} \
      --in2 !{fastq[1]} \
      --out1 !{task.process}/!{sample}_R1.fastq.gz \
      --out2 !{task.process}/!{sample}_R2.fastq.gz \
      --unpaired1 !{task.process}/!{sample}_u.fastq.gz \
      --unpaired2 !{task.process}/!{sample}_u.fastq.gz \
      2>> $err_file >> $log_file

    cp *.html !{task.process}/.
    cp *.json !{task.process}/.
  '''
}

params.filtlong_options = "--min_length 1000 --keep_percent 95"
process filtlong {
  publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/*.${workflow.sessionId}.{log,err}"
  tag "${sample}"
  cpus 1
  container 'staphb/filtlong:latest'

  input:
  tuple val(sample), file(fastq), file(short_reads) from fastq.join(clean_reads, by:0, remainder : true)

  output:
  tuple val(sample), path("${task.process}/${sample}_filtered.fastq") optional true into filtlong_fastq
  path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  if (short_reads[1] == null) {
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    filtlong --version >> $log_file

    filtlong !{params.filtlong_options} \
      !{fastq} 2>> $err_file > !{task.process}/!{sample}_filtered.fastq
  '''
  } else {
  '''
    mkdir -p !{task.process} logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    filtlong --version >> $log_file

    filtlong !{params.filtlong_options} \
      -1 !{short_reads[0]} \
      -2 !{short_reads[1]} \
      !{fastq} \
      2>> $err_file > !{task.process}/!{sample}_filtered.fastq
    '''
  }
}

process bgzip {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container 'staphb/ivar:latest'

  input:
  tuple val(sample), path(fastq) from filtlong_fastq

  output:
  tuple val(sample), path("filtlong/${fastq}.gz") optional true into filtered_fastq, filtered_fastq_medaka
  path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p filtlong logs/!{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    bgzip --version 2>> $err_file >> $log_file

    bgzip -@ {task.cpus} !{fastq}
    mv !{fastq}.gz filtlong/.
  '''
}

if ( params.assembler == 'flye' ) {
  params.flye_options = ''
  process flye {
    publishDir "${params.outdir}", mode: 'copy'
    tag "${sample}"
    cpus params.medcpus
    container 'staphb/flye:latest'
    errorStrategy 'ignore'

    input:
    tuple val(sample), path(fastq) from filtered_fastq

    output:
    tuple val(sample), path("${task.process}/${sample}.fasta") optional true into assembled_fastas
    path("${task.process}/${sample}/*")
    tuple val(sample), env(num_contigs) into num_contigs
    tuple val(sample), env(num_closed_contigs) into num_closed_contigs
    tuple val(sample), env(chr_cov) into chr_cov
    tuple val(sample), env(chr_closed) into chr_closed
    tuple val(sample), env(genome_size) into genome_size
    path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
      mkdir -p logs/!{task.process} !{task.process}/!{sample}
      log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null
      flye --version 2>> $err_file >> $log_file

      flye !{params.flye_options} \
        --nano-raw !{fastq} \
        --threads !{task.cpus} \
        --out-dir !{task.process}/!{sample} \
        2>> $err_file >> $log_file

      if [ -f "!{task.process}/!{sample}/assembly.fasta" ]
      then
        cp !{task.process}/!{sample}/assembly.fasta !{task.process}/!{sample}.fasta
        num_contigs=$(grep ">" !{task.process}/!{sample}/assembly.fasta | wc -l )
        genome_size=$(grep -v ">" !{task.process}/!{sample}/assembly.fasta | wc -c )
        num_closed_contigs=$(cut -f 4 !{task.process}/!{sample}/assembly_info.txt | grep -w "Y" | wc -l )
        chr_cov=$(head -n 2 !{task.process}/!{sample}/assembly_info.txt | tail -n 1 | cut -f 3)
        chr_closed=$(head -n 2 !{task.process}/!{sample}/assembly_info.txt | tail -n 1 | cut -f 4)
      fi

      if [ -z "$num_contigs" ]        ; then num_contigs="0"        ; fi
      if [ -z "$num_closed_contigs" ] ; then num_closed_contigs="0" ; fi
      if [ -z "$chr_cov" ]            ; then chr_cov="0"            ; fi
      if [ -z "$genome_size" ]        ; then genome_size="0"        ; fi
      if [ -z "$chr_closed" ]         ; then chr_closed="0"         ; fi
    '''
  }
} else if ( params.assembler == 'miniasm' ) {
  params.miniasm_and_minipolish_options = ''
  process miniasm_and_minipolish {
    publishDir "${params.outdir}", mode: 'copy'
    tag "${sample}"
    cpus params.medcpus
    container 'staphb/minipolish:latest'

    input:
    tuple val(sample), path(fastq) from filtered_fastq

    output:
    tuple val(sample), path("${task.process}/${sample}/*gfa") into miniasm_gfa
    path("${task.process}/${sample}/*")
    path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
      mkdir -p !{task.process}/!{sample} logs/!{task.process}
      log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null
      echo "miniasm version : $(miniasm -V)" 2>> $err_file >> $log_file
      minimap2 --version 2>> $err_file >> $log_file
      minipolish --version 2>> $err_file >> $log_file

      miniasm_and_minipolish.sh \
        !{fastq} \
        !{task.cpus} \
        2>> $err_file > !{task.process}/!{sample}/!{sample}.gfa
    '''
  }

  params.any2fasta_options = ''
  process any2fasta {
    publishDir "${params.outdir}", mode: 'copy'
    tag "${sample}"
    cpus 1
    container 'staphb/any2fasta:latest'

    input:
    tuple val(sample), path(gfa) from miniasm_gfa

    output:
    tuple val(sample), path("${task.process}/${sample}.fasta") into assembled_fastas
    path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
      mkdir -p !{task.process}/!{sample} logs/!{task.process}
      log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null
      any2fasta -v 2>> $err_file >> $log_file

      any2fasta !{params.any2fasta_options} \
        !{gfa} \
        2>> $err_file \
        > !{task.process}/!{sample}.fasta
    '''
  }
} else if ( params.assembler == 'raven' ) {
  params.raven_options = '--polishing-rounds 2'
  process raven {
    publishDir "${params.outdir}", mode: 'copy'
    tag "${sample}"
    cpus params.medcpus
    container 'staphb/raven:latest'

    input:
    tuple val(sample), path(fastq) from filtered_fastq

    output:
    tuple val(sample), path("${task.process}/${sample}/${sample}.fasta") into raven_fastas
    path("${task.process}/${sample}/*")
    path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

    shell:
    '''
      mkdir -p !{task.process}/!{sample} logs/!{task.process}
      log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
      err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

      # time stamp + capturing tool versions
      date | tee -a $log_file $err_file > /dev/null
      raven --version 2>> $err_file >> $log_file

      raven !{params.raven_options} \
        --threads !{task.cpus} \
        --graphical-fragment-assembly !{task.process}/!{sample}/!{sample}.gfa \
        !{fastq} \
        2>> $err_file \
        > !{task.process}/!{sample}/!{sample}.fasta
    '''
  }
} else {
  assembled_fastas = Channel.empty()
}

params.medaka_options = ''
process medaka {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 6
  container 'ontresearch/medaka:latest'

  input:
  tuple val(sample), path(fasta), path(fastq) from assembled_fastas.join(filtered_fastq_medaka, by:0 )

  output:
  path("${task.process}/${sample}/")
  tuple val(sample), path("${task.process}/${sample}/${sample}_medaka_consensus.fasta") into medaka_fastas
  path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process} !{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    medaka --version >> $log_file

    medaka_consensus !{params.medaka_options} \
      -i !{fastq} \
      -d !{fasta} \
      -o !{task.process}/!{sample} \
      -t 2

    cp !{task.process}/!{sample}/consensus.fasta !{task.process}/!{sample}/!{sample}_medaka_consensus.fasta
  '''
}

round = Channel.of( 0 )
changes = Channel.of(10000)

medaka_fastas
  .join(illumina_fastqs_polishing, by:0 )
  .combine(round)
  .combine(changes)
  .set{ round_1 }

new_rounds = Channel.create()

round_1
  .mix(new_rounds)
  .view()
  .set{ for_polca }

params.max_polish_rounds = 10
params.polca_options = ''
process polca {
  publishDir "${params.outdir}/${task.process}/${sample}/round_${round}", mode: 'copy', pattern: "*{err,batches,success,vcf,names,report,PolcaCorrected.fa}"
  publishDir "${params.outdir}", mode: 'copy', pattern: "logs/${task.process}/${sample}_${round}.${workflow.sessionId}.{log,err}"
  publishDir "${params.outdir}", mode: 'copy', pattern: "${task.process}/${sample}/${sample}_*.fasta"
  tag "${sample} : round ${round} : changes ${changes}"
  cpus params.medcpus
  container 'staphb/masurca:latest'

  input:
  tuple val(sample), path(fasta), path(fastq), val(round), val(changes) from for_polca

  output:
  tuple val(sample), path("next/${sample}_${round}.fasta"), path(fastq), env(next_round), env(change_test) optional true into new_rounds
  path("*")
  path("${task.process}/${sample}/${sample}_{${round},final}.fasta") optional true
  path("logs/${task.process}/${sample}_${round}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process} !{task.process}/!{sample} next
    log_file=logs/!{task.process}/!{sample}_!{round}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}_!{round}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null
    masurca --version >> $log_file

    polca.sh !{params.polca_options} \
      -r '!{fastq}' \
      -a !{fasta} \
      -t !{task.cpus}

    if [ -f "!{fasta}.report" ]
    then
      sub_err=$(grep "Substitution Errors" !{fasta}.report | awk '{print $3}')
      del_err=$(grep "Deletion Errors" !{fasta}.report | awk '{print $3}')
      if [ -z "$sub_err" ] ; then sub_err=0 ; fi
      if [ -z "$del_err" ] ; then del_err=0 ; fi
      change_test=$(( sub_err + del_err ))
    else
      echo "WARNING : report not found" >> $err_file
      sub_err=0
      del_err=0
      change_test=0
    fi

    echo "The number of changes from this round was $change_test" >> $log_file
    echo "$sub_err changes were from substitution errors" >> $log_file
    echo "$del_err changes were from insertion/deletion errors" >> $log_file

    if [ "$change_test" -lt "1" ]
    then
      next_round=1000
      cp !{fasta}.PolcaCorrected.fa !{task.process}/!{sample}/!{sample}_!{round}.fasta
      cp !{fasta}.PolcaCorrected.fa !{task.process}/!{sample}/!{sample}_final.fasta
    elif [ "$change_test" -eq "!{changes}" ] && [ "$change_test" -lt "10000" ]
    then
      next_round=1000
      cp !{fasta}.PolcaCorrected.fa !{task.process}/!{sample}/!{sample}_!{round}.fasta
      cp !{fasta}.PolcaCorrected.fa !{task.process}/!{sample}/!{sample}_final.fasta
    else
      cp !{fasta}.PolcaCorrected.fa !{task.process}/!{sample}/!{sample}_!{round}.fasta

      next_round=$(( !{round} + 1 ))
      if [ "$next_round" -le "!{params.max_polish_rounds}" ] ; then cp !{fasta}.PolcaCorrected.fa next/!{sample}_!{round}.fasta ; fi
    fi
  '''
}

num_contigs
  .join(num_closed_contigs, remainder: true, by: 0)
  .join(chr_closed, remainder: true, by: 0)
  .join(chr_cov, remainder: true, by: 0)
  .join(genome_size, remainder: true, by: 0)
  .join(nanoplot_mean_read_length, remainder: true, by: 0)
  .join(nanoplot_mean_read_quality, remainder: true, by: 0)
  .join(nanoplot_median_read_length, remainder: true, by: 0)
  .join(nanoplot_median_read_quality, remainder: true, by: 0)
  .join(nanoplot_number_of_reads, remainder: true, by: 0)
  .join(nanoplot_read_length_N, remainder: true, by: 0)
  .join(nanoplot_total_bases, remainder: true, by: 0)
  .set { results }

params.summary = true
process summary {
  publishDir "${params.outdir}", mode: 'copy'
  tag "${sample}"
  cpus 1
  container ''

  input:
  tuple val(sample),
    val(num_contigs),
    val(num_closed_contigs),
    val(chr_closed),
    val(chr_cov),
    val(genome_size),
    val(nanoplot_mean_read_length),
    val(nanoplot_mean_read_quality),
    val(nanoplot_median_read_length),
    val(nanoplot_median_read_quality),
    val(nanoplot_number_of_reads),
    val(nanoplot_read_length_N),
    val(nanoplot_total_bases) from results

  when:
    params.summary

  output:
  path("${task.process}/${sample}.summary.tsv") into summary
  path("logs/${task.process}/${sample}.${workflow.sessionId}.{log,err}")

  shell:
  '''
    mkdir -p logs/!{task.process} !{task.process}
    log_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.log
    err_file=logs/!{task.process}/!{sample}.!{workflow.sessionId}.err

    # time stamp + capturing tool versions
    date | tee -a $log_file $err_file > /dev/null

    header="sample\tnumber_of_contigs\tclosed_contigs\tchr_cov\tchr_closed\tgenome_size"
    result="!{sample}\t!{num_contigs}\t!{num_closed_contigs}\t!{chr_cov}\t!{chr_closed}\t!{genome_size}"

    if [ "!{params.nanoplot}" != "false" ]
    then
      header="$header\tmean_read_length\tmean_read_quality\tmedian_read_length\tmediean_read_quality\tnumber_of_reads\tread_length_N50\ttotal_bases"
      result="$result\t!{nanoplot_mean_read_length}\t!{nanoplot_mean_read_quality}\t!{nanoplot_median_read_length}\t!{nanoplot_median_read_quality}\t!{nanoplot_number_of_reads}\t!{nanoplot_read_length_N}\t!{nanoplot_total_bases}"
    fi

    echo -e "$header" > !{task.process}/!{sample}.summary.tsv
    echo -e "$result" >> !{task.process}/!{sample}.summary.tsv
  '''
}

summary
  .collectFile(name: "summary.tsv",
    keepHeader: true,
    sort: true,
    skip: 1,
    storeDir: "${params.outdir}")

workflow.onComplete {
  println("Donut Falls complete at : $workflow.complete")
  println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
  }
