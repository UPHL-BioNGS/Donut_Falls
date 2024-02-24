#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// read but ignored most things from
// https://carpentries-incubator.github.io/Pipeline_Training_with_Nextflow/07-Nextflow_Best_Practice/index.html

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Greetings!

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

println('')
println(' __                    ___          ')
println(' ) ) _   _      _)_    )_ _   ) ) _ ')
println('/_/ (_) ) ) (_( (_    (  (_( ( ( (  ')
println('                                 _) ')
println('')

println('Currently using the Donut Falls workflow for use with nanopore sequencing')
println('Author: Erin Young')
println('email: eriny@utah.gov')
println("Version: ${workflow.manifest.version}")
println('')

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Setting default param values

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####


params.config_file        = false
if (params.config_file) {
  def src = new File("${workflow.projectDir}/configs/donut_falls_config_template.config")
  def dst = new File("${workflow.launchDir}/edit_me.config")
  dst << src.text
  println("A config file can be found at ${workflow.launchDir}/edit_me.config")
  exit 0
}

params.sequencing_summary = ''
params.sample_sheet       = ''
params.assembler          = 'flye'
params.outdir             = 'donut_falls'
params.test               = false
params.ontime             = false

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Checking params

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

def paramCheck(keys) {
  set_keys = [
    "outdir",
    "sample_sheet",
    "sequencing_summary",
    "assembler",
    "test",
    "config_file",
    "ontime"]

  for(key in keys){
    if (key !in set_keys){
      println("FATAL: ${key} isn't a supported param!")
      println("Supported params: ${set_keys}")
      exit 1
    }
  }
}

paramCheck(params.keySet())

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Input files

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

Channel
  .fromPath(params.sequencing_summary, type:'file')
  .view { "Summary File : $it" }
  .set { ch_sequencing_summary }

// using a sample sheet with the column header of 'sample,fastq,fastq_1,fastq_2'
// sample  = meta.id
// fastq   = nanopore fastq file
// fastq_1 = illumina fastq file
// fastq_2 = illumina fastq file

Channel
  .fromPath("${params.sample_sheet}", type: "file")
  .splitCsv( header: true, sep: ',' )
  .map { it ->
    meta = [id:it.sample] 
    tuple( meta,
      "${it.fastq}",
      "${it.fastq_1}",
      "${it.fastq_2}" )
  }
  .branch { it ->
    sr:  it[2] != it[3]
    other: true
  }
  .set{ch_precheck}

ch_precheck.sr
  .map { it -> tuple(it[0], file(it[1], checkIfExists: true), [file(it[2], checkIfExists: true), file(it[3], checkIfExists: true)])}
  .mix(ch_precheck.other.map{ it -> tuple(it[0], file(it[1], checkIfExists: true), null)})
  .set{ch_input_files}

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Processes

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

process bandage {
  tag           "${meta.id}"
  label         "process_low"
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/biocontainers/bandage:0.8.1--hc9558a2_2'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'

  input:
  tuple val(meta), file(gfa)

  output:
  path "bandage/*"    , emit: files
  path "bandage/*.png", emit: png
  path "versions.yml" , emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p bandage

  Bandage image ${gfa} bandage/${prefix}.png ${args}
  Bandage image ${gfa} bandage/${prefix}.svg ${args}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    bandage: \$(echo \$(Bandage --version 2>&1) | sed 's/^.*Version: //; s/ .*\$//')
  END_VERSIONS
  exit 1
  """
}

process busco {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/busco:5.6.1-prok-bacteria_odb10_2024-01-08'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '45m'

  input:
  tuple val(meta), file(fasta)

  output:
  path("busco/*/*")
  path("busco/*/short_summary*.txt"), optional: true, emit: summary

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: '--offline'
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  busco ${args} \
    -m genome \
    -i ${fasta} \
    -o busco/${prefix} \
    --cpu ${task.cpus} \
    -l /busco_downloads/lineages/bacteria_odb10

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    busco: \$( busco --version 2>&1 | sed 's/^BUSCO //' )
  END_VERSIONS
  exit 1
  """
}

process bwa {
  tag           "${meta.id}"
  label         'process_high'
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/bwa:0.7.17'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '2h'

  input:
  tuple val(meta), file(fasta), file(fastq)

  output:
  tuple val(meta), file(fasta), file("bwa/${sample}_{1,2}.sam"), emit: sam

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p bwa

  bwa index ${fasta}
  bwa mem -t ${task.cpus} -a ${fasta} ${fastq[0]} > bwa/${prefix}_1.sam
  bwa mem -t ${task.cpus} -a ${fasta} ${fastq[1]} > bwa/${prefix}_2.sam

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
  END_VERSIONS
  exit 1
  """
}

process circulocov {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/uphl/circulocov:0.1.20240104-2024-02-21'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '1h'

  input:
  tuple val(meta), file(fastqs), file(contigs)

  output:
  tuple val(meta), file("circulocov/*/*sr.bam*"), emit: bam
  path "circulocov/*/overall_summary.txt", emit: collect
  path "circulocov/*/*", emit: everything
  path "circulocov/*/fastq/*", emit: fastq
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: '-a'
  def prefix = task.ext.prefix ?: "${meta.id}"
  def reads  = fastqs.join(" ")
  """
  mkdir -p circulocov/${prefix}
  
  circulocov ${args} \
    --threads ${task.cpus} \
    --genome ${contigs} \
    --illumina ${reads} \
    --out circulocov/${prefix} \
    --sample ${prefix}
  
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    circulocov: \$(circulocov -v)
  END_VERSIONS
  exit 1
  """
}

process copy {
  tag           "${meta.id}"
  label         "process_single"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/multiqc:1.19'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'

  input:
  tuple val(meta), file(fasta)

  output:
  path "consensus/"
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir consensus

  for fasta in !{fasta}
  do
    cat \$fasta | sed 's/_length/ length/g' | sed 's/_circular/ circular/g' | sed 's/_polypolish//g' > consensus/$fasta
  done
  exit 1
  """
}

process dnaapler {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/dnaapler:0.7.0'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '1h'

  input:
  tuple val(meta), file(reads)

  output:
  tuple val(meta), file("bbduk/*_rmphix_R{1,2}.fastq.gz"),  emit: fastq
  path "bbduk/*",                                           emit: files
  path "bbduk/*.phix.stats.txt",                            emit: stats
  path "logs/${task.process}/*.log",  emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  dnaapler --help

  dnaapler chromosome --input chromosome.fasta --output dnaapler_chr

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    dnaapler: "\$(dnaapler --version 2>&1 "
  END_VERSIONS
  exit 1
  """
}

process download {
  tag           "Downloading subset15000"
  label         "process_single"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/multiqc:1.19'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '1h'

  output:
  tuple val("subset15000"), file("nfcore_subset15000.fa.gz"), emit: fastq

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  wget -q https://github.com/nf-core/test-datasets/blob/23f5b889e4736798c8692e9b92810d9a3e37ee97/nanopore/subset15000.fq.gz?raw=true -O nfcore_subset15000.fa.gz

  exit 1
  """
}

process great_dataset {
  tag           "Downloading the great dataset"
  label         "process_single"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/multiqc:1.19'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '1h'

  output:
  tuple val("great_dataset"), file("reads.fastq.gz"), emit: fastq

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  wget -q https://bridges.monash.edu/ndownloader/files/23754659 -O great_dataset.tar.gz
  tar -xvf great_dataset.tar.gz

  exit 1
  """
}

process dragonflye {
  tag           "${meta.id}"
  label         "process_high"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/dragonflye:1.0.14'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10h'

  input:
  tuple val(meta), file(fastq)

  output:
  tuple val(meta), file("dragonflye/*/*_dragonflye.fasta"), optional: true,  emit: fasta
  tuple val(meta), file("dragonflye/*/*_dragonflye.gfa"),   optional: true,  emit: gfa
  path "dragonflye/*/*_assembly_info.tsv", emit: summary
  path "dragonflye/*/*", emit: everything
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p dragonflye

  dragonflye ${args} \
    --reads ${fastq} \
    --cpus ${task.cpus} \
    --outdir dragonflye/${prefix} \
    --prefix ${prefix}

  # renaming final files
  if [ -f "dragonflye/${prefix}/flye-unpolished.gfa" ] ; then cp dragonflye/${prefix}/flye-unpolished.gfa dragonflye/${prefix}/${prefix}_dragonflye.gfa ; fi
  if [ -f "dragonflye/${prefix}/flye.fasta" ] ; then cp dragonflye/${prefix}/flye.fasta dragonflye/${prefix}/${prefix}_dragonflye.fasta ; fi

  # getting a summary file
  head -n 1 dragonflye/${prefix}/flye-info.txt | awk '{print "sample\\t" \$0}' > dragonflye/${prefix}/${prefix}_assembly_info.tsv
  tail -n+2 dragonflye/${prefix}/flye-info.txt | awk -v sample=${prefix} '{print sample "\\t" \$0}' >> dragonflye/${prefix}/${prefix}_assembly_info.tsv
  
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    dragonflye: \$(dragonflye --version | sed 's/^.*dragonflye //' )
  END_VERSIONS

  exit 1
  """
}

process fastp {
  tag           "${meta.id}"
  label         "process_low"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/fastp:0.23.2'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'

  input:
  tuple val(meta), file(reads)

  output:
  tuple val(meta), file("fastp/*_fastp_{R1,R2}.fastq.gz"), emit: reads
  path "fastp/*", emit: everything
  path "fastp/*_fastp.json", emit: summary
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p fastp

  fastp ${args} \
    --in1 ${reads[0]} \
    --in2 ${reads[1]} \
    --out1 fastp/${prefix}_fastp_R1.fastq.gz \
    --out2 fastp/${prefix}_fastp_R2.fastq.gz \
    --unpaired1 fastp/${prefix}_u.fastq.gz \
    --unpaired2 fastp/${prefix}_u.fastq.gz \
    -h fastp/${prefix}_fastp.html \
    -j fastp/${prefix}_fastp.json

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    fastp: \$(fastp --version | sed -e 's/fastp //g')
  END_VERSIONS

  exit 1
  """
}

process filtlong {
  tag           "${meta.id}"
  label         "process_low"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/filtlong:0.2.1'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '30m'

  input:
  tuple val(meta), file(fastq), file(short_reads)

  output:
  tuple val(meta), file("filtlong/${sample}_filtered.fastq.gz"), optional: true, emit: fastq
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  if (short_reads[1] == null) {
   """
    mkdir -p filtlong

    filtlong ${args} \
      ${fastq} \
      | gzip |
      > filtlong/${prefix}_filtered.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      filtlong: \$( filtlong --version | sed -e "s/Filtlong v//g" )
    END_VERSIONS

    exit 1
    """
  } else {
    """
    mkdir -p filtlong

    filtlong ${args} \
      -1 ${short_reads[0]} \
      -2 ${short_reads[1]} \
      ${fastq} \
      | gzip \
      > filtlong/${prefix}_filtered.fastq.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      filtlong: \$( filtlong --version | sed -e "s/Filtlong v//g" )
    END_VERSIONS

    exit 1
    """
  }
}

process flye {
  tag           "${meta.id}"
  label         "process_high"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/flye:2.9.2'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10h'

  input:
  tuple val(meta), file(fastq)

  output:
  tuple val(meta), file("flye/*/*_flye.fasta"), optional: true,  emit: fasta
  tuple val(meta), file("flye/*/*_flye.gfa"), optional: true,  emit: gfa
  path "flye/*/*_assembly_info.tsv", emit: summary
  path "flye/*/*", emit: everything
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p flye/${prefix}

  flye ${args} \
    --nano-raw ${fastq} \
    --threads ${task.cpus} \
    --out-dir flye/${prefix}

  # renaming final files
  if [ -f "flye/${prefix}/assembly.fasta" ]     ; then cp flye/${prefix}/assembly.fasta     flye/${prefix}/${prefix}_flye.fasta ; fi
  if [ -f "flye/${prefix}/assembly_graph.gfa" ] ; then cp flye/${prefix}/assembly_graph.gfa flye/${prefix}/${prefix}_flye.gfa   ; fi

  # getting a summary file
  head -n 1 flye/${prefix}/assembly_info.txt | awk '{print "sample\\t" \$0}' > flye/${prefix}/${prefix}_assembly_info.tsv
  tail -n+2 flye/${prefix}/assembly_info.txt | awk -v sample=${prefix} '{print sample "\\t" \$0}' >> flye/${prefix}/${prefix}_assembly_info.tsv

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    flye: \$( flye --version )
  END_VERSIONS

  exit 1
  """
}

process gfastats {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/gfastats:1.3.6'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'

  input:
  tuple val(meta), file(gfa)

  output:
  tuple val(meta), file(gfa), file("gfastats/*_gfastats_summary.csv"), emit: stats
  path "gfastats/*_gfastats_summary.csv", emit: collect
  path "gfastats/*", emit: everything
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p gfastats

  gfastats \
    ${gfa} \
    ${args} \
    --threads ${task.cpus} \
    --tabular \
    --seq-report \
    > gfastats/${prefix}_gfastats.txt

  head -n 1 gfastats/${prefix}_gfastats.txt | tr "\\t" "," | awk '{print "sample," $0 "circular" }' > gfastats/${prefix}_gfastats_summary.csv
  tail -n+2 gfastats/${prefix}_gfastats.txt | tr "\\t" "," | awk -v sample=${prefix} '{print sample "," $0 }' >> gfastats/${prefix}_gfastats_summary.csv
  
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    gfastats: \$( gfastats -v | sed '1!d;s/.*v//' )
  END_VERSIONS

  exit 1
  """
}

process gfa_to_fasta {
  tag           "${meta.id}"
  label         "process_low"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/gfastats:1.3.6'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'

  input:
  tuple val(meta), file(gfa), file(stats)

  output:
  tuple val(meta), file("fastas/*"), emit: fasta
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p fasta

  # get gfa file and gfa stats file
  # get circular and linear fasta files
  # figure out a chromosome
  # figure out plasmids
  # fix headers

  # in python? in shell?

  exit 1
  """
}

process hybracter {
  tag           "${meta.id}"
  label         "process_high"
  publishDir    params.outdir, mode: 'copy'
  container     'quay.io/biocontainers/hybracter:0.6.0--pyhdfd78af_0'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10h'

  input:
  tuple val(meta), file(reads)

  output:
  tuple val(meta), file("bbduk/*_rmphix_R{1,2}.fastq.gz"),  emit: fastq
  path "bbduk/*",                                           emit: files
  path "bbduk/*.phix.stats.txt",                            emit: stats
  path "logs/${task.process}/*.log",  emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  hybracter -h
  
  exit 1

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    hybracter: "\$(hybracter --version 2>&1 | grep -v java | grep version | awk '{print \$NF}')"
  END_VERSIONS
  exit 1
  """
}

process pypolca {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/pypolca:0.3.1'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '30m'
  
  input:
  tuple val(meta), file(fasta), file(fastq)

  output:
  tuple val(meta), file("polca/*/*.fasta"), optional: true, emit: fasta
  path "polca/*/*", emit: everything                                                               emit: directory
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p pypolca/${prefix}

  pypolca run ${args}\
    -a ${fasta} \
    -1 ${fastq[0]} \
    -2 ${fastq[1]} \
    -t ${task.cpus} \
    -o pypolca/${prefix}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    pypolca: \$(pypolca --version)
  END_VERSIONS

  exit 1
  """
}

process medaka {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'ontresearch/medaka:v1.11.3'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '30m'

  input:
  tuple val(meta), path(fasta), path(fastq)

  output:
  tuple val(meta), path("medaka/*/*_medaka_consensus.fasta"), emit: fasta
  path "medaka/*/*", emit: everything
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p medaka

  medaka_consensus ${args} \
    -i ${fastq} \
    -d ${fasta} \
    -o medaka/${prefix} \
    -t ${task.cpus}

  if [ -f "medaka/${prefix}/consensus.fasta" ]; then cp medaka/${prefix}/consensus.fasta medaka/${prefix}/${prefix}_medaka_consensus.fasta ; fi

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    medaka: \$( medaka --version)
  END_VERSIONS

  exit 1
  """
}

process multiqc {
  tag           "combining reports"
  label         "process_low"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/multiqc:1.19'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'

  // TODO : check supported files
  //fastp
  //filtlong
  //busco

  input:
  file(input)

  output:
  path "multiqc/multiqc_report.html", emit: report
  path "multiqc/multiqc_data/*", emit: everything
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  """
  multiqc ${args} \
    --outdir multiqc \
    .

  exit 1
  """
}

process nanoplot_summary {
  tag           "${summary}"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/nanoplot:1.42.0'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'
  publishDir "${params.outdir}", mode: 'copy'

  input:
  file(summary)

  output:
  path "nanoplot/summary", emit: final_directory
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  """
  mkdir -p nanoplot/summary
  
  NanoPlot ${args} \
    --summary ${summary} \
    --threads ${task.cpus} \
    --outdir nanoplot/summary \
    --tsv_stats

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    nanoplot: \$(NanoPlot --version 2>&1)
  END_VERSIONS

  exit 1
  """
}

process nanoplot {
  tag           "${meta.id}"
  label         "process_low"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/nanoplot:1.42.0'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'

  input:
  tuple val(meta), file(fastq)

  output:
  path "nanoplot/*/*", emit: everything
  path "nanoplot/*/*_NanoStats.csv",emit: summary
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p nanoplot/${prefix}

  NanoPlot !{params.nanoplot_options} \
    --fastq !{fastq} \
    --threads ${task.cpus} \
    --tsv_stats \
    --outdir nanoplot/${prefix}

  cp nanoplot/${prefix}/NanoStats.txt nanoplot/${prefix}/${prefix}_NanoStats.txt

  echo "sample,\$(cut -f 1 nanoplot/${prefix}/${prefix}_NanoStats.txt | tr '\\n' ',' )" >  nanoplot/${prefix}/${prefix}_NanoStats.csv
  echo "${prefix},\$(cut -f 2 nanoplot/${prefix}/${prefix}_NanoStats.txt | tr '\\n' ',' )" >> nanoplot/${prefix}/${prefix}_NanoStats.csv

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    nanoplot: \$(NanoPlot --version 2>&1)
  END_VERSIONS

  exit 1

  """
}

process ontime {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/ontime:0.2.3'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '45m'
  
  input:
  tuple val(meta), file(reads)

  output:
  tuple val(meta), file("bbduk/*_rmphix_R{1,2}.fastq.gz"),  emit: fastq
  path "bbduk/*",                                           emit: files
  path "bbduk/*.phix.stats.txt",                            emit: stats
  path "logs/${task.process}/*.log",  emit: log
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  ontime --version

  ontime --help

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    ontime: "\$(ontime --version | awk '{print \$NF}')"
  END_VERSIONS

  exit 1
  """
}

process polypolish {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/polypolish:0.6.0'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '45m'

  input:
  tuple val(meta), file(fasta), file(sam)

  output:
  tuple val(meta), file("polypolish/${sample}_polypolish.fasta"), emit: fasta
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p polypolish

  polypolish -h

  polypolish -v

  polypolish filter \
    --in1 ${sam[0]} \
    --in2 ${sam[1]} \
    --out1 ${prefix}_filtered_1.sam \
    --out2 ${prefix}_filtered_2.sam

  polypolish ${args} \
    ${fasta} \
    ${prefix}_filtered_1.sam \
    ${prefix}_filtered_2.sam \
    > polypolish/${prefix}_polypolish.fasta

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    polypolish: \$(polypolish -v )
  END_VERSIONS

  exit 1
  """
}

process rasusa {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/rasusa:0.8.0'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'

  input:
  tuple val(meta), file(fastq)

  output:
  tuple val(meta), file("rasusa/${sample}/*.fastq.gz"), emit: fastq
  path "versions.yml", emit: versions 
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args       = task.ext.args   ?: ''
  def prefix     = task.ext.prefix ?: "${meta.id}"
  def replicates = task.ext.prefix ?: 1
  if (replicates == 1 ) {
    """
    mkdir -p rasusa/${prefix}

    rasusa ${args} \
      -i ${fastq} \
      -O g \
      --output rasusa/${prefix}/Whatever.txt

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      rasusa: \$(rasusa --version | sed -e "s/rasusa //g")
    END_VERSIONS

    exit 1
    """
  } else { 
    """
    mkdir -p rasusa/${prefix}

    for i in 0..${replicate}
    do
      rasusa ${args} \
        -i ${fastq} \
        -O g \
        --output rasusa/${prefix}/Whatever.${replicate}.txt
    done

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      rasusa: \$(rasusa --version | sed -e "s/rasusa //g")
    END_VERSIONS

    exit 1
    """
  }
}

process raven {
  tag           "${meta.id}"
  label         "process_high"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/raven:1.8.3'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10h'

  input:
  tuple val(meta), file(fastq)

  output:
  tuple val(meta), file("raven/*/*_raven.fasta"), emit: fasta
  tuple val(meta), file("raven/*/*_raven.gfa"), emit: gfa
  path("raven/*/*"), emit: everything
  path "versions.yml", emit: versions 
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p raven/${prefix}

  raven ${args} \
    --threads ${task.cpus} \
    --graphical-fragment-assembly raven/${prefix}/${prefix}_raven.gfa \
    ${fastq} \
    > raven/${prefix}/${prefix}_raven.fasta

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    raven: \$( raven --version )
  END_VERSIONS

  exit 1
  """
}

process summary {
  tag           "${meta.id}"
  label         "process_single"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/multiqc:1.19'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'
  
  input:
  file(input)

  output:
  path "summary/${sample}.summary.tsv",                          emit: summary
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p summary

  cat <<-END_VERSIONS > versions.yml
  "multiqc":
    multiqc: \$(multiqc --version | sed -e "s/multiqc, version //g" ))
  END_VERSIONS

  exit 1
  """
}

process unicycler {
  tag           "${meta.id}"
  label         "process_high"
  publishDir    params.outdir, mode: 'copy'
  container     'staphb/unicycler:0.5.0'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10h'

  input:
  tuple val(meta), file(nanopore), file(illumina)

  output:
  tuple val(meta), file("unicycler/*/*.fasta"), emit: fasta
  tuple val(meta), file("unicycler/*/*.gfa"), emit: gfa
  path "unicycler/*", emit: everything
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p unicycler

  unicycler ${args} \
    -1 ${illumina[0]} \
    -2 ${illumina[1]} \
    -l ${nanopore} \
    -o unicycler/${prefix} \
    -t ${task.cpus}

  if [ -f "unicycler/${prefix}/assembly.fasta" ] ; then cp unicycler/${prefix}/assembly.fasta unicycler/${prefix}/${prefix}.fasta ; fi
  if [ -f "unicycler/${prefix}/assembly.gfa" ] ; then cp unicycler/${prefix}/assembly.gfa unicycler/${prefix}/${prefix}.gfa ; fi

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    unicycler: \$(echo \$(unicycler --version 2>&1) | sed 's/^.*Unicycler v//; s/ .*\$//')
  END_VERSIONS

  exit 1
  """
}

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Subworkflow

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

workflow ASSEMBLY {
  take:
    ch_fastq_n
    ch_fastq_i

  main:
    ch_gfa       = Channel.empty()
    ch_summary   = Channel.empty()
    ch_versions  = Channel.empty()
    ch_consensus = Channel.empty()

    if (params.assembler =~ /unicycler/ ) {
      unicycler(ch_fastq_n.join(ch_fastq_i, by: 0 , remainder: false))
      
      ch_gfa       = ch_gfa.mix(unicycler.out.gfa)
      ch_consensus = ch_consensus.mix(unicycler.out.fasta)
      ch_versions  = ch_versions.mix(unicycler.out.versions.first())
    }
    
    if (params.assembler =~ /hybracter/ ) {
      hybracter(ch_fastq_n.join(ch_fastq_i, by: 0 , remainder: true))
      
      ch_gfa       = ch_gfa.mix(hybracter.out.gfa)
      ch_versions  = ch_versions.mix(hybracter.out.versions.first())
    }    
 
    if (params.assembler =~ /raven/ ) {
      raven(ch_fastq_n)
      
      ch_gfa      = ch_gfa.mix(raven.out.gfa)
      ch_versions = ch_version.mix(raven.out.versions.first())
    }
    
    if (params.assembler =~ /flye/ ) {
      flye(ch_fastq_n)

      flye.out.summary
        .collectFile(
          storeDir: "${params.outdir}/flye/",
          keepHeader: true,
          sort: { file -> file.text },
          name: "flye_summary.tsv")
        .set { flye_summary }

      ch_summary = ch_summary.mix(flye_summary)
      ch_gfa      = ch_gfa.mix(flye.out.gfa)
      ch_versions = ch_versions.mix(flye.out.versions.first())
    }
    
    if (params.assembler =~ /dragonflye/ ) {
      dragonflye(ch_fastq_n)

      dragonflye.out.summary
        .collectFile(
          storeDir: "${params.outdir}/dragonflye/",
          keepHeader: true,
          sort: { file -> file.text },
          name: "dragonflye_summary.tsv")
        .set { dragonflye_summary }

      ch_gfa       = dragonflye.out.gfa
      ch_versions  = ch_versions.mix(dragonflye.out.versions.first())
      ch_summary   = ch_summary.mix(dragonflye_summary)
      ch_consensus = ch_consensus.mix(dragonflye.out.fasta)
    }

    bandage(ch_gfa)
    gfastats(ch_gfa)
    // TODO : filter out unicycler and dragonflye
    gfa_to_fasta(ch_gfa)
    dnaapler(gfa_to_fasta.out.fasta.filter())

    gfastats.out.summary
      .collectFile(
        storeDir: "${params.outdir}/gfastats/",
        keepHeader: true,
        sort: { file -> file.text },
        name: "gfastats_summary.csv")
      .set { gfastats_summary }

    emit:
      consensus = dnaapler.out.fasta.mix(ch_consensus)
      summary   = ch_summary.mix(bandage.out.summary).mix(gfastats_summary).mix(dnaapler_summary)
      versions  = ch_versions.mix(bandage.out.versions.first()).mix(gfastats.out.versions.first()).mix(dnaapler.out.versions.first())
}

workflow FILTER {
  take:
    ch_fastq_n
    ch_fastq_i

  main:
    fastp(ch_fastq_i)

    // TODO : filter out fastp reads when there aren't enough

    if ( params.ontime ) {
      ontime(ch_fastq_n)
      ch_nanopore = ontime.out.fastq
    } else {
      ch_nanopore = ch_fastq_n
    }


    filtlong(ch_nanopore.join(fastp.out.reads, by: 0, remainder: true))

    // TODO : filter out filtlong reads when there aren't enough

    rasusa(filtlong.out.fastq)

  emit:
    fastq_n = rasusa.out.fastq.transpose()
    fastq_i = fastp.out.reads
    summary = fastp.out.summary
    version = fastp.out.versions.first().mix(filtlong.out.versions.first()).mix(rasusa.out.versions.first())
}

workflow METRICS {
  take:
    ch_fastq_n
    ch_consensus
    ch_summary
    ch_fastq_i
    ch_nanoplot_summary

  main:
    nanoplot_summary(ch_nanoplot_summary)
    nanoplot(ch_reads)
    busco(ch_consensus)

    circulocov(ch_consensus.join(ch_fastq_n, by: 0 , remainder: false).join(ch_fastq_i, by: 0, remainder: true))

    circulocov.out.summary
      .collectFile(name: "circulocov_summary.txt",
        keepHeader: true,
        storeDir: "${params.outdir}/circulocov")
      .set {circulocov_summary }

    nanoplot.out.summary
      .collectFile(name: "NanoStats.csv",
        keepHeader: true,
        storeDir: "${params.outdir}/nanoplot")
      .set { nanostats_summary }

    multiqc(ch_summary.mix(busco.out.summary).mix(nanostats_summary).mix(nanoplot_summary.out).mix(circulocov_summary).collect())
}

workflow POLISH {
  take:
    ch_assembly
    ch_fastq_n
    ch_fastq_i

  main:
    medaka(ch_fasta.join(ch_fastq_n, by:0, remainder: false))
    bwa(medaka.out.fasta.join(ch_fastq_i, by:0, remainder: false))
    polypolish(bwa.out.sam)
    pypolca(polypolish.out.fasta.join(ch_fastq_i, by:0, reads: false))
    
    medaka.out.fasta
      .mix(polypolish.out.fasta)
      .mix(pypolca.out.fasta)
      .set{ch_consensus}

    // add a process to compress fasta files for download from nf-tower?

  emit:
    fasta = ch_consensus
}

workflow DONUT_FALLS {
  take:
    ch_fastq_n
    ch_fastq_i
    ch_nanoplot_summary

  main:

    // filter out low quality and short reads to ~100X depth
    FILTER(ch_fastq_n, ch_fastq_i)

    // assemble with desired assemblers
    // ASSEMBLY(FILTER.out.fastq_n, FILTER.out.fastq_i)

    // // long-read only and short-read polishing
    // POLISH(ASSEMBLY.out.consensus, FILTER.out.fastq_n, FILTER.out.fastq_i)

    // // summary files and graphs for everything
    // METRICS(
    //   FILTER.out.fastq_n,
    //   ASSEMBLY.out.consensus.mix(POLISH.out.consensus), 
    //   FILTER.out.summary.mix(ASSEMBLY.out.summary).mix(POLISH.out.summary),
    //   FILTER.out.fastq_i,
    //   ch_nanoplot_summary)
}

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Workflow

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####


workflow {
  DONUT_FALLS(
    ch_fastq_n,
    ch_fastq_i.ifEmpty([]),
    ch_nanoplot_summary.ifEmpty([])
    )
}

workflow.onComplete {
  println("Pipeline completed at: $workflow.complete")
  println("The multiqc report can be found at ${params.outdir}/multiqc/multiqc_report.html")
  println("The consensus fasta files can be found in ${params.outdir}/consensus")
  println("The fasta files are from each phase of assembly. polca > polypolish > medaka > unpolished")
  println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}