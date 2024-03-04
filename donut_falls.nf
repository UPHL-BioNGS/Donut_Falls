#!/usr/bin/env nextflow
nextflow.enable.dsl = 2

// read but ignored most things from
// https://carpentries-incubator.github.io/Pipeline_Training_with_Nextflow/07-Nextflow_Best_Practice/index.html

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Greetings!

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

println('')
println(' __                    ___          ')
println('|  ) _   _      _)_    )_ _   ) ) _ ')
println('|_/ (_) ) ) (_( (_    (  (_( ( ( (  ')
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
params.test               = ''

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
    "config_file"]

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

if (params.sequencing_summary){
  Channel
    .fromPath("${params.sequencing_summary}", type:'file')
    .view { "Summary File : $it" }
    .set { ch_sequencing_summary }
} else {
  ch_sequencing_summary = Channel.empty()
}


// using a sample sheet with the column header of 'sample,fastq,fastq_1,fastq_2'
// sample  = meta.id
// fastq   = nanopore fastq file
// fastq_1 = illumina fastq file
// fastq_2 = illumina fastq file


if (params.sample_sheet) {
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
    .set{ ch_input_files }
} else {
  ch_input_files = Channel.empty()
}


// channel for illumina files (paired-end only)
ch_input_files
  .filter { it[2] != it[3] }
  .map { it -> tuple (it[0], [file(it[2], checkIfExists: true), file(it[3], checkIfExists: true)])}
  .set { ch_illumina_input }

// channel for nanopore files
ch_input_files
  .map { it -> tuple (it[0], file(it[1], checkIfExists: true))}
  .set { ch_nanopore_input }

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// TODO

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// process ontime {
//   tag           "${meta.id}"
//   label         "process_medium"
//   publishDir    "${params.outdir}/${meta.id}", mode: 'copy'
//   container     'staphb/ontime:0.2.3'
//   errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
//   time          '45m'
//
//   input:
//   tuple val(meta), file(reads)
//
//   output:
//   tuple val(meta), file("bbduk/*_rmphix_R{1,2}.fastq.gz"),  emit: fastq
//   path "bbduk/*",                                           emit: files
//   path "bbduk/*.phix.stats.txt",                            emit: stats
//   path "logs/${task.process}/*.log",  emit: log
//   path "versions.yml", emit: versions
//
//   when:
//   task.ext.when == null || task.ext.when
//
//   shell:
//   def args   = task.ext.args   ?: ''
//   def prefix = task.ext.prefix ?: "${meta.id}"
//   """
//   ontime --version
//
//   ontime --help
//
//   cat <<-END_VERSIONS > versions.yml
//   "${task.process}":
//     ontime: "\$(ontime --version | awk '{print \$NF}')"
//   END_VERSIONS
//
//   exit 1
//   """
// }


// someday...
// process dragonflye {
//   tag           "${meta.id}"
//   label         "process_high"
//   publishDir    "${params.outdir}/${meta.id}", mode: 'copy'
//   container     'staphb/dragonflye:1.1.2'
//   errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
//   time          '10h'
//
//   input:
//   tuple val(meta), file(fastq)
//
//   output:
//   tuple val(meta), file("dragonflye/*_dragonflye.fasta"), optional: true, emit: fasta
//   tuple val(meta), file("dragonflye/*_dragonflye.gfa"),   optional: true, emit: gfa
//   path "dragonflye/*_assembly_info.tsv", emit: summary
//   path "dragonflye/*", emit: everything
//   path "versions.yml", emit: versions
//
//   when:
//   task.ext.when == null || task.ext.when
//
//   shell:
//   def args   = task.ext.args   ?: ''
//   def prefix = task.ext.prefix ?: "${meta.id}"
//   """
//   dragonflye ${args} \
//     --reads ${fastq} \
//     --cpus ${task.cpus} \
//     --outdir dragonflye \
//     --prefix ${prefix}
//
//   # renaming final files
//   if [ -f "dragonflye/flye-unpolished.gfa" ] ; then cp dragonflye/flye-unpolished.gfa dragonflye/${prefix}_dragonflye.gfa   ; fi
//   if [ -f "dragonflye/flye.fasta" ]          ; then cp dragonflye/flye.fasta          dragonflye/${prefix}_dragonflye.fasta ; fi
//
//   # getting a summary file
//   head -n 1 dragonflye/flye-info.txt | awk '{print "sample\\t" \$0}' > dragonflye/${prefix}_assembly_info.tsv
//   tail -n+2 dragonflye/flye-info.txt | awk -v sample=${prefix} '{print sample "\\t" \$0}' >> dragonflye/${prefix}_assembly_info.tsv
//   
//   cat <<-END_VERSIONS > versions.yml
//   "${task.process}":
//     dragonflye: \$(dragonflye --version | awk '{print \$NF}' )
//   END_VERSIONS
//   """
// }

// someday...
// process hybracter {
//   tag           "${meta.id}"
//   label         "process_high"
//   publishDir    "${params.outdir}/${meta.id}", mode: 'copy'
//   container     'quay.io/biocontainers/hybracter:0.6.0--pyhdfd78af_0'
//   errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
//   time          '10h'
//
//   input:
//   tuple val(meta), file(reads), file(illumina)
//
//   output:
//   tuple val(meta), file("bbduk/*_rmphix_R{1,2}.fastq.gz"), emit: fasta
//   tuple val(meta), file("bbduk/*_rmphix_R{1,2}.fastq.gz"), emit: gfa
//   path "versions.yml", emit: versions
//
//   when:
//   task.ext.when == null || task.ext.when
//
//   shell:
//   def args   = task.ext.args   ?: ''
//   def prefix = task.ext.prefix ?: "${meta.id}"
//   """
//   hybracter -h
//
//   hybracter version
// 
//   exit 1
//
//   cat <<-END_VERSIONS > versions.yml
//   "${task.process}":
//     hybracter: "\$(hybracter --version | awk '{print \$NF}')"
//   END_VERSIONS
//   exit 1
//   """
// }

// process test_nfcore {
//   tag           "Downloading subset15000"
//   label         "process_single"
//   publishDir    "${params.outdir}/test_files/nfcore", mode: 'copy'
//   container     'staphb/multiqc:1.19'
//   errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
//   time          '1h'

//   output:
//   tuple val("nfcore-subset15000"), file("nfcore_subset15000.fa.gz"), emit: fastq

//   when:
//   task.ext.when == null || task.ext.when

//   shell:
//   """
//   wget -q https://github.com/nf-core/test-datasets/blob/23f5b889e4736798c8692e9b92810d9a3e37ee97/nanopore/subset15000.fq.gz?raw=true -O nfcore_subset15000.fa.gz
//   """
// }

// process test_great_dataset {
//   tag           "Downloading the great dataset"
//   label         "process_single"
//   publishDir    "${params.outdir}/test_files/great", mode: 'copy'
//   container     'staphb/multiqc:1.19'
//   errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
//   time          '1h'

//   output:
//   tuple val("great_dataset"), file("reads.fastq.gz"), emit: fastq

//   when:
//   task.ext.when == null || task.ext.when

//   shell:
//   """
//   wget -q https://bridges.monash.edu/ndownloader/files/23754659 -O dataset.tar.gz
//   tar -xvf dataset.tar.gz

//   exit 1
//   """
// }


// process test_good_dataset {
//   tag           "Downloading the good dataset"
//   label         "process_single"
//   publishDir    "${params.outdir}/test_files/good", mode: 'copy'
//   container     'staphb/multiqc:1.19'
//   errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
//   time          '1h'

//   output:
//   tuple val("good_dataset"), file("reads.fastq.gz"), emit: fastq

//   when:
//   task.ext.when == null || task.ext.when

//   shell:
//   """
//   wget -q https://bridges.monash.edu/ndownloader/files/23754647 -O dataset.tar.gz
//   tar -xvf dataset.tar.gz
//   """
// }

// process test_mediocre_dataset {
//   tag           "Downloading the mediocre dataset"
//   label         "process_single"
//   publishDir    "${params.outdir}/test_files/mediocre", mode: 'copy'
//   container     'staphb/multiqc:1.19'
//   errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
//   time          '1h'

//   output:
//   tuple val("mediocre_dataset"), file("reads.fastq.gz"), emit: fastq

//   when:
//   task.ext.when == null || task.ext.when

//   shell:
//   """
//   wget -q https://bridges.monash.edu/ndownloader/files/23754629 -O dataset.tar.gz
//   tar -xvf dataset.tar.gz

//   exit 1
//   """
// }

// process test_bad_dataset {
//   tag           "Downloading the bad dataset"
//   label         "process_single"
//   publishDir    "${params.outdir}/test_files/bad", mode: 'copy'
//   container     'staphb/multiqc:1.19'
//   errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
//   time          '1h'

//   output:
//   tuple val("bad_dataset"), file("reads.fastq.gz"), emit: fastq

//   when:
//   task.ext.when == null || task.ext.when

//   shell:
//   """
//   wget -q https://bridges.monash.edu/ndownloader/files/23754623 -O dataset.tar.gz
//   tar -xvf dataset.tar.gz

//   exit 1
//   """
// }

    // in DONUT FALLS WORKFLOW
    // hybracter and plassembler are on the to-do list
    // if (params.assembler =~ /hybracter/ ) {
    //   hybracter(ch_nanopore_input.join(ch_illumina_input, by: 0 , remainder: true))
    //      
    //   ch_gfa       = ch_gfa.mix(hybracter.out.gfa)
    //   // no ch_summary
    //   ch_consensus = ch_consensus.mix(hybracter.out.fasta)
    //   ch_versions  = ch_versions.mix(hybracter.out.versions.first())
    // } 

    // if (params.assembler =~ /dragonflye/ ) {
    //   dragonflye(ch_nanopore_input)
    //
    //   dragonflye.out.summary
    //     .collectFile(
    //       storeDir: "${params.outdir}/summary/",
    //       keepHeader: true,
    //       sort: { file -> file.text },
    //       name: "dragonflye_summary.tsv")
    //     .set { dragonflye_summary }
    //
    //   ch_gfa       = dragonflye.out.gfa
    //   ch_summary   = ch_summary.mix(dragonflye_summary)
    //   // no ch_consensus
    //   ch_versions  = ch_versions.mix(dragonflye.out.versions.first())
    // }

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Processes

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

process bandage {
  tag           "${meta.id}"
  label         "process_low"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy'
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
  def prefix = task.ext.prefix ?: "${gfa.baseName}"
  """
  mkdir -p bandage

  Bandage image ${gfa} bandage/${prefix}.png ${args}
  Bandage image ${gfa} bandage/${prefix}.svg ${args}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    bandage: \$(Bandage --version | awk '{print \$NF }')
  END_VERSIONS
  """
}

process busco {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy'
  container     'staphb/busco:5.6.1-prok-bacteria_odb10_2024-01-08'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '45m'

  input:
  tuple val(meta), file(fasta)

  output:
  path("busco/*/*"), emit: everything
  path("busco/*/short_summary*.txt"), optional: true, emit: summary
  path "versions.yml" , emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: '--offline -l /busco_downloads/lineages/bacteria_odb10'
  def prefix = task.ext.prefix ?: "${fasta.baseName}"
  """
  busco ${args} \
    -m genome \
    -i ${fasta} \
    -o busco/${prefix} \
    --cpu ${task.cpus}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    busco: \$( busco --version | awk '{print \$NF}' )
  END_VERSIONS
  """
}

process bwa {
  tag           "${meta.id}"
  label         'process_high'
  // no publishDir because the sam files are too big
  container     'staphb/bwa:0.7.17'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '2h'

  input:
  tuple val(meta), file(fasta), file(fastq)

  output:
  tuple val(meta), file(fasta), file("bwa/*_{1,2}.sam"), emit: sam
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${fasta.baseName}"
  """
  mkdir -p bwa

  bwa index ${fasta}
  bwa mem -t ${task.cpus} -a ${fasta} ${fastq[0]} > bwa/${prefix}_1.sam
  bwa mem -t ${task.cpus} -a ${fasta} ${fastq[1]} > bwa/${prefix}_2.sam

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    bwa: \$(bwa 2>&1 | grep -i version | awk '{print \$NF}')
  END_VERSIONS
  """
}

process circulocov {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy'
  container     'quay.io/uphl/circulocov:0.1.20240104-2024-02-21'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '1h'

  input:
  tuple val(meta), file(fasta), file(nanopore), file(illumina)

  output:
  path "circulocov/*overall_summary.txt", emit: summary
  tuple val(meta), file("circulocov/*/overall_summary.txt"), emit: results
  path "circulocov/*/*", emit: everything
  path "circulocov/*/fastq/*", emit: fastq
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: '-a'
  def prefix = task.ext.prefix ?: "${fasta.baseName}"
  def reads  = (illumina =~ /input/) ? "" : "--illumina ${illumina.join(' ')}"
  """
  mkdir -p circulocov/${prefix}

  circulocov ${args} \
    --threads ${task.cpus} \
    --genome ${fasta} \
    --nanopore ${nanopore} \
    ${reads} \
    --out circulocov/${prefix} \
    --sample ${prefix}
  
  cp circulocov/${prefix}/overall_summary.txt circulocov/${prefix}_overall_summary.txt

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    circulocov: \$(circulocov -v | awk '{print \$NF}')
  END_VERSIONS
  """
}

process copy {
  tag           "${meta.id}"
  label         "process_single"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy'
  container     'staphb/multiqc:1.19'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'

  input:
  tuple val(meta), file(fasta), file(circulocov), file(gfastats)

  output:
  path "consensus/*", emit: fastas
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  #!/usr/bin/env python3
  import glob
  import json
  import csv
  import os

  def gfastats_to_dict(header_dict):
    dict = {}
    with open("gfastats_summary.csv", mode='r', newline='') as file:
      reader = csv.DictReader(file)
      for row in reader:
        if row["sample"] == header_dict['name'] + "_" + header_dict['assembler']:
          key = row["Header"]

          dict[key] = row
    return dict

  def circulocov_to_dict(header_dict):
    dict = {}
    with open("circulocov_summary.txt", mode='r', newline='') as file:
        reader = csv.DictReader(file, delimiter="\t")
        for row in reader:
          if row["sample"].replace("_reoriented","") == header_dict['name'] + "_" + header_dict['assembler'] :
            key = row["contigs"]

            dict[key] = row
    return dict

  def copy_fasta(fasta, header_dict, gfa_dict, circulocov_dict):
    with open(fasta, 'r') as file:
      with open(f"consensus/{header_dict['fasta']}", 'w') as outfile:
          for line in file:
            line = line.strip()
            if line.startswith('>'):
              contig = line.replace(">","").split()[0]
              circular = gfa_dict[contig]['circular'].replace("Y","true").replace("N","false")
              length = gfa_dict[contig]['Total segment length']
              gc_per = gfa_dict[contig]['GC content %']
              meandepth = circulocov_dict[contig]['nanopore_meandepth']
              assembler = header_dict['assembler']
              step = header_dict['step']
              outfile.write(f">{contig} circ={circular} len={length} gc={gc_per} cov={meandepth} asmb={assembler} stp={step}\\n")
            else:
              outfile.write(f"{line}\\n")

  def main():
    os.mkdir("consensus")
    header_dict = {}
    fasta = glob.glob("*.fasta")[0]
    header_dict['fasta'] = fasta

    name = fasta.replace(".fasta", "")

    assemblers = ['dragonflye', 'flye', 'hybracter', 'raven', 'unicycler']
    steps = ["reoriented", 'polypolish', 'pypolca', 'medaka']
    for step in steps:
      if step in name:
        header_dict['step'] = step
        name = name.replace(f"_{step}","")
        break

    if 'step' not in header_dict.keys():
      header_dict['step'] = False

    for assembler in assemblers:
      if assembler in name:
        header_dict['assembler'] = assembler
        name = name.replace(f"_{assembler}","")
        break

    header_dict['name'] = name

    gfa_dict        = gfastats_to_dict(header_dict)
    circulocov_dict = circulocov_to_dict(header_dict)

    copy_fasta(fasta, header_dict, gfa_dict, circulocov_dict)

  if __name__ == "__main__":
    main()
  """
}

process dnaapler {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy'
  container     'staphb/dnaapler:0.7.0'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '1h'

  input:
  tuple val(meta), file(fasta), file(ignore)

  output:
  tuple val(meta), file("dnaapler/*_reoriented.fasta"), emit: fasta
  path "dnaapler/*", emit: files
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${fasta.baseName}"
  """
  dnaapler all ${args} \
    --input ${fasta} \
    --prefix ${prefix} \
    --output dnaapler \
    --threads ${task.cpus} \
    --ignore ${ignore}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    dnaapler: \$(dnaapler --version | awk '{print \$NF}')
  END_VERSIONS
  """
}

process fastp {
  tag           "${meta.id}"
  label         "process_low"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy'
  container     'staphb/fastp:0.23.4'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'

  input:
  tuple val(meta), file(reads), val(type)

  output:
  tuple val(meta), file("fastp/*_fastp*.fastq.gz"), val(type), emit: fastq
  path "fastp/*", emit: everything
  path "fastp/*_fastp*.json", emit: summary
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def lrargs = task.ext.lrargs ?: '--qualified_quality_phred 12 --length_required 1000'
  def prefix = task.ext.prefix ?: "${meta.id}"
  if (type == "illumina"){
    """
    mkdir -p fastp

    fastp ${args} \
      --in1 ${reads[0]} \
      --in2 ${reads[1]} \
      --out1 fastp/${prefix}_fastp_sr_R1.fastq.gz \
      --out2 fastp/${prefix}_fastp_sr_R2.fastq.gz \
      -h fastp/${prefix}_fastp_sr.html \
      -j fastp/${prefix}_fastp_sr.json

    passed_filter_reads=\$(grep passed_filter_reads fastp/${prefix}_fastp_sr.json | awk '{print \$NF}' | head -n 1 )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      fastp: \$(fastp --version 2>&1 | awk '{print \$NF}' )
    END_VERSIONS
    """
  } else {
    """
    mkdir -p fastp

    fastp ${lrargs} \
      --in1 ${reads[0]} \
      --out1 fastp/${prefix}_fastp_lr.fastq.gz \
      -h fastp/${prefix}_fastp_lr.html \
      -j fastp/${prefix}_fastp_lr.json

    passed_filter_reads=\$(grep passed_filter_reads fastp/${prefix}_fastp_sr.json | awk '{print \$NF}' | head -n 1 )

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
      fastp: \$(fastp --version 2>&1 | awk '{print \$NF}')
    END_VERSIONS
    """
  }
}

process flye {
  tag           "${meta.id}"
  label         "process_high"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy'
  container     'staphb/flye:2.9.3'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10h'

  input:
  tuple val(meta), file(fastq)

  output:
  tuple val(meta), file("flye/*_flye.fasta"), emit: fasta, optional: true
  tuple val(meta), file("flye/*_flye.gfa"), emit: gfa, optional: true
  path "flye/*_assembly_info.tsv", emit: summary
  path "flye/*", emit: everything
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p flye

  flye ${args} \
    --nano-raw ${fastq} \
    --threads ${task.cpus} \
    --out-dir flye

  # renaming final files
  if [ -f "flye/assembly.fasta" ]     ; then cp flye/assembly.fasta     flye/${prefix}_flye.fasta ; fi
  if [ -f "flye/assembly_graph.gfa" ] ; then cp flye/assembly_graph.gfa flye/${prefix}_flye.gfa   ; fi

  # getting a summary file
  head -n 1 flye/assembly_info.txt | awk '{print "sample\\t" \$0}' > flye/${prefix}_assembly_info.tsv
  tail -n+2 flye/assembly_info.txt | awk -v sample=${prefix} '{print sample "\\t" \$0}' >> flye/${prefix}_assembly_info.tsv

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    flye: \$( flye --version | awk '{print \$NF}')
  END_VERSIONS
  """
}

process gfastats {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', pattern: 'gfastats/*'
  container     'staphb/gfastats:1.3.6'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'

  input:
  tuple val(meta), file(gfa)

  output:
  tuple val(meta), file(gfa), file("gfastats/*_gfastats_summary.csv"), emit: stats
  path "gfastats/*_gfastats_summary.csv", emit: summary
  path "gfastats/*", emit: everything
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${gfa.baseName}"
  """
  mkdir -p gfastats

  gfastats \
    ${gfa} \
    ${args} \
    --threads ${task.cpus} \
    --tabular \
    --seq-report \
    > gfastats/${prefix}_gfastats.txt

  head -n 1 gfastats/${prefix}_gfastats.txt | tr "\\t" "," | awk '{print "sample," \$0 "circular" }' > gfastats/${prefix}_gfastats_summary.csv
  tail -n+2 gfastats/${prefix}_gfastats.txt | tr "\\t" "," | awk -v sample=${prefix} '{print sample "," \$0 }' >> gfastats/${prefix}_gfastats_summary.csv
  
  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    gfastats: \$( gfastats -v | head -n 1 | awk '{print \$NF}')
  END_VERSIONS
  """
}

process gfa_to_fasta {
  tag           "${meta.id}"
  label         "process_low"
  // no publishDir
  container     'staphb/multiqc:1.19'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'

  input:
  tuple val(meta), file(gfa), file(stats)

  output:
  tuple val(meta), file("*fasta"), file("noncircular.txt"), emit: fasta
  
  when:
  task.ext.when == null || task.ext.when

  """
  #!/usr/bin/env python3

  import csv
  import glob

  def convert_to_fasta(summary_dict, gfa_file):
      outfile = '_'.join(gfa_file.split('.')[:-1]) + ".fasta"
      with open(gfa_file, mode='r') as file:
          for line in file:
              parts = line.split()
              if parts and parts[0] == "S":
                  header = parts[1]
                  seq = parts[2]
                  if header in summary_dict.keys():
                      new_header = ">" + header + " length=" + summary_dict[header]['Total segment length'] + " circular=" + summary_dict[header]["circular"].replace("N","false").replace("Y","true") + " gc_per=" + summary_dict[header]["GC content %"] + "\\n"
                      with open(outfile, mode='a') as output_file:
                          output_file.write(new_header)
                          output_file.write(seq + "\\n")

  def read_summary_csv(gfastats_file):
      summary_dict = {}
      with open(gfastats_file, mode='r', newline='') as file:
          reader = csv.DictReader(file)
          for row in reader:
              key = row['Header']
              summary_dict[key] = row
              with open("noncircular.txt", mode='a') as output_file:
                  if summary_dict[key]["circular"] == "N":
                      output_file.write(key + "\\n")
      return summary_dict

  gfastats_file = glob.glob("*_gfastats_summary.csv")
  gfa_file = glob.glob("*.gfa")

  summary_dict = read_summary_csv(gfastats_file[0])
  convert_to_fasta(summary_dict, gfa_file[0])
  """
}

// From https://github.com/nanoporetech/medaka
// > It is not recommended to specify a value of --threads greater than 2 for medaka consensus since the compute scaling efficiency is poor beyond this.
// > Note also that medaka consensus may been seen to use resources equivalent to <threads> + 4 as an additional 4 threads are used for reading and preparing input data.
process medaka {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy'
  container     'ontresearch/medaka:v1.11.3'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '30m'

  input:
  tuple val(meta), path(fasta), path(fastq)

  output:
  tuple val(meta), path("medaka/*_medaka.fasta"), emit: fasta
  path "medaka/*", emit: everything
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${fasta.baseName.replaceAll('_reoriented','')}"
  """
  mkdir -p medaka

  # someday...
  # medaka tools resolve_model --auto_model consensus ${fastq}

  medaka_consensus ${args} \
    -i ${fastq} \
    -d ${fasta} \
    -o medaka \
    -t 1

  if [ -f "medaka/consensus.fasta" ]; then cp medaka/consensus.fasta medaka/${prefix}_medaka.fasta ; fi

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    medaka: \$( medaka --version | awk '{print \$NF}')
  END_VERSIONS
  """
}

process multiqc {
  tag           "combining reports"
  label         "process_low"
  publishDir    "${params.outdir}", mode: 'copy'
  container     'staphb/multiqc:1.19'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'

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
  if [ -f "pypolca_summary.tsv" ]
  then
    echo "# plot_type: 'table'" > pypolca_mqc.txt
    echo "# section_name: 'pypolca'" >> pypolca_mqc.txt
    echo "# description: 'Long read polishing'" >> pypolca_mqc.txt
    echo "# pconfig:" >> pypolca_mqc.txt
    echo "#     namespace: 'Cust Data'" >> pypolca_mqc.txt
    echo "# headers:" >> pypolca_mqc.txt
    echo "#     Substitution_Errors_Found:" >> pypolca_mqc.txt
    echo "#         title: 'Substitution Errors Found'" >> pypolca_mqc.txt
    echo "#         description: 'Substitution Errors Found'" >> pypolca_mqc.txt
    echo "#     Insertion/Deletion_Errors_Found:" >> pypolca_mqc.txt
    echo "#         title: 'Insertion/Deletion Errors Found'" >> pypolca_mqc.txt
    echo "#         description: 'Insertion/Deletion Errors Found'" >> pypolca_mqc.txt
    echo "#     Assembly_Size:" >> pypolca_mqc.txt
    echo "#         title: 'Assembly Size'" >> pypolca_mqc.txt
    echo "#         description: 'Assembly Size'" >> pypolca_mqc.txt
    echo "#     Consensus_Quality_Before_Polishing:" >> pypolca_mqc.txt
    echo "#         title: 'Consensus Quality Before Polishing'" >> pypolca_mqc.txt
    echo "#         description: 'Consensus Quality Before Polishing'" >> pypolca_mqc.txt
    echo "#     Consensus_QV_Before_Polishing:" >> pypolca_mqc.txt
    echo "#         title: 'Consensus QV Before Polishing'" >> pypolca_mqc.txt
    echo "#         description: 'Consensus QV Before Polishing'" >> pypolca_mqc.txt
    cat pypolca_summary.tsv >> pypolca_mqc.txt
  fi

  if [ -f "gfastats_summary.csv" ]
  then
    echo "# plot_type: 'table'" > gfastats_mqc.csv
    echo "# section_name: 'gfastats'" >> gfastats_mqc.csv
    echo "# description: 'Metrics for GFA files'" >> gfastats_mqc.csv
    echo "# pconfig:" >> gfastats_mqc.csv
    echo "#     namespace: 'Cust Data'" >> gfastats_mqc.csv
    echo "# headers:" >> gfastats_mqc.csv
    echo "#     sample:" >> gfastats_mqc.csv
    echo "#         title: 'Sample and analysis'" >> gfastats_mqc.csv
    echo "#         description: 'Sample and analysis that generated contig'" >> gfastats_mqc.csv
    echo "#     Header:" >> gfastats_mqc.csv
    echo "#         title: 'Header'" >> gfastats_mqc.csv
    echo "#         description: 'Name of contig'" >> gfastats_mqc.csv
    echo "#     Total segment length:" >> gfastats_mqc.csv
    echo "#         title: 'Total segment length'" >> gfastats_mqc.csv
    echo "#         description: 'Total segment length'" >> gfastats_mqc.csv
    echo "#     A:" >> gfastats_mqc.csv
    echo "#         title: 'A'" >> gfastats_mqc.csv
    echo "#         description: 'Number of A'" >> gfastats_mqc.csv
    echo "#     C:" >> gfastats_mqc.csv
    echo "#         title: 'C'" >> gfastats_mqc.csv
    echo "#         description: 'Number of C'" >> gfastats_mqc.csv
    echo "#     G:" >> gfastats_mqc.csv
    echo "#         title: 'G'" >> gfastats_mqc.csv
    echo "#         description: 'Number of G'" >> gfastats_mqc.csv
    echo "#     T:" >> gfastats_mqc.csv
    echo "#         title: 'T'" >> gfastats_mqc.csv
    echo "#         description: 'Number of T'" >> gfastats_mqc.csv
    echo "#     GC content %:" >> gfastats_mqc.csv
    echo "#         title: 'GC content %'" >> gfastats_mqc.csv
    echo "#         description: 'GC content %'" >> gfastats_mqc.csv
    echo "#     # soft-masked bases:" >> gfastats_mqc.csv
    echo "#         title: '# soft-masked bases'" >> gfastats_mqc.csv
    echo "#         description: '# soft-masked bases'" >> gfastats_mqc.csv
    echo "#     circular:" >> gfastats_mqc.csv
    echo "#         title: 'circular'" >> gfastats_mqc.csv
    echo "#         description: 'circular'" >> gfastats_mqc.csv
    cat gfastats_summary.csv | awk '{print NR ',' \$0}' >> gfastats_mqc.csv
  fi

  if [ -f "flye_summary.tsv" ]
  then
    echo "# plot_type: 'table'" > flye_mqc.csv
    echo "# section_name: 'flye'" >> flye_mqc.csv
    echo "# description: 'Assembly Info'" >> flye_mqc.csv
    echo "# pconfig:" >> flye_mqc.csv
    echo "#     namespace: 'Cust Data'" >> flye_mqc.csv
    echo "# headers:" >> flye_mqc.csv
    echo "#     sample:" >> flye_mqc.csv
    echo "#         title: 'Sample'" >> flye_mqc.csv
    echo "#         description: 'Sample that generated contig'" >> flye_mqc.csv
    echo "#     #seq_name:" >> flye_mqc.csv
    echo "#         title: '#seq_name'" >> flye_mqc.csv
    echo "#         description: 'Name of contig'" >> flye_mqc.csv
    echo "#     length:" >> flye_mqc.csv
    echo "#         title: 'length'" >> flye_mqc.csv
    echo "#         description: 'length'" >> flye_mqc.csv
    echo "#     cov.:" >> flye_mqc.csv
    echo "#         title: 'cov'" >> flye_mqc.csv
    echo "#         description: 'Coverage'" >> flye_mqc.csv
    echo "#     circ:" >> flye_mqc.csv
    echo "#         title: 'circ'" >> flye_mqc.csv
    echo "#         description: 'Whether contig is circular'" >> flye_mqc.csv
    echo "#     repeat:" >> flye_mqc.csv
    echo "#         title: 'repeat'" >> flye_mqc.csv
    echo "#         description: 'repeat'" >> flye_mqc.csv
    echo "#     mult.:" >> flye_mqc.csv
    echo "#         title: 'mult'" >> flye_mqc.csv
    echo "#         description: 'mult.'" >> flye_mqc.csv
    echo "#     alt_group:" >> flye_mqc.csv
    echo "#         title: 'alt_group'" >> flye_mqc.csv
    echo "#         description: 'alt_group'" >> flye_mqc.csv
    echo "#     graph_path:" >> flye_mqc.csv
    echo "#         title: 'graph_path'" >> flye_mqc.csv
    echo "#         description: 'graph_path'" >> flye_mqc.csv
    cat flye_summary.tsv | awk '{print NR '\\t' \$0}' >> flye_mqc.csv
  fi

  circulocov_check=\$(ls * | grep overall_summary.txt | head -n 1)
  if [ -n "\$circulocov_check" ]
  then
    illumina_check=\$(grep -h illumina *overall_summary.txt | head -n 1)
    if [ -n "\$illumina_check" ]
    then
      circulocov_summary_header=\$illumina_check
    else
      circulocov_summary_header=\$(grep -h nanopore_numreads *overall_summary.txt | head -n 1)
    fi

    echo \$circulocov_summary_header | awk '{print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" \$6 "\\t" \$7 "\\t" \$8 "\\t" \$9 "\\t" \$10 "\\t" \$11 "\\t" \$12}' > circulocov_summary.txt
    cat *overall_summary.txt | grep -v nanopore_numreads | awk '{print \$1 "\\t" \$2 "\\t" \$3 "\\t" \$4 "\\t" \$5 "\\t" \$6 "\\t" \$7 "\\t" \$8 "\\t" \$9 "\\t" \$10 "\\t" \$11 "\\t" \$12}' >> circulocov_summary.txt

    echo "# plot_type: 'table'" > circulocov_mqc.txt
    echo "# section_name: 'CirculoCov'" >> circulocov_mqc.txt
    echo "# description: 'Coverage estimates for circular sequences'" >> circulocov_mqc.txt
    echo "# pconfig:" >> circulocov_mqc.txt
    echo "#     namespace: 'Cust Data'" >> circulocov_mqc.txt
    echo "# headers:" >> circulocov_mqc.txt
    echo "#     sample:" >> circulocov_mqc.txt
    echo "#         title: 'Sample'" >> circulocov_mqc.txt
    echo "#         description: 'Sample that generated contig'" >> circulocov_mqc.txt
    echo "#     circ:" >> circulocov_mqc.txt
    echo "#         title: 'circ'" >> circulocov_mqc.txt
    echo "#         description: 'Whether contig was circular'" >> circulocov_mqc.txt
    echo "#     contigs:" >> circulocov_mqc.txt
    echo "#         title: 'contigs'" >> circulocov_mqc.txt
    echo "#         description: 'name of contig'" >> circulocov_mqc.txt
    echo "#     length:" >> circulocov_mqc.txt
    echo "#         title: 'length'" >> circulocov_mqc.txt
    echo "#         description: 'length of contig'" >> circulocov_mqc.txt
    echo "#     nanopore_numreads:" >> circulocov_mqc.txt
    echo "#         title: 'numreads'" >> circulocov_mqc.txt
    echo "#         description: 'number of nanopore reads mapping to contig'" >> circulocov_mqc.txt
    echo "#     nanopore_covbases:" >> circulocov_mqc.txt
    echo "#         title: 'covbases'" >> circulocov_mqc.txt
    echo "#         description: 'nanopore covbases of contig'" >> circulocov_mqc.txt
    echo "#     nanopore_coverage:" >> circulocov_mqc.txt
    echo "#         title: 'coverage'" >> circulocov_mqc.txt
    echo "#         description: 'nanopore coverage of contig'" >> circulocov_mqc.txt
    echo "#     nanopore_meandepth:" >> circulocov_mqc.txt
    echo "#         title: 'meandepth'" >> circulocov_mqc.txt
    echo "#         description: 'nanopore meandepth of contig'" >> circulocov_mqc.txt
    cat circulocov_summary.txt | awk '{print NR '\\t' \$0}' >> circulocov_mqc.txt
  fi

  touch whatever.png

  pngs=\$(ls *png)
  for png in \${pngs[@]}
  do
    new_name=\$(echo \$png | sed 's/.png\$/_mqc.png/g')
    cp \$png \$new_name
  done

  rm whatever_mqc.png

  multiqc ${args} \
    --outdir multiqc \
    .
  """
}

process nanoplot_summary {
  tag           "${summary}"
  label         "process_low"
  publishDir    "${params.outdir}/summary", mode: 'copy'
  container     'staphb/nanoplot:1.42.0'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'

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
  mkdir -p nanoplot
  
  NanoPlot ${args} \
    --summary ${summary} \
    --threads ${task.cpus} \
    --outdir nanoplot \
    --tsv_stats

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    nanoplot: \$(NanoPlot --version | awk '{print \$NF}'))
  END_VERSIONS

  exit 1
  """
}

process nanoplot {
  tag           "${meta.id}"
  label         "process_low"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy'
  container     'staphb/nanoplot:1.42.0'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'

  input:
  tuple val(meta), file(fastq)

  output:
  path "nanoplot/*", emit: everything
  path "nanoplot/${meta.id}_NanoStats.txt", emit: stats
  path "nanoplot/${meta.id}_NanoStats.csv", emit: summary
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p nanoplot

  NanoPlot ${args} \
    --fastq ${fastq} \
    --threads ${task.cpus} \
    --tsv_stats \
    --outdir nanoplot

  cp nanoplot/NanoStats.txt nanoplot/${prefix}_NanoStats.txt

  echo "sample,\$(   cut -f 1 nanoplot/${prefix}_NanoStats.txt | tr '\\n' ',' )" >  nanoplot/${prefix}_NanoStats.csv
  echo "${prefix},\$(cut -f 2 nanoplot/${prefix}_NanoStats.txt | tr '\\n' ',' )" >> nanoplot/${prefix}_NanoStats.csv

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    nanoplot: \$(NanoPlot --version | awk '{print \$NF}')
  END_VERSIONS
  """
}

process polypolish {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy'
  container     'staphb/polypolish:0.6.0'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '45m'

  input:
  tuple val(meta), file(fasta), file(sam)

  output:
  tuple val(meta), file("polypolish/*_polypolish.fasta"), emit: fasta
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def filarg = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${fasta.baseName.replaceAll('_medaka','')}"
  """
  mkdir -p polypolish

  polypolish filter \
    ${filarg} \
    --in1 ${sam[0]} \
    --in2 ${sam[1]} \
    --out1 ${prefix}_filtered_1.sam \
    --out2 ${prefix}_filtered_2.sam

  polypolish polish \
    ${args} \
    ${fasta} \
    ${prefix}_filtered_1.sam \
    ${prefix}_filtered_2.sam \
    > polypolish/${prefix}_polypolish.fasta

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    polypolish: \$(polypolish --version | awk '{print \$NF}')
  END_VERSIONS
  """
}

process pypolca {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy'
  container     'staphb/pypolca:0.3.1'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '30m'
  
  input:
  tuple val(meta), file(fasta), file(fastq)

  output:
  tuple val(meta), file("pypolca/*_pypolca.fasta"), optional: true, emit: fasta
  path "pypolca/*pypolca_summary.tsv", optional: true, emit: summary
  path "pypolca/*/*", emit: everything
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${fasta.baseName.replaceAll('_polypolish','')}"
  """
  pypolca run ${args}\
    -a ${fasta} \
    -1 ${fastq[0]} \
    -2 ${fastq[1]} \
    -t ${task.cpus} \
    -o pypolca/${prefix}

  if [ -f "pypolca/${prefix}/pypolca.report" ]
  then
    cut -f 1 -d : pypolca/${prefix}/pypolca.report | \
      sed 's/ /_/g' | \
      tr "\\n" "\\t" | \
      awk '{print "sample\\t" \$0 }' \
      > pypolca/${prefix}_pypolca_summary.tsv

    cut -f 2 -d : pypolca/${prefix}/pypolca.report | \
      awk '{( \$1 = \$1 ) ; print \$0 }' | \
      sed 's/ /_/g' | \
      tr "\\n" "\\t" | \
      awk '{print "${prefix}\\t" \$0 }' \
      >> pypolca/${prefix}_pypolca_summary.tsv
  fi  

  if [ -f "pypolca/${prefix}/pypolca_corrected.fasta" ]; then cp pypolca/${prefix}/pypolca_corrected.fasta pypolca/${prefix}_pypolca.fasta ; fi

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    pypolca: \$(pypolca --version | awk '{print \$NF}')
  END_VERSIONS
  """
}

process rasusa {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy'
  container     'staphb/rasusa:0.8.0'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'

  input:
  tuple val(meta), file(fastq)

  output:
  tuple val(meta), file("rasusa/*.fastq.gz"), emit: fastq
  path "versions.yml", emit: versions 
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args       = task.ext.args   ?: '--genome-size 5mb --coverage 150'
  def prefix     = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p rasusa

  rasusa ${args} \
    --input ${fastq} \
    --output rasusa/${prefix}_rasusa.fastq.gz

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    rasusa: \$(rasusa --version | awk '{print \$NF}' )
  END_VERSIONS
  """
}

process raven {
  tag           "${meta.id}"
  label         "process_high"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy'
  container     'staphb/raven:1.8.3'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10h'

  input:
  tuple val(meta), file(fastq)

  output:
  tuple val(meta), file("raven/*_raven.fasta"), emit: fasta, optional: true
  tuple val(meta), file("raven/*_raven.gfa"), emit: gfa, optional: true
  path("raven/*"), emit: everything
  path "versions.yml", emit: versions 
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: '--polishing-rounds 2'
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p raven

  raven ${args} \
    --threads ${task.cpus} \
    --graphical-fragment-assembly raven/${prefix}_raven.gfa \
    ${fastq} \
    > raven/${prefix}_raven.fasta

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    raven: \$( raven --version | awk '{print \$NF}' )
  END_VERSIONS
  """
}

process summary {
  tag           "Creating summary"
  label         "process_single"
  publishDir    "${params.outdir}/summary", mode: 'copy'
  container     'staphb/multiqc:1.19'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10m'
  
  input:
  file(input)

  output:
  path "donut_falls_summary.json", emit: summary
  
  when:
  task.ext.when == null || task.ext.when

  """
  #!/usr/bin/env python3

  import glob
  import json
  import csv
  from os.path import exists

  def file_to_dict(file, header, delim):
    dict = {}
    with open(file, mode='r', newline='') as file:
      reader = csv.DictReader(file, delimiter=delim)
      for row in reader:
        key = row[header]
        dict[key] = row
    return dict

  def file_to_dict_uniq(file, header, header2, delim):
    dict = {}
    with open(file, mode='r', newline='') as file:
      reader = csv.DictReader(file, delimiter=delim)
      for row in reader:
        if row[header] not in dict.keys():
          dict[row[header]] = {}
        key = row[header] + "_" + row[header2]
        dict[row[header]][key] = row
    return dict

  def final_file(dict):
    with open('donut_falls_summary.json', 'w') as json_file:
      json.dump(dict, json_file, indent=4)

  def main():
    if exists('nanoplot_summary.csv') :
      nanoplot_dict = file_to_dict('nanoplot_summary.csv', 'sample', ',')
    else:
      nanoplot_dict = {}

    if exists('pypolca_summary.tsv') :
      pypolca_dict  = file_to_dict('pypolca_summary.tsv', 'sample', '\t')
    else:
      pypolca_dict = {}

    if exists('gfastats_summary.csv') :
      gfastats_dict = file_to_dict_uniq('gfastats_summary.csv', 'sample', 'Header', ',')
    else:
      gfastats_dict = {}

    busco_dict = {}
    busco_files = glob.glob("short_summary*txt")
    for file in busco_files:
      sample_analysis = file.split(".")[-2]
      with open(file, 'r') as f:
        for line in f:
          if "C:" and "S:" and "D:" and "F:" and "M:" and "n:" in line:
            busco_dict[sample_analysis] = line.strip()
            break

    circulocov_dict = {}
    circulocov_files = glob.glob("*overall_summary.txt")
    for file in circulocov_files:
      sample_analysis = file.replace("_overall_summary.txt", "").replace("_reoriented", "")
      circulocov_dict[sample_analysis] = {}
      with open(file, 'r') as f:
        for line in f:
          parts = line.split()
          if parts[2] == "all":
            circulocov_dict[sample_analysis]["coverage"] = parts[7]

          if parts[2] == "missing":
            if len(parts) > 8:
              unmapped_illumina = parts[8]
            else:
              unmapped_illumina = 0
          
            circulocov_dict[sample_analysis]["unmapped_nanopore"] = parts[4]
            circulocov_dict[sample_analysis]["unmapped_illumina"] = unmapped_illumina
      
    final_results = {}
    assemblers = ['dragonflye', 'flye', 'hybracter', 'raven', 'unicycler']
    for key in nanoplot_dict.keys():
      final_results[key] = {}
      final_results[key]['name'] = key

      # from nanostas
      final_results[key]['number_of_reads'] = nanoplot_dict[key]['number_of_reads']
      final_results[key]['mean_read_length'] = nanoplot_dict[key]['mean_read_length']
      final_results[key]['mean_qual'] = nanoplot_dict[key]['mean_qual']
      for assembler in assemblers:
        if key + "_" + assembler in gfastats_dict.keys():
          final_results[key][assembler] = {}

          # gfastats results
          total_length  = 0
          num_circular = 0
          for contig in gfastats_dict[key + "_" + assembler].keys():
            total_length = total_length + int(gfastats_dict[key + "_" + assembler][contig]["Total segment length"])
            if gfastats_dict[key + "_" + assembler][contig]["circular"] == "Y":
              num_circular = num_circular + 1
          
          final_results[key][assembler]['total_length'] = total_length
          final_results[key][assembler]['num_contigs'] = len(gfastats_dict[key + "_" + assembler].keys())
          final_results[key][assembler]['circ_contigs'] = num_circular
                  
          # circulocov results
          if key + "_" + assembler in circulocov_dict.keys():
            final_results[key][assembler]['coverage'] = circulocov_dict[key + '_' + assembler]['coverage']
            final_results[key][assembler]['unmapped_nanopore'] = circulocov_dict[key + '_' + assembler]['unmapped_nanopore']
            final_results[key][assembler]['unmapped_illumina'] = circulocov_dict[key + '_' + assembler]['unmapped_illumina']

          # busco results
          if key + "_" + assembler in busco_dict.keys():
            final_results[key][assembler]['busco'] = busco_dict[key + "_" + assembler]
          if key + "_" + assembler + '_reoriented' in busco_dict.keys():                
            final_results[key][assembler]['busco'] = busco_dict[key + "_" + assembler + '_reoriented']
          for step in ['polypolish', 'pypolca', 'medaka']:
            if key + "_" + assembler + '_' + step in busco_dict.keys():                
              final_results[key][assembler]['busco_' + step ] = busco_dict[key + "_" + assembler + '_' + step]
            else:
              final_results[key][assembler]['busco_' + step ] = 'NF'

          # pypolca results
            if key + "_" + assembler in pypolca_dict.keys():
              final_results[key][assembler]['Consensus_Quality_Before_Polishing'] = pypolca_dict[key + "_" + assembler]['Consensus_Quality_Before_Polishing']
              final_results[key][assembler]['Consensus_QV_Before_Polishing'] = pypolca_dict[key + "_" + assembler]['Consensus_QV_Before_Polishing']
            else:
              final_results[key][assembler]['Consensus_Quality_Before_Polishing'] = 0
              final_results[key][assembler]['Consensus_QV_Before_Polishing'] = 0

      final_file(final_results)

  if __name__ == "__main__":
      main()

  """
}

process unicycler {
  tag           "${meta.id}"
  label         "process_high"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy'
  container     'staphb/unicycler:0.5.0'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '10h'

  input:
  tuple val(meta), file(illumina), file(nanopore)

  output:
  tuple val(meta), file("unicycler/*_unicycler.fasta"), emit: fasta, optional: true
  tuple val(meta), file("unicycler/*_unicycler.gfa"), emit: gfa, optional: true
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
    -o unicycler/ \
    -t ${task.cpus}

  if [ -f "unicycler/assembly.fasta" ] ; then cp unicycler/assembly.fasta unicycler/${prefix}_unicycler.fasta ; fi
  if [ -f "unicycler/assembly.gfa" ]   ; then cp unicycler/assembly.gfa   unicycler/${prefix}_unicycler.gfa   ; fi

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    unicycler: \$(unicycler --version | awk '{print \$NF }' )
  END_VERSIONS
  """
}

process versions {
  tag           "extracting versions"
  label         "process_single"
  publishDir    "${params.outdir}/summary", mode: 'copy'
  container     'staphb/multiqc:1.19'
  time          '10m'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}

  input:
  file(input)

  output:
  path "software_versions_mqc.yml", emit: versions
  path "software_versions.yml", emit: yml

  when:
  task.ext.when == null || task.ext.when

  """
  #!/usr/bin/env python3

  # Stolen and modified from 
  # https://github.com/nf-core/rnaseq/blob/b89fac32650aacc86fcda9ee77e00612a1d77066/modules/nf-core/custom/dumpsoftwareversions/templates/dumpsoftwareversions.py#L4

  import yaml
  from textwrap import dedent

  def _make_versions_html(versions):

      html = [
          dedent(
              \"\"\"\\
              <style>
              #nf-core-versions tbody:nth-child(even) {
                  background-color: #f2f2f2;
              }
              </style>
              <table class="table" style="width:100%" id="nf-core-versions">
                  <thead>
                      <tr>
                          <th> Process Name </th>
                          <th> Software </th>
                          <th> Version  </th>
                      </tr>
                  </thead>
              \"\"\"
          )
      ]
      for process, tmp_versions in sorted(versions.items()):
          html.append("<tbody>")
          for i, (tool, version) in enumerate(sorted(tmp_versions.items())):
              html.append(
                  dedent(
                      f\"\"\"\\
                      <tr>
                          <td><samp>{process if (i == 0) else ''}</samp></td>
                          <td><samp>{tool}</samp></td>
                          <td><samp>{version}</samp></td>
                      </tr>
                      \"\"\"
                  )
              )
          html.append("</tbody>")
      html.append("</table>")
      return "\\n".join(html)

  def main():

      with open("versions.yml") as f:
          versions_by_process = yaml.load(f, Loader=yaml.BaseLoader) 

      versions_by_module = {}
      for process, process_versions in versions_by_process.items():
          module = process.split(":")[-1]
          try:
              if versions_by_module[module] != process_versions:
                  raise AssertionError(
                      "There's something wrong with the designated containers of this workflow"
                  )
          except KeyError:
              versions_by_module[module] = process_versions

      versions_mqc = {
          "id": "software_versions",
          "section_name": "Donut Falls Software Versions",
          "section_href": "https://github.com/UPHL-BioNGS/Donut_Falls",
          "plot_type": "html",
          "description": "Collected at run time from the software output.",
          "data": _make_versions_html(versions_by_module),
      }

      with open("software_versions.yml", "w") as f:
          yaml.dump(versions_by_module, f, default_flow_style=False)

      with open("software_versions_mqc.yml", "w") as f:
          yaml.dump(versions_mqc, f, default_flow_style=False)

  if __name__ == "__main__":
      main()

  """
}

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Downloading files for testing

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

process test_unicycler {
  tag           "Downloading Unicycler test files"
  label         "process_single"
  publishDir    "${params.outdir}/test_files/unicycler", mode: 'copy'
  container     'staphb/multiqc:1.19'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '1h'

  output:
  tuple val("unicycler"), file("long_reads_low_depth.fastq.gz"), file("short_reads*.fastq.gz"), emit: fastq

  when:
  task.ext.when == null || task.ext.when

  shell:
  """
  wget --quiet https://github.com/rrwick/Unicycler/raw/69e712eb95c4b9f8a46aade467260260a9ce7a91/sample_data/short_reads_1.fastq.gz
  wget --quiet https://github.com/rrwick/Unicycler/raw/69e712eb95c4b9f8a46aade467260260a9ce7a91/sample_data/short_reads_2.fastq.gz
  wget --quiet https://github.com/rrwick/Unicycler/raw/69e712eb95c4b9f8a46aade467260260a9ce7a91/sample_data/long_reads_low_depth.fastq.gz
  """
}

process test_donut_falls {
  tag           "Downloading R10.4 reads"
  label         "process_single"
  publishDir    "${params.outdir}/test_files/df", mode: 'copy'
  container     'staphb/multiqc:1.19'
  errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
  time          '1h'

  output:
  tuple val("df"), file("test_files/test_nanopore.fastq.gz"), file("test_files/test_illumina_{1,2}.fastq.gz"), emit: fastq
  tuple val("df_lr"), file("test_files/test_nanopore_only.fastq.gz"), emit: lrfastq

  when:
  task.ext.when == null || task.ext.when

  shell:
  """
  wget --quiet https://zenodo.org/records/10779911/files/df_test_files.tar.gz?download=1 -O dataset.tar.gz
  tar -xvf dataset.tar.gz

  cp test_files/test_nanopore.fastq.gz test_files/test_nanopore_only.fastq.gz
  """
}

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Donut Falls

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

workflow DONUT_FALLS {
  take:
    ch_nanopore_input
    ch_illumina_input

  main:
    // channel for gfa files for gfa stats
    ch_gfa       = Channel.empty()
    // channel for files for multiqc or workflow summary
    ch_summary   = Channel.empty()
    // channel for assembled genomes
    ch_consensus = Channel.empty()
    // versions channel
    ch_versions  = Channel.empty()

    if (params.assembler =~ /unicycler/ ) {
      unicycler(ch_illumina_input.join(ch_nanopore_input, by: 0, remainder: false))
      
      ch_gfa       = ch_gfa.mix(unicycler.out.gfa)
      // no ch_summary
      ch_consensus = ch_consensus.mix(unicycler.out.fasta)
      ch_versions  = ch_versions.mix(unicycler.out.versions.first())
    }


    if (params.assembler.replaceAll('dragonflye','dragon') =~ /flye/ || params.assembler =~ /raven/ ) {
      // quality filter
      ch_illumina_input
        .map { it -> [it[0], it[1], "illumina"]}
        .mix(ch_nanopore_input.map { it -> [it[0], it[1], "nanopore"]})
        .filter{it[0]}
        .set { ch_input }

      fastp(ch_input)

      ch_versions = ch_versions.mix(fastp.out.versions)
      ch_summary  = ch_summary.mix(fastp.out.summary)

      fastp.out.fastq
        .filter { it[1].size() > 200 }
        .branch { it ->
          nanopore: it[2] == 'nanopore'
          illumina: it[2] == 'illumina'
        }
      .set { ch_filter }

      rasusa(ch_filter.nanopore.map {it -> tuple(it[0], it[1])})

      ch_versions = ch_versions.mix(rasusa.out.versions)

      if (params.assembler =~ /raven/ ) {
        raven(rasusa.out.fastq)
        
        ch_gfa      = ch_gfa.mix(raven.out.gfa)
        // no ch_summary
        // no ch_consensus
        ch_versions = ch_versions.mix(raven.out.versions.first())
      }
      
      if (params.assembler =~ /flye/ ) {
        flye(rasusa.out.fastq)

        flye.out.summary
          .collectFile(
            storeDir: "${params.outdir}/summary/",
            keepHeader: true,
            sort: { file -> file.text },
            name: "flye_summary.tsv")
          .set { flye_summary }

        ch_gfa      = ch_gfa.mix(flye.out.gfa)
        ch_summary  = ch_summary.mix(flye_summary)
        // no ch_consensus
        ch_versions = ch_versions.mix(flye.out.versions.first())
      }
    }

    bandage(ch_gfa)
    gfastats(ch_gfa)

    gfastats.out.summary
      .collectFile(
        storeDir: "${params.outdir}/summary/",
        keepHeader: true,
        sort: { file -> file.text },
        name: "gfastats_summary.csv")
      .set { gfastats_summary }

    ch_versions = ch_versions.mix(bandage.out.versions).mix(gfastats.out.versions)
    ch_summary  = ch_summary.mix(gfastats_summary).mix(bandage.out.png)

    if (params.assembler.replaceAll('dragonflye','dragon') =~ /flye/ || params.assembler =~ /raven/ ) {
      gfa_to_fasta(gfastats.out.stats.filter { it -> !(it[1] =~ /unicycler/ )} )

      dnaapler(gfa_to_fasta.out.fasta)

      ch_versions = ch_versions.mix(dnaapler.out.versions)

      dnaapler.out.fasta
        .branch {
          dragonflye: it =~ /dragonflye/
          raven: it =~ /raven/
          flye: it =~ /flye/
        }
        .set { ch_dnaapler_out }

      ch_dnaapler_out.flye
        .join(ch_nanopore_input, by:0, remainder: false)
        .mix(ch_dnaapler_out.raven.join(ch_nanopore_input, by:0, remainder: false))
        .set { ch_reoriented }

      medaka(ch_reoriented)

      ch_versions = ch_versions.mix(medaka.out.versions)

      medaka.out.fasta
        .branch {
          dragonflye: it =~ /dragonflye/
          raven: it =~ /raven/
          flye: it =~ /flye/
        }
        .set { ch_medaka_out }   

      ch_medaka_out.flye
        .join(ch_illumina_input, by:0, remainder: false)
        .mix(ch_medaka_out.raven.join(ch_illumina_input, by:0, remainder: false))
        .set { ch_medaka_polished }  

      bwa(ch_medaka_polished)
      polypolish(bwa.out.sam)

      ch_versions = ch_versions.mix(bwa.out.versions).mix(polypolish.out.versions)

      polypolish.out.fasta
        .branch {
          dragonflye: it =~ /dragonflye/
          raven: it =~ /raven/
          flye: it =~ /flye/
        }
        .set { ch_polypolish_out }   

      ch_polypolish_out.flye
        .join(ch_illumina_input, by:0, remainder: false)
        .mix(ch_polypolish_out.raven.join(ch_illumina_input, by:0, remainder: false))
        .set { ch_polypolish_polished }
           
      pypolca(ch_polypolish_polished)

      pypolca.out.summary
        .collectFile(name: "pypolca_summary.tsv",
          keepHeader: true,
          sort: { file -> file.text },
          storeDir: "${params.outdir}/summary")
        .set { pypolca_summary }

      ch_summary = ch_summary.mix(pypolca_summary)
      ch_versions = ch_versions.mix(pypolca.out.versions)
    
      ch_consensus = ch_consensus.mix(dnaapler.out.fasta).mix(medaka.out.fasta).mix(polypolish.out.fasta).mix(pypolca.out.fasta)
    }

    nanoplot(ch_nanopore_input)

    nanoplot.out.summary
      .collectFile(name: "nanoplot_summary.csv",
        keepHeader: true,
        storeDir: "${params.outdir}/summary")
      .set { nanostats_summary }

    ch_summary = ch_summary.mix(nanostats_summary).mix(nanoplot.out.stats)
    ch_versions = ch_versions.mix(nanoplot.out.versions)

    busco(ch_consensus)

    ch_summary = ch_summary.mix(busco.out.summary)
    ch_versions = ch_versions.mix(busco.out.versions)

    ch_consensus
      .filter{ it -> !(it[1] =~ /pypolca/ )}
      .filter{ it -> !(it[1] =~ /medaka/ )}
      .filter{ it -> !(it[1] =~ /poylpolish/ )}
      .branch {
        dragonflye: it =~ /dragonflye/
        raven: it =~ /raven/
        flye: it =~ /flye/
        unicycler: it =~ /unicycler/
      }
      .set { ch_assemblies }

    ch_assemblies.dragonflye
      .join(ch_nanopore_input, by: 0 , remainder: false).join(ch_illumina_input, by: 0, remainder: true)
      .mix(ch_assemblies.flye.join(ch_nanopore_input, by: 0 , remainder: false).join(ch_illumina_input, by: 0, remainder: true))
      .mix(ch_assemblies.unicycler.join(ch_nanopore_input, by: 0 , remainder: false).join(ch_illumina_input, by: 0, remainder: true))
      .mix(ch_assemblies.raven.join(ch_nanopore_input, by: 0 , remainder: false).join(ch_illumina_input, by: 0, remainder: true))
      .filter{ it -> if (it) {it[1]}}
      .set{ch_assembly_reads}

    circulocov(ch_assembly_reads)

    circulocov.out.summary
      .collectFile(name: "circulocov_summary.txt",
        keepHeader: true,
        storeDir: "${params.outdir}/summary")
      .set { circulocov_summary }

    ch_versions = ch_versions.mix(circulocov.out.versions)
    ch_summary = ch_summary.mix(circulocov.out.summary)

    ch_versions
      .collectFile(
        keepHeader: false,
        name: "versions.yml")
      .set { ch_collated_versions }

    versions(ch_collated_versions)
    ch_summary = ch_summary.mix(versions.out.versions)

    summary(ch_summary.unique().collect())

    multiqc(ch_summary.unique().collect())

    ch_consensus
      .combine(circulocov_summary)
      .combine(gfastats_summary)
      .set { ch_fasta_info }

    copy(ch_fasta_info)
  emit:
    fasta = ch_consensus
}


// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Workflow

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

workflow {

  if (params.test) {

    test_unicycler()

    test_unicycler.out.fastq
      .map { it ->
        meta = [id:it[0]] 
        tuple( meta,
          file("${it[1]}", checkIfExists: true),
          [file("${it[2][0]}", checkIfExists: true), file("${it[2][1]}", checkIfExists: true)])
      }
      .set{ ch_unicycler_out }

    test_donut_falls()

    test_donut_falls.out.fastq
      .map { it ->
        meta = [id:it[0]] 
        tuple( meta,
          file("${it[1]}", checkIfExists: true),
          [file("${it[2][0]}", checkIfExists: true), file("${it[2][1]}", checkIfExists: true)])
      }
      .set{ ch_test_df_out }

    test_donut_falls.out.lrfastq
      .map { it ->
        meta = [id:it[0]] 
        tuple( meta,
          file("${it[1]}", checkIfExists: true),
          null )
      }
      .set{ ch_test_dflr_out }

    ch_unicycler_out
      .mix(ch_test_df_out)
      .mix(ch_test_dflr_out)
      .set { ch_test }

    ch_test
      .map{it -> tuple(it[0], it[1])}
      .set { ch_test_nanopore }

    ch_test
      .filter{ it[2] }
      .map{it -> tuple(it[0], it[2])}
      .set { ch_test_illumina }

    ch_nanopore_input = ch_nanopore_input.mix(ch_test_nanopore)
    ch_illumina_input = ch_illumina_input.mix(ch_test_illumina)
  }

  if (params.sequencing_summary) {
    nanoplot_summary(ch_nanoplot_summary)
  }

  DONUT_FALLS(ch_nanopore_input, ch_illumina_input.ifEmpty([]))
}

workflow.onComplete {
  println("Pipeline completed at: $workflow.complete")
  println("The multiqc report can be found at ${params.outdir}/multiqc/multiqc_report.html")
  println("The consensus fasta files can be found in ${params.outdir}/consensus")
  println("The fasta files are from each phase of assembly. polca > polypolish > medaka > unpolished")
  println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}