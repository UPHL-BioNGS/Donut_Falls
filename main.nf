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

if (params.config_file) {
  def src = new File("${workflow.projectDir}/configs/donut_falls_config_template.config")
  def dst = new File("${workflow.launchDir}/edit_me.config")
  dst << src.text
  println("A config file can be found at ${workflow.launchDir}/edit_me.config")
  exit 0
}

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Checking params

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

def paramCheck(keys) {
  set_keys = [
    "input",
    "outdir",
    "sample_sheet",
    "sequencing_summary",
    "assembler",
    "test",
    "config_file"]

  for(key in keys){
    if (key !in set_keys){
      println("WARNING: ${key} isn't a supported param!")
      println("Supported params: ${set_keys}")
    }
  }
}

paramCheck(params.keySet())

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Input files

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

if (params.sequencing_summary){
  Channel
    .fromPath("${params.sequencing_summary}", type: 'file', checkIfExists: true)
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
        file("${it.fastq}", checkIfExists: true),
        "${it.fastq_1}",
        "${it.fastq_2}")
    }
    .set{ ch_input_files }
} else {
  ch_input_files = Channel.empty()
}

// channel for illumina files (paired-end only)
ch_input_files
  .filter { it[2] != it[3] }
  .map { it -> tuple(it[0], [file(it[2], checkIfExists: true), file(it[3], checkIfExists: true)])}
  .set { ch_illumina_input }

// channel for nanopore files
ch_input_files
  .map { it -> tuple (it[0], file(it[1], checkIfExists: true))}
  .set { ch_nanopore_input }

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Processes

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

process bandage {
  tag           "${meta.id}"
  label         'process_low'
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/bandage:0.8.1'
  time          '10m'

  input:
  tuple val(meta), file(gfa)

  output:
  path "bandage/*", emit: files
  tuple val(meta), file("bandage/*.png"), emit: png
  path "versions.yml", emit: versions
  
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
    bandage: \$(Bandage --version | awk '{print \$NF}')
  END_VERSIONS
  """
}

process busco {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/busco:5.8.2-prok-bacteria_odb12_2024-11-14'
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
  container     'staphb/bwa:0.7.18'
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
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/circulocov:0.1.20240104'
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
  label         'process_low'
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/multiqc:1.25'
  time          '10m'

  input:
  tuple val(meta), file(fasta), file(circulocov), file(gfastats)

  output:
  path "consensus/*", emit: fastas
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  """
  #!/usr/bin/env python3
  import glob
  import json
  import csv
  import os

  def gfastats_to_dict(header_dict):
    dict = {}
    with open('gfastats_summary.csv', mode='r') as file:
      reader = csv.DictReader(file)
      for row in reader:
        if row['sample'] == header_dict['name'] + '_' + header_dict['assembler']:
          key = row['Header']
          dict[key] = row
    return dict

  def circulocov_to_dict(header_dict):
    dict = {}
    with open('circulocov_summary.txt', mode='r', newline='') as file:
        reader = csv.DictReader(file, delimiter='\\t')
        for row in reader:
          if row['sample'].replace('_reoriented','') == header_dict['name'] + '_' + header_dict['assembler'] :
            key = row['contigs']

            dict[key] = row
    return dict

  def copy_fasta(fasta, header_dict, gfa_dict, circulocov_dict):
    with open(fasta, 'r') as file:
      fasta_dict = {}
      for line in file:
        line = line.strip()
        if line.startswith('>'):
          contig    = str(line.replace('>','').split()[0])
          circular  = gfa_dict[contig]['Is circular'].replace('Y','true').replace('N','false')
          length    = gfa_dict[contig]['Total segment length']
          gc_per    = gfa_dict[contig]['GC content %']
          meandepth = circulocov_dict[contig]['nanopore_meandepth']
          assembler = header_dict['assembler']
          step      = header_dict['step']
          
          # creating the header
          header = '>' + contig
          header = header + ' circ=' + circular
          header = header + ' len=' + length
          header = header + ' gc=' + gc_per
          header = header + ' cov=' + meandepth
          header = header + ' asmb=' + assembler 
          if assembler != 'unicycler':
            header = header + ' stp=' + step
          header = header + '\\n'

          # creating the dict
          fasta_dict[contig] = {}
          fasta_dict[contig]['seq']    = ''
          fasta_dict[contig]['header'] = header
          fasta_dict[contig]['length'] = int(length)

        else:
          fasta_dict[contig]['seq'] = fasta_dict[contig]['seq'] + line
    
    sorted_dict = dict(sorted(fasta_dict.items(), key=lambda item: item[1]['length'], reverse = True))

    with open(f"consensus/{header_dict['fasta']}", 'w') as outfile:    
      for contig in sorted_dict:
        seq = '\\n'.join([fasta_dict[contig]['seq'][i:i+70] for i in range(0, len(fasta_dict[contig]['seq']), 70)])
        outfile.write(fasta_dict[contig]['header'])
        outfile.write(seq + '\\n')

  def main():
    os.mkdir('consensus')
    header_dict = {}
    fasta = glob.glob('*.fasta')[0]
    header_dict['fasta'] = fasta

    name = fasta.replace('.fasta', '')

    assemblers = ['dragonflye', 'flye', 'hybracter', 'raven', 'unicycler']
    steps = ['reoriented', 'polypolish', 'pypolca', 'medaka']
    for step in steps:
      if step in name:
        header_dict['step'] = step
        name = name.replace(f"_{step}",'')
        break

    if 'step' not in header_dict.keys():
      header_dict['step'] = False

    for assembler in assemblers:
      if assembler in name:
        header_dict['assembler'] = assembler
        name = name.replace(f"_{assembler}",'')
        break

    header_dict['name'] = name

    gfa_dict        = gfastats_to_dict(header_dict)
    circulocov_dict = circulocov_to_dict(header_dict)

    copy_fasta(fasta, header_dict, gfa_dict, circulocov_dict)

  if __name__ == '__main__':
    main()
  """
}

process dnaapler {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/dnaapler:1.0.1'
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
    --ignore ignore.txt

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    dnaapler: \$(dnaapler --version | awk '{print \$NF}')
  END_VERSIONS
  """
}

process fastp {
  tag           "${meta.id}"
  label         "process_low"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/fastp:0.24.0'

  input:
  tuple val(meta), file(reads)

  output:
  tuple val(meta), file("fastp/*_fastp*.fastq.gz"), emit: fastq
  path "fastp/*", emit: everything
  path "fastp/*_fastp*.json", emit: summary
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
    -h fastp/${prefix}_fastp.html \
    -j fastp/${prefix}_fastp.json

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    fastp: \$(fastp --version 2>&1 | awk '{print \$NF}' )
  END_VERSIONS
  """
}

process fastplong {
  tag           "${meta.id}"
  label         "process_low"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/fastplong:0.2.2'

  input:
  tuple val(meta), file(reads)

  output:
  tuple val(meta), file("fastplong/*_fastp*.fastq.gz"), emit: fastq
  path "fastplong/*", emit: everything
  path "fastplong/*_fastp*.json", emit: summary
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p fastplong

  fastplong ${args} \
    --in ${reads} \
    --out fastplong/${prefix}_fastplong.fastq.gz \
    --html fastplong/${prefix}_fastplong.html \
    --json fastplong/${prefix}_fastplong.json

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    fastp: \$(fastp --version 2>&1 | awk '{print \$NF}')
  END_VERSIONS
  """
}

process flye {
  tag           "${meta.id}"
  label         "process_high"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/flye:2.9.5'
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
  def args      = task.ext.args   ?: ''
  def read_type = task.ext.read_type ?: '--nano-hq'
  def prefix    = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p flye

  flye ${args} \
    ${read_type} ${fastq} \
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
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/gfastats:1.3.7'
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

  head -n 1 gfastats/${prefix}_gfastats.txt | tr "\\t" "," | awk '{print "sample," \$0 }' > gfastats/${prefix}_gfastats_summary.csv
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
  container     'staphb/multiqc:1.26'
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
                      new_header = ">" + header + " length=" + summary_dict[header]['Total segment length'] + " circular=" + summary_dict[header]["Is circular"].replace("N","false").replace("Y","true") + " gc_per=" + summary_dict[header]["GC content %"] + "\\n"
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
                  if summary_dict[key]["Is circular"] == "N":
                      output_file.write(key + "\\n")
      return summary_dict

  gfastats_file = glob.glob("*_gfastats_summary.csv")
  gfa_file = glob.glob("*.gfa")

  summary_dict = read_summary_csv(gfastats_file[0])
  convert_to_fasta(summary_dict, gfa_file[0])
  """
}

// mash results : Reference-ID, Query-ID, Mash-distance, P-value, and Matching-hashes
process mash {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/mash:2.3'
  time          '30m'
  
  input:
  tuple val(meta), file(illumina), file(nanopore)

  output:
  tuple val(meta), env(dist), optional: true, emit: dist
  path "mash/*", optional: true, emit: txt
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args     = task.ext.args     ?: ''
  def ont_args = task.ext.ont_args ?: '-m 2'
  def ill_args = task.ext.ill_args ?: '-m 2'
  def short_re = "${illumina.join(' ')}"
  def prefix   = task.ext.prefix   ?: "${meta.id}"
  """
  mkdir mash

  cat ${short_re} | \
    mash sketch ${ill_args} \
    -o ${prefix}.illumina - 

  mash sketch ${ont_args} \
    -o ${prefix}.nanopore ${nanopore}

  mash dist ${args} \
      -p ${task.cpus} \
      ${prefix}.illumina.msh \
      ${prefix}.nanopore.msh | \
      awk -v prefix=${prefix} '{print prefix "\\t" \$0 }' \
      > mash/${prefix}.mashdist.txt

  dist=\$(head -n 1 mash/${prefix}.mashdist.txt | awk '{print \$4}')

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    mash: \$( mash --version )
  END_VERSIONS
  """
}

// From https://github.com/nanoporetech/medaka
// > It is not recommended to specify a value of --threads greater than 2 for medaka consensus since the compute scaling efficiency is poor beyond this.
// > Note also that medaka consensus may been seen to use resources equivalent to <threads> + 4 as an additional 4 threads are used for reading and preparing input data.
process medaka {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/medaka:2.0.1'
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
  publishDir    "${params.outdir}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/multiqc:1.25'
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
    echo "#     Is circular:" >> gfastats_mqc.csv
    echo "#         title: 'Is circular'" >> gfastats_mqc.csv
    echo "#         description: 'Is circular'" >> gfastats_mqc.csv
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

  multiqc ${args} \
    --outdir multiqc \
    .
  """
}

process nanoplot_summary {
  tag           "${summary}"
  label         "process_low"
  publishDir    "${params.outdir}/summary", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/nanoplot:1.42.0'
  time          '30m'

  input:
  file(summary)

  output:
  path "nanoplot/*", emit: final_directory
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
    nanoplot: \$(NanoPlot --version | awk '{print \$NF}')
  END_VERSIONS
  """
}

process nanoplot {
  tag           "${meta.id}"
  label         "process_low"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/nanoplot:1.42.0'
  time          '10m'

  input:
  tuple val(meta), file(fastq)

  output:
  path "nanoplot/*", emit: everything
  path "nanoplot/${meta.id}*_NanoStats.txt", emit: stats
  path "nanoplot/${meta.id}*_NanoStats.csv", emit: summary
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  def readtype = fastq.toList().size() > 1 ? '_illumina' : ''
  """
  mkdir -p nanoplot

  NanoPlot ${args} \
    --fastq ${fastq.join(' ')} \
    --threads ${task.cpus} \
    --prefix ${prefix}${readtype}_ \
    --tsv_stats \
    --outdir nanoplot

  echo "sample,\$(   cut -f 1 nanoplot/${prefix}${readtype}_NanoStats.txt | tr '\\n' ',' )" >  nanoplot/${prefix}${readtype}_NanoStats.csv
  echo "${prefix}${readtype},\$(cut -f 2 nanoplot/${prefix}${readtype}_NanoStats.txt | tr '\\n' ',' )" >> nanoplot/${prefix}${readtype}_NanoStats.csv

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    nanoplot: \$(NanoPlot --version | awk '{print \$NF}')
  END_VERSIONS
  """
}

process png {
  tag           "${meta.id}"
  label         "process_low"
  publishDir    "${params.outdir}/${meta.id}/bandage", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/multiqc:1.25'
  time          '10m'

  input:
  tuple val(meta), file(png)

  output:
  path("bandage_*_mqc.png"), emit: png

  when:
  task.ext.when == null || task.ext.when

  shell:
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  #!/usr/bin/env python3

  from PIL import Image, ImageDraw, ImageFont
  import glob
  import shutil
  import os

  def main():
    png_files = glob.glob("*.png")
    png_files.sort()

    if len(png_files) >= 2:

      images_with_titles = []
      for file in png_files:
        analysis = str(file).split('_')[-1].split('.')[0]
        image = Image.open(file)
        draw = ImageDraw.Draw(image)
        draw.text((10, 10), analysis, fill="black", font_size=100)
        images_with_titles.append(image)

      total_width = sum(image.width for image in images_with_titles)
      max_height  = max(image.height for image in images_with_titles)

      combined_image = Image.new("RGB", (total_width, max_height), color="white")

      offset = 0
      for image in images_with_titles:
        combined_image.paste(image, (offset, 0))
        offset += image.width

      combined_image.save("bandage_${prefix}_mqc.png")

      for image in images_with_titles:
          image.close()

    else:
      shutil.copy(png_files[0], "bandage_${prefix}_mqc.png")

  if __name__ == "__main__":
      main()

  """
}

process polypolish {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/polypolish:0.6.0'
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
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/pypolca:0.3.1'
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
  def args   = task.ext.args   ?: '--careful'
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
    pypolca: \$(pypolca --version | head -n 1 | awk '{print \$NF}')
  END_VERSIONS
  """
}

process rasusa {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/rasusa:2.1.0'
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

  rasusa \
    reads \
    ${args} \
    --output rasusa/${prefix}_rasusa.fastq.gz \
    ${fastq}

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    rasusa: \$(rasusa --version | awk '{print \$NF}' )
  END_VERSIONS
  """
}

process raven {
  tag           "${meta.id}"
  label         "process_high"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/raven:1.8.3'
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
  label         "process_low"
  publishDir    "${params.outdir}/summary", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/multiqc:1.25'
  time          '30m'
  
  input:
  file(input)

  output:
  path "donut_falls_summary.json", emit: summary
  path "donut_falls_summary.tsv", emit: csv
  
  when:
  task.ext.when == null || task.ext.when

  """
  #!/usr/bin/env python3

  import glob
  import json
  import csv
  from os.path import exists

  def busco_results():
      dict = {}
      files = glob.glob('short_summary*txt')
      for file in files:
        sample_analysis = file.split('.')[-2]
        with open(file, 'r') as f:
          for line in f:
            if 'C:' and 'S:' and 'D:' and 'F:' and 'M:' and 'n:' in line:
              dict[sample_analysis] = line.strip()
              break
      return dict

  def circulocov_results():
      dict = {}
      files = glob.glob('*overall_summary.txt')
      for file in files:
        sample_analysis = file.replace('_overall_summary.txt', '').replace('_reoriented', '')
        dict[sample_analysis] = {}
        with open(file, 'r') as f:
          for line in f:
            parts = line.split()
            if parts[2] == 'all':
              dict[sample_analysis]['coverage'] = parts[7]
              dict[sample_analysis]['nanopore_numreads'] = parts[5]
              if len(parts) >= 9:
                dict[sample_analysis]['illumina_numreads'] = parts[9]

            if parts[2] == 'missing':
              if len(parts) > 8:
                unmapped_illumina = parts[8]
              else:
                unmapped_illumina = 0

              dict[sample_analysis]['unmapped_nanopore'] = parts[4]
              dict[sample_analysis]['unmapped_illumina'] = unmapped_illumina
      return dict

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
          key = row[header] + '_' + row[header2]
          dict[row[header]][key] = row
      return dict

  def final_file(dict):
      with open('donut_falls_summary.json', 'w') as json_file:
        json.dump(dict, json_file, indent=4)

  def mash_file(file):
      dict = {}
      with open(file, mode = 'r') as file:
        reader = csv.DictReader(file, delimiter='\\t', fieldnames=['sample', 'illumina','nanopore', 'dist', 'pvalue', 'hash'])
        for row in reader:
          key = row['sample']
          dict[key] = row
      return dict

  def tsv_file(dict):
      final_dict = {}
      with open('donut_falls_summary.tsv', 'w') as tsv:
        i = 0
        sorted_keys = sorted(dict.keys())
        for key in sorted_keys:
          final_dict[key] = {}
          for result_key in ['name', 'number_of_reads', 'mean_read_length', 'mean_qual', 'total_illumina_reads', 'nanopore_illumina_mash_distance', 'assemblers']:
            final_dict[key][result_key] = dict[key][result_key]
          
          results = ['total_length', 'num_contigs', 'circ_contigs', 'coverage', 'unmapped_nanopore', 'unmapped_nanopore_pc', 'unmapped_illumina', 'unmapped_illumina_pc']

          for assembler in ['flye', 'raven']:
            if assembler in dict[key]['assemblers'].replace('dragonflye','dragon'):
              if assembler in dict[key].keys():
                for result in results + ['busco']:
                  final_dict[key][assembler + '_' + result] = dict[key][assembler][result]
                final_dict[key][assembler + '_busco_polished'] = dict[key][assembler]['busco_pypolca']
                final_dict[key][assembler + '_quality_before_polishing'] = dict[key][assembler]['Consensus_Quality_Before_Polishing']
                final_dict[key][assembler + '_QV_before_polishing'] = dict[key][assembler]['Consensus_QV_Before_Polishing']
              else:
                for result in results + ['quality_before_polishing', 'QV_before_polishing' ]:
                  final_dict[key][assembler + '_' + result] = 0

                for result in ['busco', 'busco_polished']:
                  final_dict[key][assembler + '_' + result] = 'NF'

          if 'unicycler' in dict[key]['assemblers']:
            if 'unicycler' in dict[key].keys():
              for result in results + [ 'busco']:
                final_dict[key]['unicycler_' + result] = dict[key]['unicycler'][result]
            else:
              for result in results:
                final_dict[key]['unicycler_' + result] = 0
              final_dict[key]['unicycler_busco'] = 'NF'

          w = csv.DictWriter(tsv, final_dict[key].keys(), delimiter='\\t')
          if i < 1 :
            w.writeheader()
            i = i+1
          w.writerow(final_dict[key])

  def main():
      if exists('nanoplot_summary.csv') :
        nanoplot_dict = file_to_dict('nanoplot_summary.csv', 'sample', ',')
      else:
        nanoplot_dict = {}

      if exists('mash_summary.tsv') :
        mash_dict = mash_file('mash_summary.tsv')
      else:
        mash_dict = {}

      if exists('pypolca_summary.tsv') :
        pypolca_dict  = file_to_dict('pypolca_summary.tsv', 'sample', '\t')
      else:
        pypolca_dict = {}

      if exists('gfastats_summary.csv') :
        gfastats_dict = file_to_dict_uniq('gfastats_summary.csv', 'sample', 'Header', ',')
      else:
        gfastats_dict = {}

      busco_dict = busco_results()

      circulocov_dict = circulocov_results()

      final_results = {}
      assemblers = ['dragonflye', 'flye', 'hybracter', 'raven', 'unicycler']
      for key in nanoplot_dict.keys():

        if key.endswith('_illumina'):
          actual_key = key.replace('_illumina', '')
          if actual_key not in final_results.keys():
            final_results[actual_key] = {}
          final_results[actual_key]['total_illumina_reads'] = int(nanoplot_dict[key]['number_of_reads'])
        else:
          if key not in final_results.keys():
            final_results[key] = {}

          if 'total_illumina_reads' not in final_results[key].keys():
            final_results[key]['total_illumina_reads'] = 0

          final_results[key]['name'] = key

          # from nanostas
          final_results[key]['number_of_reads']  = int(nanoplot_dict[key]['number_of_reads'])
          final_results[key]['mean_read_length'] = float(nanoplot_dict[key]['mean_read_length'])
          final_results[key]['mean_qual']        = float(nanoplot_dict[key]['mean_qual'])

          # from mash
          if key in mash_dict.keys():
            final_results[key]['nanopore_illumina_mash_distance'] = float(mash_dict[key]['dist'])
          else:
            final_results[key]['nanopore_illumina_mash_distance'] = 'NF'

          final_results[key]['assemblers'] = '${params.assembler}'

          # for each assembler
          for assembler in assemblers:
            if key + '_' + assembler in gfastats_dict.keys():
              final_results[key][assembler] = {}
              final_results[key][assembler]['assembler'] = assembler

              # gfastats results
              total_length = 0
              num_circular = 0
              for contig in gfastats_dict[key + '_' + assembler].keys():
                total_length = total_length + int(gfastats_dict[key + '_' + assembler][contig]['Total segment length'])
                if gfastats_dict[key + '_' + assembler][contig]['Is circular'] == 'Y':
                  num_circular = num_circular + 1

              final_results[key][assembler]['total_length'] = total_length
              final_results[key][assembler]['num_contigs']  = len(gfastats_dict[key + '_' + assembler].keys())
              final_results[key][assembler]['circ_contigs'] = num_circular

              # circulocov results
              if key + '_' + assembler in circulocov_dict.keys():
                if 'coverage' in circulocov_dict[key + '_' + assembler].keys():
                  final_results[key][assembler]['coverage']          = float(circulocov_dict[key + '_' + assembler]['coverage'])
                else:
                  final_results[key][assembler]['coverage']          = 'NF'

                if 'unmapped_nanopore' in circulocov_dict[key + '_' + assembler].keys():
                  final_results[key][assembler]['unmapped_nanopore']    = int(circulocov_dict[key + '_' + assembler]['unmapped_nanopore'])
                  final_results[key][assembler]['unmapped_nanopore_pc'] = float('{:.2f}'.format(int(final_results[key][assembler]['unmapped_nanopore']) / int(nanoplot_dict[key]['number_of_reads']) * 100))
                else:
                  final_results[key][assembler]['unmapped_nanopore']    = 'NF'
                  final_results[key][assembler]['unmapped_nanopore_pc'] = 'NF'

                if 'unmapped_illumina' in circulocov_dict[key + '_' + assembler].keys():
                  final_results[key][assembler]['unmapped_illumina'] = int(circulocov_dict[key + '_' + assembler]['unmapped_illumina'])
                  if 'total_illumina_reads' in final_results[key].keys() and final_results[key]['total_illumina_reads'] > 0:
                    final_results[key][assembler]['unmapped_illumina_pc'] = float('{:.2f}'.format(int(final_results[key][assembler]['unmapped_illumina']) / int(final_results[key]['total_illumina_reads']) * 100 ))
                  else:
                    final_results[key][assembler]['unmapped_illumina_pc'] = 0.0
                else:
                  final_results[key][assembler]['unmapped_illumina'] = 'NF'
                  final_results[key][assembler]['unmapped_illumina_pc'] = 'NF'

              # busco results
              if key + '_' + assembler in busco_dict.keys():
                final_results[key][assembler]['busco'] = busco_dict[key + '_' + assembler]
              elif key + '_' + assembler + '_reoriented' in busco_dict.keys():
                final_results[key][assembler]['busco'] = busco_dict[key + '_' + assembler + '_reoriented']
              else:
                final_results[key][assembler]['busco'] = 'NF'

              if assembler != 'unicycler':
                for step in ['polypolish', 'pypolca', 'medaka']:
                  if key + '_' + assembler + '_' + step in busco_dict.keys():
                    final_results[key][assembler]['busco_' + step ] = busco_dict[key + '_' + assembler + '_' + step]
                  else:
                    final_results[key][assembler]['busco_' + step ] = 'NF'

              # pypolca results
              if key + '_' + assembler in pypolca_dict.keys():
                if 'Consensus_Quality_Before_Polishing' in pypolca_dict[key + '_' + assembler].keys():
                  final_results[key][assembler]['Consensus_Quality_Before_Polishing'] = float(pypolca_dict[key + '_' + assembler]['Consensus_Quality_Before_Polishing'])
                else:
                  final_results[key][assembler]['Consensus_Quality_Before_Polishing'] = 'NF'
                if 'Consensus_QV_Before_Polishing' in pypolca_dict[key + '_' + assembler].keys():
                  final_results[key][assembler]['Consensus_QV_Before_Polishing']      = float(pypolca_dict[key + '_' + assembler]['Consensus_QV_Before_Polishing'])
                else:
                  final_results[key][assembler]['Consensus_QV_Before_Polishing']      = 'NF'

              elif assembler != 'unicycler':
                final_results[key][assembler]['Consensus_Quality_Before_Polishing'] = 0.0
                final_results[key][assembler]['Consensus_QV_Before_Polishing']      = 0.0

      final_file(final_results)
      tsv_file(final_results)

  if __name__ == '__main__':
    main()

  """
}

process unicycler {
  tag           "${meta.id}"
  label         'process_high'
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/unicycler:0.5.0'
  time          '10h'

  input:
  tuple val(meta), file(illumina), file(nanopore)

  output:
  tuple val(meta), file('unicycler/*_unicycler.fasta'), emit: fasta, optional: true
  tuple val(meta), file('unicycler/*_unicycler.gfa'), emit: gfa, optional: true
  path 'unicycler/*', emit: everything
  path 'versions.yml', emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  shell:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
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
    unicycler: \$(unicycler --version | awk '{print \$NF}' )
  END_VERSIONS
  """
}

process versions {
  tag           "extracting versions"
  label         "process_low"
  publishDir    "${params.outdir}/summary", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/multiqc:1.25'
  time          '30m'
  
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

process test {
  tag           "Downloading R10.4 reads"
  label         "process_low"
  publishDir    "${params.outdir}/test_files/", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/gfastats:1.3.7'
  time          '1h'

  output:
  tuple val("df"), file("test_files/test_nanopore.fastq.gz"), file("test_files/test_illumina_{1,2}.fastq.gz"), emit: fastq

  when:
  task.ext.when == null || task.ext.when

  shell:
  """
  wget --quiet https://zenodo.org/records/10779911/files/df_test_files.tar.gz?download=1 -O dataset.tar.gz
  tar -xvf dataset.tar.gz
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

    mash(ch_illumina_input.join(ch_nanopore_input, by: 0, remainder: false))

    ch_versions = ch_versions.mix(mash.out.versions)

    mash.out.txt
      .collectFile(
        storeDir: "${params.outdir}/summary/",
        keepHeader: false,
        sort: { file -> file.text },
        name: "mash_summary.tsv")
      .set { mash_summary }

    ch_summary  = ch_summary.mix(mash_summary)

    ch_illumina_input
      .join(mash.out.dist, by: 0)
      .filter{it[2] as float < 0.5}
      .map{it -> tuple(it[0], it[1])}
      .set {ch_dist_filter}

    if (params.assembler =~ /unicycler/ ) {
      unicycler(ch_dist_filter.join(ch_nanopore_input, by: 0, remainder: false))
      
      ch_gfa       = ch_gfa.mix(unicycler.out.gfa)
      // no ch_summary
      ch_consensus = ch_consensus.mix(unicycler.out.fasta)
      ch_versions  = ch_versions.mix(unicycler.out.versions.first())
    }

    if (params.assembler.replaceAll('dragonflye','dragon') =~ /flye/ || params.assembler =~ /raven/ ) {
      // quality filter

      fastp(ch_dist_filter.map { it -> [it[0], it[1]]}.filter{it[0]})

      ch_versions = ch_versions.mix(fastp.out.versions)
      ch_summary  = ch_summary.mix(fastp.out.summary)

      fastplong(ch_nanopore_input.map { it -> [it[0], it[1]]})

      rasusa(fastplong.out.fastq)

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
    ch_summary  = ch_summary.mix(gfastats_summary)

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
        .join(ch_dist_filter, by:0, remainder: false)
        .mix(ch_medaka_out.raven.join(ch_dist_filter, by:0, remainder: false))
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
        .join(ch_dist_filter, by:0, remainder: false)
        .mix(ch_polypolish_out.raven.join(ch_dist_filter, by:0, remainder: false))
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

    nanoplot(ch_nanopore_input.mix(ch_illumina_input.filter{it[1]}))

    nanoplot.out.summary
      .collectFile(name: "nanoplot_summary.csv",
        keepHeader: true,
        storeDir: "${params.outdir}/summary")
      .set { nanostats_summary }

    ch_summary = ch_summary.mix(nanostats_summary).mix(nanoplot.out.stats)
    ch_versions = ch_versions.mix(nanoplot.out.versions)

    busco(ch_consensus)

    ch_summary = ch_summary.mix(busco.out.summary)
    ch_versions = ch_versions.mix(busco.out.versions.first())

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
      .join(ch_nanopore_input, by: 0 , remainder: false).join(ch_dist_filter, by: 0, remainder: true)
      .mix(ch_assemblies.flye.join(ch_nanopore_input, by: 0 , remainder: false).join(ch_dist_filter, by: 0, remainder: true))
      .mix(ch_assemblies.unicycler.join(ch_nanopore_input, by: 0 , remainder: false).join(ch_dist_filter, by: 0, remainder: true))
      .mix(ch_assemblies.raven.join(ch_nanopore_input, by: 0 , remainder: false).join(ch_dist_filter, by: 0, remainder: true))
      .filter{ it -> if (it) {it[1]}}
      .set{ch_assembly_reads}

    png(bandage.out.png.groupTuple(by:0))
    ch_summary = ch_summary.mix(png.out.png)

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

    test()

    test.out.fastq
      .map { it ->
        meta = [id:it[0]] 
        tuple( meta,
          file("${it[1]}", checkIfExists: true),
          [file("${it[2][0]}", checkIfExists: true), file("${it[2][1]}", checkIfExists: true)])
      }
      .set{ ch_test_out }

    ch_test_out
      .map{it -> tuple(it[0], it[1])}
      .set { ch_test_nanopore }

    ch_test_out
      .filter{ it[2] }
      .map{it -> tuple(it[0], it[2])}
      .set { ch_test_illumina }

    ch_nanopore_input = ch_nanopore_input.mix(ch_test_nanopore)
    ch_illumina_input = ch_illumina_input.mix(ch_test_illumina)
  }

  if (params.sequencing_summary) {
    nanoplot_summary(ch_sequencing_summary)
  }

  DONUT_FALLS(ch_nanopore_input, ch_illumina_input.ifEmpty([]))
}

workflow.onComplete {
  println("Pipeline completed at: $workflow.complete")
  println("The multiqc report can be found at ${params.outdir}/multiqc/multiqc_report.html")
  println("The consensus fasta files can be found in ${params.outdir}/sample/consensus")
  println("The fasta files are from each phase of assembly. polca > polypolish > medaka > unpolished")
  println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}
