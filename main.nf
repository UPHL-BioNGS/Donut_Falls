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
  def src = new File("${workflow.projectDir}/configs/donut_falls_template.config")
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
    "config_file",
    "custom_config_version",
    "custom_config_base",
    "config_profile_name",
    "config_profile_description",
    "config_profile_contact",
    "config_profile_url"]

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

  script:
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

process bcftools {
  tag           "${meta.id}"
  label         'process_medium'
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/bcftools:1.22'
  time          '10m'

  input:
  tuple val(meta), file(fasta), file(vcf)

  output:
  tuple val(meta), file("clair3/*_clair3.fasta"), emit: fasta, optional: true
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ""
  def prefix = task.ext.prefix ?: "${fasta.baseName.replaceAll('_reoriented','')}"
  """
  mkdir -p clair3 logs/${task.process}

  bcftools index ${vcf}
  
  bcftools consensus \
    ${args} \
    -f ${fasta} \
    ${vcf} \
    -o clair3/${prefix}_clair3.fasta

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    bcftools: \$(bcftools --version | head -n 1 | awk '{print \$NF}')
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

  script:
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
  container     'staphb/bwa:0.7.19'
  time          '2h'

  input:
  tuple val(meta), file(fasta), file(fastq)

  output:
  tuple val(meta), file(fasta), file("bwa/*_{1,2}.sam"), emit: sam
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
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
  label         "process_high"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', pattern: "circulocov/*/*"
  container     'staphb/circulocov:0.1.20240104'
  time          '1h'

  input:
  tuple val(meta), file(fasta), file(nanopore), file(illumina)

  output:
  path "circulocov/*overall_summary.txt", emit: summary
  tuple val(meta), file("circulocov/*/overall_summary.txt"), emit: results
  tuple val(meta), file(fasta), file("circulocov/*/*map-ont.bam"), file("circulocov/*/*map-ont.bam.bai"), emit: bam
  path "circulocov/*/*", emit: everything
  path "circulocov/*/fastq/*", emit: fastq
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
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

process clair3 {
  tag           "${meta.id}"
  label         'process_medium'
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', pattern: "clair3/*"
  container     'staphb/clair3:1.1.0'
  time          '10m'

  input:
  tuple val(meta), file(fasta), file(bam), file(bai)

  output:
  tuple val(meta), file(fasta), path("clair3/merge_output.vcf.gz"), emit: vcf
  path("clair3/*"), emit: everything
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  script:
  def args = task.ext.args ?: '--include_all_ctgs'
  def prefix = task.ext.prefix ?: "${meta.id}"
  """
  mkdir -p clair3/${prefix}

  samtools faidx ${fasta}

  run_clair3.sh ${args} \
    --bam_fn=${bam[0]} \
    --ref_fn=${fasta} \
    --threads=${task.cpus} \
    --output=clair3 \
    --platform=ont \
    --model_path /clair3/models/ont

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    clair3: \$(run_clair3.sh  --version | awk '{print \$NF}')
  END_VERSIONS
  """
}

process copy {
  tag           "${meta.id}"
  label         'process_medium'
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/multiqc:1.28'
  time          '10m'

  input:
  tuple val(meta), file(fasta), file(gfastats)

  output:
  path "consensus/*", emit: fastas
  
  when:
  task.ext.when == null || task.ext.when

  script:
  """
  #!/usr/bin/env python3 
import glob
import os
import shutil

def check_circ(file, name):
      for x in [".fasta", "_reoriented", "_clair3", "_polypolish", "_pypolca"]:
          name = name.replace(x, "")

      with open(file, "r") as f:
          for line in f:
              line = line.strip()
              searching_line = line.split(",")
              if searching_line[0] == name and searching_line[-1] == "N":
                  return False
      return True


def sub_fasta(fasta):
      with open(fasta, 'r') as file:
          i = 0
          with open(f"consensus/sub_{fasta}", 'w') as outfile:
              for line in file:
                  line = line.strip()
                  if line.startswith('>') and i < 1:
                      outfile.write(f"{line.split()[0]} [location=chromosome][topology=circular][completeness=complete]\\n")
                      i += 1
                  elif line.startswith('>') and i >= 1:
                      outfile.write(f"{line.split()[0]} [plasmid-name=unnamed{i}][topology=circular][completeness=complete]\\n")
                  else:
                      outfile.write(f"{line}\\n")


def main():
    
      os.mkdir('consensus')
    
      fasta = glob.glob('*.fasta')[0]
    
      shutil.copy(fasta, f"consensus/{fasta}")
    
      if check_circ('gfastats_summary.csv', fasta):
    
          sub_fasta(fasta)

if __name__ == '__main__':
      main()
  """
}

process dnaapler {
  tag           "${meta.id}"
  label         "process_medium"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/dnaapler:1.2.0'
  time          '1h'

  input:
  tuple val(meta), file(fasta), file(ignore)

  output:
  tuple val(meta), file("dnaapler/*_reoriented.fasta"), emit: fasta
  path "dnaapler/*", emit: files
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
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
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/fastp:1.0.1'

  input:
  tuple val(meta), file(reads)

  output:
  tuple val(meta), file("fastp/*_fastp*.fastq.gz"), emit: fastq
  path "fastp/*", emit: everything
  path "fastp/*_fastp*.json", emit: summary
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  script:
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

  script:
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
    fastplong: \$(fastplong --version 2>&1 | awk '{print \$NF}')
  END_VERSIONS
  """
}

process flye {
  tag           "${meta.id}"
  label         "process_high"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/flye:2.9.6'
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

  script:
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
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', pattern: "gfastats/*"
  container     'staphb/gfastats:1.3.11'
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

  script:
  def args   = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${gfa.baseName}"
  """
  mkdir -p gfastats

  gfastats \
    ${gfa} \
    ${args} \
    --threads ${task.cpus} \
    --tabular \
    --segment-report \
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
  container     'staphb/multiqc:1.28'
  time          '10m'

  input:
  tuple val(meta), file(gfa), file(stats)

  output:
  tuple val(meta), file("*fasta"), file("noncircular.txt"), emit: fasta
  
  when:
  task.ext.when == null || task.ext.when

  script:
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
                      new_header = f">{header} length={summary_dict[header]['Total segment length']} circular={summary_dict[header]['Is circular'].replace('N','false').replace('Y','true')} gc_per={summary_dict[header]['GC content %']}\\n"
                      with open(outfile, mode='a') as output_file:
                          output_file.write(new_header)
                          output_file.write(f"{seq}\\n")

  def read_summary_csv(gfastats_file):
      summary_dict = {}
      with open(gfastats_file, mode='r', newline='') as file:
          reader = csv.DictReader(file)
          for row in reader:
              key = row['Header']
              summary_dict[key] = row
              with open("noncircular.txt", mode='a') as output_file:
                  if summary_dict[key]["Is circular"] == "N":
                      output_file.write(f"{key}\\n")
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
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', pattern: "mash/*.mashdist.txt"
  container     'staphb/mash:2.3'
  time          '30m'
  
  input:
  tuple val(meta), file(illumina), file(nanopore)

  output:
  tuple val(meta), file("*.head.mashdist.txt"), optional: true, emit: dist
  path "mash/*", optional: true, emit: txt
  path "versions.yml", emit: versions
  
  when:
  task.ext.when == null || task.ext.when

  script:
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

  head -n 1 mash/${prefix}.mashdist.txt > ${prefix}.head.mashdist.txt

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    mash: \$( mash --version )
  END_VERSIONS
  """
}

process multiqc {
  tag           "combining reports"
  label         "process_low"
  publishDir    "${params.outdir}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/multiqc:1.28'
  time          '10m'

  input:
  file(input)

  output:
  path "multiqc/multiqc_report.html", emit: report
  path "multiqc/multiqc_data/*", emit: everything
  
  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: ''
  """
  multiqc ${args} \
    --outdir multiqc \
    .
  """
}

process png {
  tag           "${meta.id}"
  label         "process_low"
  publishDir    "${params.outdir}/${meta.id}/bandage", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/multiqc:1.28'
  time          '10m'

  input:
  tuple val(meta), file(png)

  output:
  path("bandage_*_mqc.png"), emit: png

  when:
  task.ext.when == null || task.ext.when

  script:
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

  script:
  def args   = task.ext.args   ?: ''
  def filarg = task.ext.args   ?: ''
  def prefix = task.ext.prefix ?: "${fasta.baseName.replaceAll('_clair3','')}"
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

  script:
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

  script:
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

  script:
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

process seqkit {
  tag           "${meta.id}"
  label         "process_low"
  publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/seqkit:2.10.0'
  time          '10m'

  input:
  tuple val(meta), file(fastq)

  output:
  path("seqkit/*seqkit_stats.tsv"), emit: stats
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: '--all'
  def prefix = task.ext.prefix ?: "${meta.id}"
  def ill_rds = fastq[1] ? "${fastq[1]} ${fastq[2]}" : ""
  """
  mkdir -p seqkit

  seqkit stats \
    ${args} \
    --tabular \
    --threads ${task.cpus} \
    ${fastq[0]} \
    ${ill_rds} | \
    awk '{print "${prefix}\\t" \$0}' | \
    sed 's/${prefix}\\tfile/sample\\tfile/g' > \
    seqkit/${prefix}_seqkit_stats.tsv

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    seqkit: \$(seqkit version | sed 's/v//g' | awk '{print \$NF}')
  END_VERSIONS
  """  
}

process sort {
  tag           "${meta.id}"
  label         "process_low"
  // no publishdir by default
  container     'staphb/seqkit:2.10.0'
  time          '10m'

  input:
  tuple val(meta), file(fasta), file(stats)

  output:
  tuple val(meta), file("*/*.fasta"), file(stats), emit: fasta
  path "versions.yml", emit: versions

  when:
  task.ext.when == null || task.ext.when

  script:
  def args   = task.ext.args   ?: '--two-pass'
  def prefix = task.ext.prefix ?: "${fasta.baseName}"
  """
  mkdir -p seqkit

  seqkit sort \
    ${args} \
    --by-length \
    --reverse \
    ${fasta} > \
    seqkit/${prefix}.fasta

  cat <<-END_VERSIONS > versions.yml
  "${task.process}":
    seqkit: \$(seqkit version | sed 's/v//g' | awk '{print \$NF}')
  END_VERSIONS
  """  
}


process summary {
  tag           "Creating summary"
  label         "process_low"
  publishDir    "${params.outdir}/summary", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/multiqc:1.28'
  time          '30m'
  
  input:
  file(input)

  output:
  path "donut_falls_summary.json", emit: summary
  path "donut_falls_summary.tsv", emit: csv
  
  when:
  task.ext.when == null || task.ext.when

  script:
  """
  #!/usr/bin/env python3
  import glob
  import json
  import csv
  from os.path import exists

  def busco_results():
      results = {}
      files = glob.glob('short_summary*txt')
      for file in files:
        if len(file.split('.')[-2].split('_')) > 2:
          sample, assembler, step = file.split('.')[-2].rsplit("_", 2)
        else:
          sample, assembler = file.split('.')[-2].rsplit("_", 1)
          step = 'unicycler'
        if sample not in results.keys():
          results[sample] = {}
        if assembler not in results[sample].keys():
          results[sample][assembler] = {}      
        with open(file, 'r') as f:
          for line in f:
            if 'C:' and 'S:' and 'D:' and 'F:' and 'M:' and 'n:' in line:
              results[sample][assembler][step] = line.strip()
              break
      return results

  def circulocov_results():
      results = {}
      files = glob.glob('*_reoriented_overall_summary.txt') + glob.glob('*_unicycler_overall_summary.txt')
      for file in files:
        sample, assembler = file.replace('_reoriented', '').replace('_overall_summary.txt','').rsplit("_", 1)
        if sample not in results.keys():
          results[sample] = {}
        if assembler not in results[sample].keys():
          results[sample][assembler] = {}

        with open(file, 'r') as f:
          reader = csv.DictReader(f, delimiter="\\t")
          for row in reader:
            results[sample][assembler][row['contigs']] = row

        results[sample][assembler]['all']['unmapped_nanopore'] = results[sample][assembler]['missing']['nanopore_numreads'] if 'nanopore_numreads' in results[sample][assembler]['missing'].keys() else 0
        results[sample][assembler]['all']['unmapped_nanopore_pc'] = round(float(results[sample][assembler]['all']['unmapped_nanopore']) / float(results[sample][assembler]['all']['nanopore_numreads']), 2)

        if 'illumina_numreads' in results[sample][assembler]['missing'].keys():
          results[sample][assembler]['all']['unmapped_illumina'] = results[sample][assembler]['missing']['illumina_numreads']
          results[sample][assembler]['all']['unmapped_illumina_pc'] = round(float(results[sample][assembler]['all']['unmapped_illumina']) / float(results[sample][assembler]['all']['illumina_numreads']), 2)

      return results

  def combine_results(seqkit_dict, mash_dict, pypolca_dict, gfastats_dict, busco_dict, circulocov_dict):
        final_results = seqkit_dict

        for key in final_results.keys():

          if key in mash_dict.keys():
            final_results[key]['mash'] = mash_dict[key]

          if key in pypolca_dict.keys():
            final_results[key]['pypolca'] = pypolca_dict[key]

          if key in gfastats_dict.keys():
            final_results[key]['gfastats'] = gfastats_dict[key]

          if key in busco_dict.keys():
            final_results[key]['busco'] = busco_dict[key]

          if key in circulocov_dict.keys():
            final_results[key]['circulocov'] = circulocov_dict[key]

          final_results[key]['assemblers'] = "${params.assembler}"

        return final_results
          
  def final_file(results):
      with open('donut_falls_summary.json', 'w') as json_file:
        json.dump(results, json_file, indent=4)

  def gfastats_file(file):
      results = {}
      with open(file, mode='r', newline='') as file:
        reader = csv.DictReader(file, delimiter=",")
        for row in reader:
          sample, assembler = row['sample'].rsplit("_", 1)
          if sample not in results.keys():
            results[sample] = {}
          if assembler not in results[sample].keys():
            results[sample][assembler] = {}
          results[sample][assembler][row['Header']] = row

      for sample in results.keys():
        for assembler in results[sample].keys():
          total_length = 0
          num_circular = 0
          contigs = list(results[sample][assembler].keys())
          results[sample][assembler]['num_contigs'] = len(contigs)
          for contig in contigs:
            total_length = total_length + int(results[sample][assembler][contig]['Total segment length'])
            if results[sample][assembler][contig]['Is circular'] == 'Y':
              num_circular = num_circular + 1
          results[sample][assembler]['total_length'] = total_length
          results[sample][assembler]['circ_contigs'] = num_circular
      return results

  def mash_file(file):
      results = {}
      with open(file, mode = 'r') as file:
        reader = csv.DictReader(file, delimiter='\\t', fieldnames=['sample', 'illumina','nanopore', 'dist', 'pvalue', 'hash'])
        for row in reader:
          key = row['sample']
          results[key] = row
      return results

  def pypolca_file(file):
      results = {}
      with open(file, mode='r', newline='') as file:
        reader = csv.DictReader(file, delimiter='\\t')
        for row in reader:
          sample, assembler = row['sample'].rsplit("_", 1)
          results[sample] = { assembler : row }
      return results

  def seqkit_file(file):
      results = {}
      with open(file, mode='r', newline='') as file:
        reader = csv.DictReader(file, delimiter='\\t')
        for row in reader:
          key = row['sample']
          fastq = row['file']
          if key not in results.keys():
            results[key] = {}
          results[key][fastq] = row

      for key in results.keys():
        sorted_results = dict(sorted(results[key].items(), key=lambda fastq: float(fastq[1]['avg_len']), reverse=True))
        results[key] = sorted_results
        nanopore_fastq = list(results[key].keys())[0]
        results[key][nanopore_fastq]['platform'] = 'nanopore'
        results[key]['nanopore'] = results[key][nanopore_fastq]
        del results[key][nanopore_fastq]

        ill_fastq = []
        results_keys = list(results[key].keys())      
        for fastq in results_keys:
          if 'platform' not in results[key][fastq].keys():
            results[key][fastq]['platform'] = 'illumina'
            ill_fastq.append(fastq)
            results[key]['illumina'] = { 'platform' : 'illumina' }

        if len(ill_fastq) > 1 :
          for seqkit_val in ['num_seqs', 'sum_len', 'sum_gap', 'N50_num', 'sum_n']:
            results[key]['illumina'][seqkit_val] = round(float(results[key][ill_fastq[0]][seqkit_val]) + float(results[key][ill_fastq[1]][seqkit_val]), 0)

          for seqkit_val in ['avg_len', 'Q1', 'Q2', 'Q3', 'N50', 'Q20(%)', 'Q30(%)', 'AvgQual', 'GC(%)']:
            results[key]['illumina'][seqkit_val] = round((float(results[key][ill_fastq[0]][seqkit_val]) + float(results[key][ill_fastq[1]][seqkit_val]))/2, 2)

          results[key]['illumina']['file'] = f"{results[key][ill_fastq[0]]['file']},{results[key][ill_fastq[1]]['file']}"

          results[key]['illumina']['min_len'] = min(
            float(results[key][ill_fastq[0]]['min_len']), 
            float(results[key][ill_fastq[1]]['min_len']))
          
          results[key]['illumina']['max_len'] = max(
            float(results[key][ill_fastq[0]]['max_len']), 
            float(results[key][ill_fastq[1]]['max_len']))
          
          del results[key][ill_fastq[0]]
          del results[key][ill_fastq[1]]
      return results

  def tsv_file(results_dict):
      final_results_dict = {}
      
      # converting the final dict to something tsv friendly
      sorted_keys = list(sorted(results_dict.keys()))
      all_keys = []
      for sample in sorted_keys:
        final_results_dict[sample] = {}
        if 'nanopore' in results_dict[sample].keys():
          for result in results_dict[sample]['nanopore'].keys():
            final_results_dict[sample][f"seqkit_{result}"] = results_dict[sample]['nanopore'][result]

        if 'illumina' in results_dict[sample].keys():
          for result in results_dict[sample]['illumina'].keys():
            final_results_dict[sample][f"seqkit_illumina_{result}"] = results_dict[sample]['illumina'][result]

        if 'mash' in results_dict[sample].keys():
          for result in results_dict[sample]['mash'].keys():
            final_results_dict[sample][f"mash_{result}"] = results_dict[sample]['mash'][result]

        for analysis in ['pypolca', 'busco']:
          if analysis in results_dict[sample].keys():
            for assembler in results_dict[sample][analysis].keys():
              for result in results_dict[sample][analysis][assembler].keys():
                final_results_dict[sample][f"{assembler}_{analysis}_{result}"] = results_dict[sample][analysis][assembler][result]

        if 'gfastats' in results_dict[sample].keys():
          for assembler in results_dict[sample]['gfastats'].keys():
            for result in ['num_contigs', 'total_length', 'circ_contigs']:
              final_results_dict[sample][f"{assembler}_gfastats_{result}"] = results_dict[sample]['gfastats'][assembler][result]

        if 'circulocov' in results_dict[sample].keys():
          for assembler in results_dict[sample]['circulocov'].keys():
            for result in results_dict[sample]['circulocov'][assembler]['all'].keys():
              final_results_dict[sample][f"{assembler}_circulocov_{result}"] = results_dict[sample]['circulocov'][assembler]['all'][result]
        
        all_keys += list(final_results_dict[sample].keys())

      unique_key = list(set(all_keys))

      with open('donut_falls_summary.tsv', 'w') as tsv:
        i = 0
        for sample in sorted_keys:
          for key in unique_key:
            if key not in final_results_dict[sample].keys():
              final_results_dict[sample][key] = ""

          final_results_dict[sample] = { "sample": sample, **dict(sorted(final_results_dict[sample].items()))}

          possible_fieldnames = [
            "sample",
            "seqkit_num_seqs",
            "seqkit_avg_len",
            "seqkit_AvgQual",
            "seqkit_GC(%)",
            "mash_dist",
            "flye_gfastats_total_length",
            "flye_gfastats_num_contigs",
            "flye_gfastats_circ_contigs",
            "flye_circulocov_nanopore_meandepth",
            "flye_circulocov_unmapped_nanopore_pc",
            "flye_circulocov_illumina_meandepth",
            "flye_circulocov_unmapped_illumina_pc",
            "flye_busco_reoriented",
            "flye_busco_clair3",
            "flye_busco_polypolish",
            "flye_busco_pypolca",
            "flye_pypolca_Insertion/Deletion_Errors_Found",
            "flye_pypolca_Substitution_Errors_Found",
            "raven_gfastats_total_length",
            "raven_gfastats_num_contigs",
            "raven_gfastats_circ_contigs",
            "raven_circulocov_nanopore_meandepth",
            "raven_circulocov_unmapped_nanopore_pc",
            "raven_circulocov_illumina_meandepth",
            "raven_circulocov_unmapped_illumina_pc",
            "raven_busco_reoriented",
            "raven_busco_clair3",
            "raven_busco_polypolish",
            "raven_busco_pypolca",
            "raven_pypolca_Insertion/Deletion_Errors_Found",
            "raven_pypolca_Substitution_Errors_Found",
            "unicycler_gfastats_total_length",
            "unicycler_gfastats_num_contigs",
            "unicycler_gfastats_circ_contigs",
            "unicycler_circulocov_nanopore_meandepth",
            "unicycler_circulocov_unmapped_nanopore_pc",
            "unicycler_circulocov_illumina_meandepth",
            "unicycler_circulocov_unmapped_illumina_pc",
            "unicycler_busco_unicycler"
          ]

          fieldnames = []
          for name in possible_fieldnames:
            if name in final_results_dict[sample].keys():
              fieldnames.append(name)

          w = csv.DictWriter(tsv, fieldnames=fieldnames, delimiter='\\t')
          if i < 1 :
            w.writeheader()
            i += 1
          filtered_row = {k: final_results_dict[sample][k] for k in fieldnames}
          w.writerow(filtered_row)

  def main():

      seqkit_dict     = seqkit_file('seqkit_summary.tsv') if exists('seqkit_summary.tsv') else {}
      
      if not seqkit_dict:
        print('FATAL : Something is wrong and seqkit results were not located.')
        exit(1)

      mash_dict       = mash_file('mash_summary.tsv') if exists('mash_summary.tsv') else {}
      pypolca_dict    = pypolca_file('pypolca_summary.tsv') if exists('pypolca_summary.tsv') else {}
      gfastats_dict   = gfastats_file('gfastats_summary.csv') if exists('gfastats_summary.csv') else {}
      busco_dict      = busco_results()
      circulocov_dict = circulocov_results()

      final_results   = combine_results(seqkit_dict, mash_dict, pypolca_dict, gfastats_dict, busco_dict, circulocov_dict)

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
  container     'staphb/unicycler:0.5.1'
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

  script:
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
  container     'staphb/multiqc:1.28'
  time          '30m'
  
  input:
  file(input)

  output:
  path "software_versions_mqc.yml", emit: versions
  path "software_versions.yml", emit: yml

  when:
  task.ext.when == null || task.ext.when

  script:
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

// TODO : Add seqkit watch for visualizations
// process watch {
//   tag           "${meta.id}"
//   label         "process_low"
//   publishDir    "${params.outdir}/${meta.id}", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
//   container     'staphb/seqkit:2.10.0'
//   time          '10m'

//   input:
//   file(fastq)

//   output:
//   path("seqkit/*.tsv"), emit: stats
//   path("seqkit/*hist"), emit: histogram
//   path "versions.yml", emit: versions

//   when:
//   task.ext.when == null || task.ext.when

//   script:
//   def args   = task.ext.args   ?: '--all'
//   def prefix = task.ext.prefix ?: ""
//   """
//   mkdir -p seqkit

//   seqkit watch -p 500 --fields ReadLen *.fastq.gz -y -B 50

//   cat <<-END_VERSIONS > versions.yml
//   "${task.process}":
//     seqkit: \$(seqkit version | sed 's/v//g' | awk '{print \$NF}')
//   END_VERSIONS

//   exit 1
//   """  
// }

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Downloading files for testing

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

process test {
  tag           "Downloading R10.4 reads"
  label         "process_low"
  publishDir    "${params.outdir}/test_files/", mode: 'copy', saveAs: { filename -> filename.equals('versions.yml') ? null : filename }
  container     'staphb/gfastats:1.3.10'
  time          '1h'

  output:
  tuple val("df"), file("test_files/test_nanopore.fastq.gz"), file("test_files/test_illumina_{1,2}.fastq.gz"), emit: fastq

  when:
  task.ext.when == null || task.ext.when

  script:
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
    ch_gfa           = Channel.empty()
    // channel for files for multiqc or workflow summary
    ch_summary       = Channel.empty()
    // channel for all assembled genomes at all phases
    ch_consensus     = Channel.empty()
    // channel for de novo assembled genomes after reorienting
    ch_reoriented    = Channel.empty()
    // channel for de novo assembled genomes after clair3 polishing
    ch_clair3_fa     = Channel.empty()
    // channel for de novo assembled genomes after polypolish polishing
    ch_polypolish_fa = Channel.empty()
    // channel for de novo assembled genomes after pypolca polishing
    ch_pypolca_fa    = Channel.empty()
    // channel for hybrid assembled genomes with unicycler
    ch_unicycler_fa  = Channel.empty()
    // versions channel
    ch_versions      = Channel.empty()

    // get the mash distance between illumina and nanopore reads
    // this is helpful because these should be very close
    mash(ch_illumina_input.join(ch_nanopore_input, by: 0, remainder: false))

    mash.out.txt
      .collectFile(
        storeDir: "${params.outdir}/summary/",
        keepHeader: false,
        sort: { file -> file.text },
        name: "mash_summary.tsv")
      .set { mash_summary }

    mash.out.dist
      .splitCsv(sep: '\t')
      .map { it -> tuple(it[0], it[1][3]) }
      .set {ch_mash_dist}

    ch_summary  = ch_summary.mix(mash_summary)
    ch_versions = ch_versions.mix(mash.out.versions.first())

    // general qc information
    seqkit(ch_nanopore_input.mix(ch_illumina_input.transpose()).groupTuple().filter{it})

    seqkit.out.stats
      .collectFile(
        storeDir: "${params.outdir}/summary/",
        keepHeader: true,
        sort: { file -> file.text },
        name: "seqkit_summary.tsv")
      .set { seqkit_summary }

    ch_summary = ch_summary.mix(seqkit_summary)
    ch_versions = ch_versions.mix(seqkit.out.versions)

    ch_illumina_input
      .join(ch_mash_dist, by: 0)
      .filter{it[2] as float < 0.5}
      .map{it -> tuple(it[0], it[1])}
      .set {ch_dist_filter}

    if (params.assembler =~ /unicycler/ ) {
      unicycler(ch_dist_filter.join(ch_nanopore_input, by: 0, remainder: false))
      
      ch_gfa          = ch_gfa.mix(unicycler.out.gfa)
      // no ch_summary
      ch_consensus    = ch_consensus.mix(unicycler.out.fasta)
      ch_unicycler_fa = ch_unicycler_fa.mix(unicycler.out.fasta)
      ch_versions     = ch_versions.mix(unicycler.out.versions.first())
    }

    if (params.assembler =~ /flye/ || params.assembler =~ /myloasm/ || params.assembler =~ /raven/ ) {
      // quality filter

      fastp(ch_dist_filter.map { it -> [it[0], it[1]]}.filter{it[0]})
      ch_versions = ch_versions.mix(fastp.out.versions.first())
      ch_summary  = ch_summary.mix(fastp.out.summary)

      fastplong(ch_nanopore_input.map { it -> [it[0], it[1]]})
      ch_versions = ch_versions.mix(fastplong.out.versions.first())
      ch_summary  = ch_summary.mix(fastplong.out.summary)

      rasusa(fastplong.out.fastq)

      ch_versions = ch_versions.mix(rasusa.out.versions)

      if (params.assembler =~ /raven/ ) {
        raven(rasusa.out.fastq)
        
        ch_gfa      = ch_gfa.mix(raven.out.gfa)
        // no ch_summary
        // no ch_consensus
        ch_versions = ch_versions.mix(raven.out.versions.first())
      }

      // if (params.assembler =~ /myloasm/ ) {
      //   myloasm(rasusa.out.fastq)
        
      //   ch_gfa      = ch_gfa.mix(myloasm.out.gfa)
      //   ch_versions = ch_versions.mix(myloasm.out.versions.first())
      // }

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
    ch_versions = ch_versions.mix(bandage.out.versions.first())

    gfastats(ch_gfa)

    gfastats.out.summary
      .collectFile(
        storeDir: "${params.outdir}/summary/",
        keepHeader: true,
        sort: { file -> file.text },
        name: "gfastats_summary.csv")
      .set { gfastats_summary }

    ch_versions = ch_versions.mix(gfastats.out.versions.first())
    ch_summary  = ch_summary.mix(gfastats_summary)

    if (params.assembler =~ /flye/ || params.assembler =~ /myloasm/ || params.assembler =~ /raven/ ) {
      gfa_to_fasta(gfastats.out.stats.filter { it -> !(it[1] =~ /unicycler/ )} )

      sort(gfa_to_fasta.out.fasta)
      ch_versions = ch_versions.mix(sort.out.versions)

      dnaapler(sort.out.fasta)
      ch_versions = ch_versions.mix(dnaapler.out.versions)

      ch_consensus  = ch_consensus.mix(dnaapler.out.fasta)
      ch_reoriented = ch_reoriented.mix(dnaapler.out.fasta)
    }

    ch_reoriented
      .branch {
        myloasm: it =~ /myloasm/
        raven: it =~ /raven/
        flye: it =~ /flye/
      }
      .set { ch_reoriented_out }

    ch_reoriented_out.flye
      .join(ch_nanopore_input, by: 0 , remainder: false).join(ch_dist_filter, by: 0, remainder: true)
      .mix(ch_reoriented_out.raven.join(ch_nanopore_input, by: 0 , remainder: false).join(ch_dist_filter, by: 0, remainder: true))
      .mix(ch_unicycler_fa.join(ch_nanopore_input, by: 0 , remainder: false).join(ch_dist_filter, by: 0, remainder: true))
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

    if (params.assembler =~ /flye/ || params.assembler =~ /myloasm/ || params.assembler =~ /raven/ ) {
      clair3(circulocov.out.bam.filter { it -> !(it[1] =~ /unicycler/ )} )
      ch_versions = ch_versions.mix(clair3.out.versions.first())

      bcftools(clair3.out.vcf)
      ch_versions = ch_versions.mix(bcftools.out.versions)
      ch_consensus = ch_consensus.mix(bcftools.out.fasta)

      ch_clair3_fa = ch_clair3_fa.mix(bcftools.out.fasta)

      ch_clair3_fa
        .branch {
          myloasm: it =~ /myloasm/
          raven: it =~ /raven/
          flye: it =~ /flye/
        }
        .set { ch_clair3_out }   

      ch_clair3_out.flye
        .join(ch_dist_filter, by:0, remainder: false)
        .mix(ch_clair3_out.raven.join(ch_dist_filter, by:0, remainder: false))
        .set { ch_clair3_polished }

      bwa(ch_clair3_polished)
      ch_versions = ch_versions.mix(bwa.out.versions.first())

      polypolish(bwa.out.sam)
      ch_versions = ch_versions.mix(polypolish.out.versions.first())
      ch_polypolish_fa = ch_polypolish_fa.mix(polypolish.out.fasta)
      ch_consensus = ch_consensus.mix(ch_polypolish_fa)

      polypolish.out.fasta
          .branch {
            myloasm: it =~ /myloasm/
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
    
      ch_consensus = ch_consensus.mix(pypolca.out.fasta)
      ch_pypolca_fa = ch_pypolca_fa.mix(pypolca.out.fasta)
    }

    busco(ch_consensus)

    ch_summary = ch_summary.mix(busco.out.summary)
    ch_versions = ch_versions.mix(busco.out.versions.first())

    png(bandage.out.png.groupTuple(by:0))
    ch_summary = ch_summary.mix(png.out.png)

    ch_versions
      .collectFile(
        keepHeader: false,
        name: "versions.yml")
      .set { ch_collated_versions }

    versions(ch_collated_versions)
    ch_summary = ch_summary.mix(versions.out.versions)

    summary(ch_summary.unique().collect())

    multiqc(ch_summary.unique().collect())

    copy(ch_consensus.combine(gfastats_summary))

  emit:
    gfa              = ch_gfa
    consensus        = ch_consensus
    reoriented       = ch_reoriented
    clair3_polished  = ch_clair3_fa
    polypolished     = ch_polypolish_fa
    pypolca_polished = ch_pypolca_fa
    unicycler        = ch_unicycler_fa
    versions         = ch_versions
}

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

// Workflow

// ##### ##### ##### ##### ##### ##### ##### ##### ##### #####

workflow {

  if (params.test) {

    test()

    test.out.fastq
      .map { it ->
        def meta = [id:it[0]] 
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

  DONUT_FALLS(ch_nanopore_input, ch_illumina_input.ifEmpty([]))
}

workflow.onComplete {
  println("Pipeline completed at: $workflow.complete")
  println("The multiqc report can be found at ${params.outdir}/multiqc/multiqc_report.html")
  println("The consensus fasta files can be found in ${params.outdir}/sample/consensus")
  println("The fasta files are from each phase of assembly: unpolished/reoriented -> clair3 -> polypolish (if illumina reads are supplied) -> pypolca")
  println("Execution status: ${ workflow.success ? 'OK' : 'failed' }")
}
