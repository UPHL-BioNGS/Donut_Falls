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

// process test_unicycler {
//   tag           "Downloading Unicycler test files"
//   label         "process_single"
//   publishDir    "${params.outdir}/test_files/unicycler", mode: 'copy'
//   container     'staphb/multiqc:1.19'
//   errorStrategy { task.attempt < 2 ? 'retry' : 'ignore'}
//   time          '1h'

//   output:
//   tuple val("unicycler"), file("long_reads_low_depth.fastq.gz"), file("short_reads*.fastq.gz"), emit: fastq

//   when:
//   task.ext.when == null || task.ext.when

//   shell:
//   """
//   wget --quiet https://github.com/rrwick/Unicycler/raw/69e712eb95c4b9f8a46aade467260260a9ce7a91/sample_data/short_reads_1.fastq.gz
//   wget --quiet https://github.com/rrwick/Unicycler/raw/69e712eb95c4b9f8a46aade467260260a9ce7a91/sample_data/short_reads_2.fastq.gz
//   wget --quiet https://github.com/rrwick/Unicycler/raw/69e712eb95c4b9f8a46aade467260260a9ce7a91/sample_data/long_reads_low_depth.fastq.gz
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