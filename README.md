# Donut Falls

<img src="https://upload.wikimedia.org/wikipedia/commons/d/db/Donut_falls_in_utah.jpg" width="250" align="left" />

Named after the beautiful [Donut Falls](https://en.wikipedia.org/wiki/Doughnut_Falls)

Location: 40.630°N 111.655°W , 7,942 ft (2,421 m) elevation

More information about the trail leading up to this landmark can be found at [utah.com/hiking/donut-falls](https://utah.com/hiking/donut-falls)

Donut Falls is a [Nextflow](https://www.nextflow.io/) workflow developed by [@erinyoung](https://github.com/erinyoung) at the [Utah Public Health Laborotory](https://uphl.utah.gov/) for long-read [nanopore](https://nanoporetech.com) sequencing of microbial isolates. Built to work on linux-based operating systems. Additional config options are needed for cloud batch usage.

Donut Falls will also probably be a workflow of the [staphb-toolkit](https://github.com/StaPH-B/staphb_toolkit) once [@erinyoung]("https://github.com/erinyoung") gets around to it. Donut Falls, admittedly, is just a temporary stop-gap until nf-core's [genomeassembler](https://github.com/nf-core/genomeassembler) is released, so do not be suprised if this repository gets archived.

## Getting started

### Install dependencies
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
- [Singularity](https://singularity.lbl.gov/install-linux) or [Docker](https://docs.docker.com/get-docker/)

It is highly recommended to use [bandage](https://rrwick.github.io/Bandage/) to visualize the assemblies, but this is optional. 

### combine all fastq or fastq files into one file per barcode and rename file

Example with fastq.gz for barcode 01
```
cat fastq_pass/barcode01/*fastq.gz > reads/sample.fastq.gz
```

## Usage

```
nextflow run UPHL-BioNGS/Donut_Falls -profile singularity --reads <path to reads>
```

### If there are Illumina short - reads that can used for polishing

```
nextflow run UPHL-BioNGS/Donut_Falls -profile singularity --reads <path to reads> --illumina <path to illumina reads>
```

Illumina reads need to match the same naming convention as the nanopore reads (i.e. `12345.fastq.gz` for nanopore and `12345_R1.fastq.gz` and `12345_R2.fastq.gz` for Illumina)

### Switching assemblers
There are currently four options of assemblers
- [flye](https://github.com/fenderglass/Flye)
- [miniasm](https://github.com/rrwick/Minipolish)
- [raven](https://github.com/lbcb-sci/raven)
- [unicycler](https://github.com/rrwick/Unicycler) (requires short reads for hybrid assembly)

These are specified with the `assembler` paramater.

```
# flye (default)
nextflow run UPHL-BioNGS/Donut_Falls -profile singularity --reads reads --assembler flye
# miniasm and minipolish
nextflow run UPHL-BioNGS/Donut_Falls -profile singularity --reads reads --assembler miniasm
# raven
nextflow run UPHL-BioNGS/Donut_Falls -profile singularity --reads reads --assembler raven
# unicycler
nextflow run UPHL-BioNGS/Donut_Falls -profile singularity --reads reads --illumina illumina --assembler unicycler
```

### Reading the sequencing summary file
Although not used for anything else, the sequencing summary file can be read in and put through nanoplot to assess the quality of a sequencing run. This is an optional file and can be set with 'params.sequencing_summary'. 
```
nextflow run UPHL-BioNGS/Donut_Falls -profile singularity --reads <path to reads> --sequencing_summary <sequencing summary file>
```
* WARNING : Does not work with _older_ versions of the summary file.

## Diagram of workflow
<img src="./images/Donut Falls.png">


## Final file structure

<details>
   <summary>Final File Tree</summary>

```
donut_falls/
├── fastp
│   ├── sample_R1.fastq.gz
│   ├── sample_R2.fastq.gz
│   ├── fastp.html
│   └── fastp.json
├── filtlong
│   ├── sample_filtered.fastq
│   └── sample_filtered.fastq.gz
├── flye
│   ├── sample
│   │   ├── 00-assembly
│   │   │   ├── draft_assembly.fasta
│   │   │   └── draft_assembly.fasta.fai
│   │   ├── 10-consensus
│   │   │   ├── consensus.fasta
│   │   │   ├── minimap.bam.bai
│   │   │   └── minimap.stderr
│   │   ├── 20-repeat
│   │   │   ├── graph_after_rr.gv
│   │   │   ├── graph_before_rr.fasta
│   │   │   ├── graph_before_rr.gv
│   │   │   ├── read_alignment_dump
│   │   │   ├── repeat_graph_dump
│   │   │   └── repeat_graph_edges.fasta
│   │   ├── sample.fasta
│   │   ├── sample.fasta.fai
│   │   ├── sample.fasta.map-ont.mmi
│   │   ├── 30-contigger
│   │   │   ├── contigs.fasta
│   │   │   ├── contigs.fasta.fai
│   │   │   ├── contigs_stats.txt
│   │   │   ├── graph_final.fasta
│   │   │   ├── graph_final.gfa
│   │   │   ├── graph_final.gv
│   │   │   └── scaffolds_links.txt
│   │   ├── 40-polishing
│   │   │   ├── contigs_stats.txt
│   │   │   ├── edges_aln.bam.bai
│   │   │   ├── filtered_contigs.fasta
│   │   │   ├── filtered_contigs.fasta.fai
│   │   │   ├── filtered_stats.txt
│   │   │   ├── minimap_1.bam.bai
│   │   │   ├── minimap.stderr
│   │   │   └── polished_edges.gfa
│   │   ├── assembly.fasta
│   │   ├── assembly_graph.gfa
│   │   ├── assembly_graph.gv
│   │   ├── assembly_info.txt
│   │   ├── flye.log
│   │   ├── graph.png
│   │   └── params.json
├── logs
├── medaka
│   └── sample
│       ├── sample_medaka_consensus.fasta
│       ├── calls_to_draft.bam
│       ├── calls_to_draft.bam.bai
│       ├── consensus.fasta
│       ├── consensus.fasta.gaps_in_draft_coords.bed
│       └── consensus_probs.hdf
├── miniasm
│   └── sample
│       ├── sample.fasta
│       └── sample.gfa
├── nanoplot
│   ├── sample_il
│   │   ├── Dynamic_Histogram_Read_length.html
│   │   ├── Dynamic_Histogram_Read_length.png
│   │   ├── HistogramReadlength.png
│   │   ├── LengthvsQualityScatterPlot_dot.png
│   │   ├── LengthvsQualityScatterPlot_kde.png
│   │   ├── LogTransformed_HistogramReadlength.png
│   │   ├── NanoPlot_20220728_1344.log
│   │   ├── NanoPlot-report.html
│   │   ├── NanoStats.txt
│   │   ├── Weighted_HistogramReadlength.png
│   │   ├── Weighted_LogTransformed_HistogramReadlength.png
│   │   └── Yield_By_Length.png
│   ├── sample
│   │   ├── 2735071_NanoStats.csv
│   │   ├── 2735071_NanoStats.txt
│   │   ├── Dynamic_Histogram_Read_length.html
│   │   ├── Dynamic_Histogram_Read_length.png
│   │   ├── HistogramReadlength.png
│   │   ├── LengthvsQualityScatterPlot_dot.png
│   │   ├── LengthvsQualityScatterPlot_kde.png
│   │   ├── LogTransformed_HistogramReadlength.png
│   │   ├── NanoPlot_20220728_1223.log
│   │   ├── NanoPlot-report.html
│   │   ├── NanoStats.txt
│   │   ├── Weighted_HistogramReadlength.png
│   │   ├── Weighted_LogTransformed_HistogramReadlength.png
│   │   └── Yield_By_Length.png
│   ├── NanoStats.csv
│   └── summary
│       ├── NanoPlot_20220728_1344.log
│       ├── NanoStats_barcoded.txt
│       ├── NanoStats.txt
├── polca
│   ├── sample
│   │   ├── round_1
│   │   │   ├── sample_medaka_consensus.fasta.alignSorted.bam
│   │   │   ├── sample_medaka_consensus.fasta.alignSorted.bam.bai
│   │   │   ├── sample_medaka_consensus.fasta.batches
│   │   │   ├── sample_medaka_consensus.fasta.bwa.amb
│   │   │   ├── sample_medaka_consensus.fasta.bwa.ann
│   │   │   ├── sample_medaka_consensus.fasta.bwa.bwt
│   │   │   ├── sample_medaka_consensus.fasta.bwa.pac
│   │   │   ├── sample_medaka_consensus.fasta.bwa.sa
│   │   │   ├── sample_medaka_consensus.fasta.fai
│   │   │   ├── sample_medaka_consensus.fasta.fix.success
│   │   │   ├── sample_medaka_consensus.fasta.index.success
│   │   │   ├── sample_medaka_consensus.fasta.map.success
│   │   │   ├── sample_medaka_consensus.fasta.names
│   │   │   ├── sample_medaka_consensus.fasta.PolcaCorrected.fa
│   │   │   ├── sample_medaka_consensus.fasta.report
│   │   │   ├── sample_medaka_consensus.fasta.report.success
│   │   │   ├── sample_medaka_consensus.fasta.sort.success
│   │   │   ├── sample_medaka_consensus.fasta.unSorted.sam
│   │   │   ├── sample_medaka_consensus.fasta.vcf
│   │   │   └── sample_medaka_consensus.fasta.vc.success
│   │   └── round_{2,3,4,5}
│   │       ├── sample_round1.fa.alignSorted.bam
│   │       ├── sample_round1.fa.alignSorted.bam.bai
│   │       ├── sample_round1.fa.batches
│   │       ├── sample_round1.fa.bwa.amb
│   │       ├── sample_round1.fa.bwa.ann
│   │       ├── sample_round1.fa.bwa.bwt
│   │       ├── sample_round1.fa.bwa.pac
│   │       ├── sample_round1.fa.bwa.sa
│   │       ├── sample_round1.fa.fai
│   │       ├── sample_round1.fa.fix.success
│   │       ├── sample_round1.fa.index.success
│   │       ├── sample_round1.fa.map.success
│   │       ├── sample_round1.fa.names
│   │       ├── sample_round1.fa.PolcaCorrected.fa
│   │       ├── sample_round1.fa.report
│   │       ├── sample_round1.fa.report.success
│   │       ├── sample_round1.fa.sort.success
│   │       ├── sample_round1.fa.unSorted.sam
│   │       ├── sample_round1.fa.vcf
│   │       └── sample_round1.fa.vc.success
│   └── sample_final.fa
└── raven
    ├── 2735071
    │   ├── 2735071.fasta
    │   ├── 2735071.gfa
    │   └── graph.png
    ├── 2953647
    │   ├── 2953647.fasta
    │   ├── 2953647.gfa
    │   └── graph.png
    ├── 3302496
    │   ├── 3302496.fasta
    │   ├── 3302496.gfa
    │   └── graph.png
    └── 3313737
        ├── 3313737.fasta
        ├── 3313737.gfa
        └── graph.png    
```
</details>

## Credits

Donut Falls would not be possible without
- [flye](https://github.com/fenderglass/Flye)
- [raven](https://github.com/lbcb-sci/raven)
- [miniasm and minipolish](https://github.com/rrwick/Minipolish)
- [filtlong](https://github.com/rrwick/Filtlong)
- [fastp](https://github.com/OpenGene/fastp)
- [nanoplot](https://github.com/wdecoster/NanoPlot)
- [polca](https://github.com/alekseyzimin/masurca)
- [medaka](https://github.com/nanoporetech/medaka)
- [any2fasta](https://github.com/tseemann/any2fasta)
- [bgzip](https://github.com/samtools/htslib)
- [unicycler](https://github.com/rrwick/Unicycler)
- [circlator](https://github.com/sanger-pathogens/circlator)
- [trycycler](https://github.com/rrwick/Trycycler)

## Frequently Asked Questions (aka FAQ)
### What do I do if I encounter an error?

**TELL ME ABOUT IT!!!**
* [Github issue](https://github.com/UPHL-BioNGS/Donut_Falls/issues)
* [Email me](eriny@utah.gov)
* Send me a message on slack

Be sure to include the command used, what config file was used, and what the **nextflow** error was. 

### Which assembler is best?
I am flattered that there are those out there that think I can summarize a response into something that would fit into a github readme. In summary, the three aligners do different things and have different issues. The default is flye due to its popularity.

I took four samples through each of the three assemblers and visualized their gfa file in bandage. For the most part, these all look very similar. There are instances, however, where the genome is closed from one assembler but could not be circularized in another. Sometimes there are different numbers of plasmids as well. 

| | flye | miniasm and minipolish | raven |
|:-------------------------:|:-------------------------:|:-------------------------:|:-------------------------:|
| sample 1 | <img width="1604" src="./images/flye_2735071.png"> | <img width="1604" src="./images/miniasm_2735071.png"> | <img width="1604" src="./images/raven_2735071.png"> |
| sample 2 | <img width="1604" src="./images/flye_2953647.png"> | <img width="1604" src="./images/miniasm_2953647.png"> | <img width="1604" src="./images/raven_2953647.png"> |
| sample 3 | <img width="1604" src="./images/flye_3302496.png"> | <img width="1604" src="./images/miniasm_3302496.png"> | <img width="1604" src="./images/raven_3302496.png"> |
| sample 4 | <img width="1604" src="./images/flye_3313737.png"> | <img width="1604" src="./images/miniasm_3313737.png"> | <img width="1604" src="./images/raven_3313737.png"> |

In the future, it would be nice to have something [Trycycler](https://github.com/rrwick/Trycycler)-esq to include that would resolve or provide confidence in assembler choices.

### Is there anything for quality evalution?
Currenlty the brunt of QC is being undertaken by [NanoPlot](https://github.com/wdecoster/NanoPlot). 

Graphs for the sequencing run are generated by NanoPlot like these:
<img src="./images/barcode01_ActivityMap_ReadsPerChannel.png" width="1604" align="left" />

And read quality is visualized for both long nanopore reads and short illumina reads.
| Illumina "short" reads | Nanopore "long" reads|
|:-------------------------:|:-------------------------:|
| <img width="1604" src="./images/ILL_LVQ.png"> | <img width="1604" src="./images/ONT_LVQ.png"> |

In general, longer reads are better.


### Where is an example config file?
To get a copy of the [template file](.configs/donut_falls_template.config) that Donut Falls uses by default, run 
```
nextflow run UPHL-BioNGS/Donut_Falls --config_file true
```
This creates an `edit_me.config` file in the current directory that the **End User** can edit for their own purpose. This file can be renamed with no penalty.

To use this edited config file, simply use `-c` on the command line.

```
nextflow run UPHL-BioNGS/Donut_Falls -profile singularity --reads reads -c edit_me.config
```

### Why does flye keep failing?
Most of this has to do with the quality of the nanpore reads and flye's internal workings. Right now, flye errors are set to be ignored by default, but this can be changed by the **End User** in a config file. A common error and ways around it are found in [this issue thread](https://github.com/fenderglass/Flye/issues/128). This means that the **End User** will need a config file with the specified flye parameters.

For example
```
params.flye_options = "--asm-coverage 50"
# OR
params.flye_options = "--meta"
```

### What is currently in development for this workflow?
Hopefully these ideas will be added _soon_, but due to limited time and unlimited responsibilites, these have been put on the backburner.

- Circlator(https://github.com/sanger-pathogens/circlator) to rotate complete plasmids
  - Was not included in the current version because there is not an easy way in the fasta to indicate if the sequence is closed or not.
- Something [Trycycler](https://github.com/rrwick/Trycycler)-esq
  - Currently there are three dsl1 workflows that need to be converted to dsl-2 in [this repo](./trycycler_dsl1). There is still a mandatory manual curation step inherent to Trycycler.
- Some QC on assembled sequences
  - I just need some time some tools and build some containers.

### How were raven, flye, and miniasm chosen for this workflow?
They perform _well_, their containers were easy to create, and they were all part of the [Trycycler walkthrough](https://github.com/rrwick/Trycycler/wiki/Clustering-contigs). 

I did attempt adding [canu](https://github.com/marbl/canu), but the assembly took forever for my tests. 

If the **End User** prefers other assemblers, please let me know and we'll work in some options. 

**Warning** : If there's not a relaible container of the suggested tool, I'll request the **End User** create a container for that tool and contribute to [StaPH-B's docker repositories](https://github.com/StaPH-B/docker-builds).

### What if I'm starting with fast5 files?
Then the **End User** needs to do basecalling and demultiplexing with [guppy](https://nanoporetech.com/nanopore-sequencing-data-analysis) first. 

Something like the following to get the fastq files.
```
# With config file
guppy_basecaller -i <input path> -s <save path> -c <config file> [options]
# With flowcell and kit name
guppy_basecaller -i <input path> -s <save path> --flowcell <flowcell name> --kit <kit name>
```

Something like the following to get the demultiplexed fastq files.
```
guppy_barcoder -i <input fastq path> -s <save path> --barcode_kits <kit name>
```
