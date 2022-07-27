# Donut Falls

<img src="https://upload.wikimedia.org/wikipedia/commons/d/db/Donut_falls_in_utah.jpg" width="250" align="left" />

Named after the beautiful [Donut Falls](https://en.wikipedia.org/wiki/Doughnut_Falls)

Location: 40.630°N 111.655°W , 7,942 ft (2,421 m) elevation

More information about the trail leading up to this landmark can be found at [utah.com/hiking/donut-falls](https://utah.com/hiking/donut-falls)

Donut Falls is a [Nextflow](https://www.nextflow.io/) workflow developed by [@erinyoung](https://github.com/erinyoung) at the [Utah Public Health Laborotory](https://uphl.utah.gov/) for long-read [nanopore](https://nanoporetech.com) sequencing of microbial isolates. Built to work on linux-based operating systems. Additional config options are needed for cloud batch usage.

Donut Falls will also probably be a workflow of the [staphb-toolkit](https://github.com/StaPH-B/staphb_toolkit) once [@erinyoung]("https://github.com/erinyoung") gets around to it. Donut Falls, admittedly, is just a temporary stop-gap until nf-core's [genomeassembler](https://github.com/nf-core/genomeassembler) is released, so do not be suprised if this repository gets archived.

# Getting started

## Install dependencies
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
   - Nextflow version 20+ is required (`nextflow -v` to check)
- [Singularity](https://singularity.lbl.gov/install-linux) or [Docker](https://docs.docker.com/get-docker/)

It is highly recommended to use [bandage](https://rrwick.github.io/Bandage/) to visualize the assemblies, but this is optional. 

## combine all fastq or fastq files into one file per barcode and rename file

Example with fastq.gz for barcode 01
```
cat fastq_pass/barcode01/*fastq.gz > reads/sample.fastq.gz
```

# Usage

```
nextflow run UPHL-BioNGS/Donut_Falls -profile singularity --reads <path to reads>
```

## If there are Illumina short - reads that can used for polishing

```
nextflow run UPHL-BioNGS/Donut_Falls -profile singularity --reads <path to reads> --illumina <path to illumina reads>
```

Illumina reads need to match the same naming convention as the nanopore reads (i.e. `12345.fastq.gz` for nanopore and `12345_R1.fastq.gz` and `12345_R2.fastq.gz` for Illumina)

## Switching assemblers
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

This workflow wouldn't be possible without
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

# Frequently Asked Questions (aka FAQ)
### What do I do if I encounter an error?

**TELL ME ABOUT IT!!!**
* [Github issue](https://github.com/UPHL-BioNGS/Donut_Falls/issues)
* [Email me](eriny@utah.gov)
* Send me a message on slack

Be sure to include the command used, what config file was used, and what the **nextflow** error was. 

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
