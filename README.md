# Donut Falls

# Warning : What you currently see is under development. It is likely incomplete and incorrect. 

### (We're working on it.)

<img src="https://upload.wikimedia.org/wikipedia/commons/d/db/Donut_falls_in_utah.jpg" width="250" align="left" />

Named after the beautiful [Donut Falls](https://en.wikipedia.org/wiki/Doughnut_Falls)

Location: 40.630°N 111.655°W , 7,942 ft (2,421 m) elevation

More information about the trail leading up to this landmark can be found at [utah.com/hiking/donut-falls](https://utah.com/hiking/donut-falls)

Donut Falls is a [Nextflow](https://www.nextflow.io/) workflow developed by [@erinyoung](https://github.com/erinyoung) at the [Utah Public Health Laborotory](https://uphl.utah.gov/) for long-read [nanopore](https://nanoporetech.com) sequencing of microbial isolates. Built to work on linux-based operating systems. Additional config options are needed for cloud batch usage.

Donut Falls will also probably be a workflow of the [staphb-toolkit](https://github.com/StaPH-B/staphb_toolkit) once [@erinyoung]("https://github.com/erinyoung") gets around to it and all the containers are ready.

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
nextflow run UPHL-BioNGS/Donut_falls -profile singularity --reads <path to reads> --reads reads
```

## If there are Illumina short - reads that can used for polishing

```
nextflow run UPHL-BioNGS/Donut_falls -profile singularity --reads <path to reads> --reads reads --illumina
```

Illumina reads need to match the same naming convention as the nanopore reads (i.e. 12345.fastq.gz for nanopore 12345_R1.fastq.gz and 12345_R2.fastq.gz for Illumina)


# Frequently Asked Questions (aka FAQ)
### What do I do if I encounter an error?

**TELL ME ABOUT IT!!!**
* [Github issue](https://github.com/UPHL-BioNGS/Donut_Falls/issues)
* [Email me](eriny@utah.gov)
* Send me a message on slack

Be sure to include the command used, what config file was used, and what the **nextflow** error was. 

### Where is an example config file?
There is a template file with all the variables [here](.configs/donut_falls_template.config) that the **End User** can copy and edit.

There's also a config file what we use here at UPHL [here](./configs/UPHL.config).

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

![alt text](https://uphl.utah.gov/wp-content/uploads/New-UPHL-Logo.png)
