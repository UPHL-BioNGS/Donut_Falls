# Donut Falls

# Warning : What you currently see is under development. It is likely incomplete and incorrect. 

### (We're working on it.)

Named after the beautiful [Donut Falls](https://en.wikipedia.org/wiki/Doughnut_Falls)

Location: 40.630°N 111.655°W , 7,942 ft (2,421 m) elevation

![alt text](https://upload.wikimedia.org/wikipedia/commons/d/db/Donut_falls_in_utah.jpg)

Donut Falls is a [Nextflow](https://www.nextflow.io/) workflow developed by [@erinyoung]("https://github.com/erinyoung") at the [Utah Public Health Laborotory](https://uphl.utah.gov/) for long-read [nanopore](https://nanoporetech.com) sequencing of microbial isolates. Built to work on linux-based operating systems. Additional config options are needed for cloud batch usage.

The heavy lifter of this workflow is [trycycler](https://github.com/rrwick/Trycycler), and many of the processes are **heavily inspired** by the [trycycler wiki](https://github.com/rrwick/Trycycler/wiki).

Donut Falls will also probably be a workflow of the [staphb-toolkit](https://github.com/StaPH-B/staphb_toolkit) once [@erinyoung]("https://github.com/erinyoung") gets around to it.

# Getting started

```
git clone https://github.com/UPHL-BioNGS/Donut_Falls.git
```

To make your life easier, follow with

```
cd Donut_Falls
git init
```

so that you can use `git pull` for updates.

## Prior to starting the workflow

### Install dependencies
- [Nextflow](https://www.nextflow.io/docs/latest/getstarted.html)
   - Nextflow version 20+ is required (`nextflow -v` to check your installation)
- [Singularity](https://singularity.lbl.gov/install-linux)

or
- [Docker](https://docs.docker.com/get-docker/) (*with the caveat that the creator and maintainer uses singularity and may not be able to troubleshoot all docker issues*)

# Usage
```
nextflow run Donut_Falls.nf -c configs/singularity.config
```

# Usage including guppy basecalling

```
TODO
```

# Usage for if you have Illumina short-read sequencing
```
nextflow run Donus_Fall.nf -c configs/singularity.config --short-reads fastq
```

## The main components of Donut Falls are:

- [trycycler](https://github.com/rrwick/Trycycler)
- [guppy](https://nanoporetech.com/nanopore-sequencing-data-analysis) - for basecalling, debarcoding, and adapter trimming
- [filtlong](https://github.com/rrwick/Filtlong) - for filtering reads
- [flye](https://github.com/fenderglass/Flye) - de novo assembly of long reads
- [medaka](https://andersen-lab.github.io/ivar/html/manualpage.html) - calling variants and creating a consensus fasta; optional primer trimmer
- [pilon](http://www.htslib.org/) - for QC metrics and sorting; optional primer trimmer; optional converting bam to fastq files
- [unicycler](https://github.com/s-andrews/FastQC) - for hybrid assembly
- [quast](https://bedtools.readthedocs.io/en/latest/) - for depth estimation over amplicons

# Frequently Asked Questions (aka FAQ)
### What do I do if I encounter an error?

**TELL ME ABOUT IT!!!**
* [Github issue](https://github.com/UPHL-BioNGS/Donut_Falls/issues)
* Email me
* Send me a message on slack

Be sure to include the command that you used, what config file you used, and what the **nextflow** error was. 

### Where is an example config file?
You're more than welcome to look at what we use here at UPHL [here](./configs/UPHL.config).

### This workflow has too many bells and whistles. I really only care about generating a consensus fasta. How do I get rid of all the extras?

Change the parameters in your config file and set most of them to false. 

```
```
And, yes, this means I added some bells and whistles so you could turn off the bells and whistles. /irony

# Directed Acyclic Diagrams (DAG)
### Full workflow
![alt text](images/Donut_falls_workflow.png)

![alt text](https://uphl.utah.gov/wp-content/uploads/New-UPHL-Logo.png)
