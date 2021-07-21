# Donut Falls

# Warning : What you currently see is under development. It is likely incomplete and incorrect. 

### (We're working on it.)

<img src="https://upload.wikimedia.org/wikipedia/commons/d/db/Donut_falls_in_utah.jpg" width="250" align="left" />

Named after the beautiful [Donut Falls](https://en.wikipedia.org/wiki/Doughnut_Falls)

Location: 40.630°N 111.655°W , 7,942 ft (2,421 m) elevation

More information about the trail leading up to this landmark can be found at [utah.com/hiking/donut-falls](https://utah.com/hiking/donut-falls)

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

## The setup:

The default file structure for the MinKnow on a GridIon/MinIon with basecalling and demultiplexing is `<directory>/fastq_pass`. This workflow assumes **the end user** to be in that directory, but this can also be specified with `'params.fastq_pass'`.

<details>
   <summary>Initial File Tree with 12 barcodes</summary>

```
├── sample_key.csv    # optional file that connects a barcode to additional identifiers and (even more optional) paired-end Illumina reads
├── fast5_fail
│   ├── barcode01
│   ├── barcode02
│   ├── barcode03
│   ├── barcode04
│   ├── barcode05
│   ├── barcode06
│   ├── barcode07
│   ├── barcode08
│   ├── barcode09
│   ├── barcode10
│   ├── barcode11
│   ├── barcode12
│   └── unclassified
├── fast5_pass
│   ├── barcode01
│   ├── barcode02
│   ├── barcode03
│   ├── barcode04
│   ├── barcode05
│   ├── barcode06
│   ├── barcode07
│   ├── barcode08
│   ├── barcode09
│   ├── barcode10
│   ├── barcode11
│   ├── barcode12
│   └── unclassified
├── fastq_fail
│   ├── barcode01
│   ├── barcode02
│   ├── barcode03
│   ├── barcode04
│   ├── barcode05
│   ├── barcode06
│   ├── barcode07
│   ├── barcode08
│   ├── barcode09
│   ├── barcode10
│   ├── barcode11
│   ├── barcode12
│   └── unclassified
├── fastq_pass
│   ├── barcode01
│   ├── barcode02
│   ├── barcode03
│   ├── barcode04
│   ├── barcode05
│   ├── barcode06
│   ├── barcode07
│   ├── barcode08
│   ├── barcode09
│   ├── barcode10
│   ├── barcode11
│   ├── barcode12
│   └── unclassified
├── illumina_fastq
├── other_reports
```

</details>

# Usage for `Phase 1`

```
nextflow run Donut_Falls.nf -c configs/singularity.config
```

This workflow is meant to run in **TWO** phases. The first phase takes the multiple fastq files from the sequencing run and combines them into one file. If a `'sample key'` is supplied, the file will be renamed (more details on the sample key below). This fastq file is filtered with [filtlong](https://github.com/rrwick/Filtlong), split into groups by [trycycler](https://github.com/rrwick/Trycycler), and undergoes *de novo* alignment with three different aligners : [canu](https://github.com/marbl/canu), [raven](https://github.com/lbcb-sci/raven), and [flye](https://github.com/fenderglass/Flye). Once the different assemblies have been completed, [trycycler](https://github.com/rrwick/Trycycler) clusters the resultant contigs. 

There is then a break in the workflow that requires manual input from the **End User**. The **End User** should look at the newick files in `donut_falls/trycycler_cluster/sample/newick`. This can be done multiple ways, but here I will suggest dragging and dropping the file into [ETE's Phylogenetic tree (newick) viewer](http://etetoolkit.org/treeview/). The **End User** must remove any clusters or contigs that do not appear to *cluster well*. This can be done by acutally removing the files and directory with `rm -rf` (or whatever tool your prefer) or they can just be moved out the way into a different directory. 

#### Clusters can be moved
```
mkdir -p donut_falls/trycycler_cluster/sample/failed_clusters
mv donut_falls/trycycler_cluster/sample/cluster_002 donut_falls/trycycler_cluster/sample/failed_clusters/cluster_002
```
#### Or individual files can be moved
```
mkdir -p donut_falls/trycycler_cluster/sample/failed_clusters/cluster_002/1_contigs
mv donut_falls/trycycler_cluster/sample/cluster_002/1_contigs/A_contig_3.fasta donut_falls/trycycler_cluster/sample/failed_clusters/cluster_002/1_contigs/A_contig_3.fasta
```

As for what is meant by *"clustering well"*, please read Trycycler's wiki at [https://github.com/rrwick/Trycycler/wiki/Clustering-contigs](https://github.com/rrwick/Trycycler/wiki/Clustering-contigs#choose-your-clusters) for a good definition and some great examples.


# Usage for `Phase 2`
Once the **End User** has manually inspected the contig clustering and removed abberant contigs and clusters, `phase 2` can begin. The initial directories for the filtered fastq files in `phase 1` are still required by the workflow, so the `-resume` paramater **must** now be included in the command or the workflow will start from the beginning (which is really annoying). The `params.phase2` must also be updated to `true` in a config file (`'params.phase2 = true'`) or on the command line as follows:

```
nextflow run Donut_Falls.nf -c configs/singularity.config -resume --phase2 true
```

The second phase of the workflow takes the clusters that were generated in `phase 1` and [trycycler](https://github.com/rrwick/Trycycler) will reconcile differences between the assemblies, align them, parition the reads, and generate a consensus fasta. This fasta is then polished with [medaka](https://github.com/nanoporetech/medaka). [Trycycler](https://github.com/rrwick/Trycycler) will fail at the reconcile step if there is a contig that is causing problems. A dotplot at `donut_falls/trycycler_dotplot/sample/` should help you identify which contigs are creating issues. More information can be found in [trycycler's wiki](https://github.com/rrwick/Trycycler/wiki/Reconciling-contigs#dotplots). The problematic contigs can then be removed and the workflow can be resumed just like before.

Optionally, `phase 2` can include *an extra polishing "step"* with Illumina paired-end short-read sequencing via [pilon](https://github.com/broadinstitute/pilon). This requires the `sample key` csv file (specify with `'params.sample_key'`), an associated parameter must be changed to true, `'params.illumina' = true`, and the directory with paired-end fastq files must be specified with `'params.illumina_fastq'` in a config file or on the command line as follows:

# Usage for `Phase 2` if you have Illumina short-read sequencing for polishing
```
nextflow run Donut_Falls.nf -c configs/singularity.config -resume --phase2 true --illumina true --illumina_fastq illumina_fastq
```
This will add a quick QC step for the Illumina reads with [fastp](https://github.com/OpenGene/fastp), alignment with [bwa](http://bio-bwa.sourceforge.net/), sorting with [samtools](http://www.htslib.org/), and polishing with [pilon](https://github.com/broadinstitute/pilon) until no changes are being made during the polishing step.

# Examining your sequences
In a perfect world, you now have a consensus sequence of complete genome of whatever you were sequencing. To check how well this was done, I recommend looking at your `X` file with [bandage](https://rrwick.github.io/Bandage/). This fasta sequence can now undergo any other bioinformatic process that a fasta file could go through: 
- Gene annotation with [Prokka](https://github.com/tseemann/prokka) or [PGAP](https://www.ncbi.nlm.nih.gov/genome/annotation_prok/) (and submission to NCBI!)
- Antimicrobial gene prediction with [AMRFinderPlus](https://github.com/ncbi/amr)
- Serotyping with tools with [SeqSero2](https://github.com/denglab/SeqSero2)
- And More!

# Creating the sample key
The `sample key` a csv file with the barcode in the first column, and the id in the second column. Common bioinformatic constraints apply. Do not include weird special characters or spaces in your ids ("-", "_", and "." are okay), and the ids must be unique. The `sample key` file is set with `'params.sample_key'`.

Example:
```
$ cat sample_key.csv 
#barcode,id
barcode01,2820899
barcode02,2788440
```

If you have Illumina paired-end files, the third column becomes the forward reads, and the fourth column is the reverse reads. The Illumina to nanopore sequences MUST be a 1:1 ratio. If your paired-end files can be used for more than one sample on your nanopore run, create a copy of your Illumina fastq reads with a slightly different name.

Example:
```
$ cat sample_key.csv 
#barcode,id,R1,R2
barcode01,2820899,2820899_R1.fastq.gz,2820899_R2.fastq.gz
barcode02,2788440,2788440_R1.fastq.gz,2788440_R2.fastq.gz
barcode03,2788440A,2788440A_R1.fastq.gz,2788440A_R2.fastq.gz
```

## Donut Falls wouldn't be possible without:

- [trycycler](https://github.com/rrwick/Trycycler)
- [filtlong](https://github.com/rrwick/Filtlong)     - for filtering reads
- [flye](https://github.com/fenderglass/Flye)        - de novo assembly of long reads
- [canu](https://github.com/marbl/canu)              - de novo assembly of long reads
- [raven](https://github.com/lbcb-sci/raven)         - de novo assembly of long reads
- [medaka](https://github.com/nanoporetech/medaka)   - polisher with nanopore reads
- [pilon](http://www.htslib.org/)                    - polisher with Illumina reads
- [fastp](https://github.com/OpenGene/fastp)         - quick QC of Illumina reads
- [bwa](http://bio-bwa.sourceforge.net/)             - alignment of Illumina reads
- [samtools](http://www.htslib.org/)                 - sorts bam files generated by bwa

# Frequently Asked Questions (aka FAQ)
### What do I do if I encounter an error?

**TELL ME ABOUT IT!!!**
* [Github issue](https://github.com/UPHL-BioNGS/Donut_Falls/issues)
* [Email me](eriny@utah.gov)
* Send me a message on slack

Be sure to include the command that you used, what config file you used, and what the **nextflow** error was. 

### Where is an example config file?
There is a template file with all the variables [here](.configs/donut_falls_template.config) that you can copy and edit for your purposes.

You're also more than welcome to look at what we use here at UPHL [here](./configs/UPHL.config).

### How were raven, flye, and canu chosen for this workflow?
They perform **well** and their containers were easy to create.

If you prefer other assemblers, please let me know and we'll work in some options. Warning : If there's not a relaible container of the tool that you suggest, I'll request that you create a container for that tool and contribute to [StaPH-B's docker repositories](https://github.com/StaPH-B/docker-builds).

### What if I'm starting with fast5 files?
Then you need to do basecalling and demultiplexing with [guppy](https://nanoporetech.com/nanopore-sequencing-data-analysis) first. 

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

### What if I'm starting with a different file structure?
I do have a goal to make this more flexible as time goes on. For right now, you'll need to adjust your file structure so that `'params.fastq_pass'` is a directory with `'barcode*/*fastq*'` files. 

```
├── directory
    ├── barcode01
    │    └── *.fastq or *.fastq.gz
    ├── barcode02
    │    └── *.fastq or *.fastq.gz
```

# Directed Acyclic Diagrams (DAG)
### Full workflow
![alt text](images/Donut_falls_workflow.png)

![alt text](https://uphl.utah.gov/wp-content/uploads/New-UPHL-Logo.png)
