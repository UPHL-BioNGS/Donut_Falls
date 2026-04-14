---
title: 'Donut Falls: A Nextflow workflow for automated de novo assembly and polishing of microbial long-read sequencing data'
tags:
  - bioinformatics
  - nanopore
  - nextflow
  - microbial genomics
  - assembly
authors:
  - name: Erin Young
    orcid: 0000-0002-7535-006X
    affiliation: 1
affiliations:
  - name: Utah Public Health Laboratory, Department of Health and Human Services, State of Utah, Salt Lake City, UT, USA
    index: 1
date: 14 April 2026
bibliography: paper.bib
---

# Summary
Donut Falls is a Nextflow-based pipeline designed for the assembly of microbial isolates using Oxford Nanopore Technologies (ONT) long reads. It automates the process of read QC, subsampling, de novo assembly, rotation, and polishing.

# Statement of Need
While several assemblers exist for ONT data (e.g., Flye, Raven), a standardized, reproducible workflow is required in public health settings to ensure consistent results across different laboratories. Donut Falls addresses this by providing a multi-assembler framework that supports both long-read-only and hybrid assembly (using Illumina reads) while performing rigorous post-assembly quality assessment with BUSCO and CirculoCov.

# Workflow Architecture
[Insert a brief description of the logic: Fastp -> Rasusa -> Assembler -> Dnaapler -> Polishing -> MultiQC]

# Acknowledgements
We thank the Utah Public Health Laboratory BioNGS team...

# References
