#!/usr/bin/env nextflow

/*
========================================================================================
                  Virus Genome Mapping Pipeline
========================================================================================
 #### Homepage / Documentation
 
 @#### Authors
Paul RK Cruz <kurtisc@uw.edu>
----------------------------------------------------------------------------------------
Pipeline overview:
 - 1. : Preprocessing
 	- 1.1: FastQC - sequence read quality control.
 	- 1.2: Trimmomatic - sequence trimming of adaptors and low quality reads.
 - 2. : Genome Mapping
 	- 2.1 : Bowtie2 - Remove host genome and map to reference Virus genome.
 	- 2.2 : Samtools - SAM and BAM file processing.
 - 3. : Variant calling, annotation and consensus:
  - 3.1 : Bcftools - Consensus in *.fasta format
 ----------------------------------------------------------------------------------------
*/


