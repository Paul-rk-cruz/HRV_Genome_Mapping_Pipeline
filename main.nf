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


def helpmsg() {

    log.info""""
    =========================================
     Virus Genome Mapping Pipeline :  v${version}
    =========================================
    Pipeline Usage:

    To run the pipeline, enter the following in the command line:

        nextflow run Virus_Genome_Mapping_Pipeline/main.nf --reads PATH_TO_FASTQ --viral_fasta .PATH_TO_VIR_FASTA --viral_index PATH_TO_VIR_INDEX --host_fasta PATH_TO_HOST_FASTA --host_index PATH_TO_HOST_INDEX --outdir ./output


    Arguments:
      --reads                       Path to input fastq.gz folder).
      --viral_fasta                 Path to fasta reference sequences (concatenated)
      --viral_index                 Path to indexed virus reference databases
      --host_fasta                  Path to host Fasta sequence
      --host_index                  Path to host fasta index

    Optional Commands:
      --singleEnd                   Specifies that the input is single end reads
      --notrim                      Specifying --notrim will skip the adapter trimming step.
      --saveTrimmed                 Save the trimmed Fastq files in the the Results directory.
      --trimmomatic_adapters_file   Adapters index for adapter removal
      --trimmomatic_mininum_length  Minimum length of reads
      --outdir                      The output directory where the results will be saved

    """.stripIndent()
}

/*
 * Configuration Setup
 */
params.help = false


// Pipeline version
version = '1.0'

// Show help msg
if (params.help){
    helpmsg()
    exit 0
}

