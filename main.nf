#!/usr/bin/env nextflow

/*
========================================================================================
                  Virus Genome Mapping Pipeline
========================================================================================
 Github Repo:
 Greninger Lab
 
 Author:
 Paul RK Cruz <kurtisc@uw.edu>
----------------------------------------------------------------------------------------
Pipeline overview:
 - 1. : Fastq File Processing
 		-Trimmomatic - sequence trimming of adaptors and low quality reads.
 - 2. : Genome Mapping
 		-Bowtie2 - Remove host genome and map to reference Virus genome.
 		-Samtools - SAM and BAM file processing.
 - 3. : Variant calling, annotation and consensus:
  		-Bcftools - Consensus in *.fasta format

  Dependencies:
  
  trimmomatic
  bowtie2
  samtools
  bedtools
  Bcftools

	PATHS FOR EASTLAKE KC iMac - TESTING & DEBUGGING
	PATH to Virus Reference Fasta:
	/Users/Kurtisc/Downloads/CURRENT/Virus_Genome_Mapping_Pipeline/Virus_Genome_Mapping_Pipeline/virus_ref_db/rhv_abc_sars2.fasta 
	
	PATH to indexed Reference Fasta Database Files:
	/Users/Kurtisc/Downloads/CURRENT/Virus_Genome_Mapping_Pipeline/Virus_Genome_Mapping_Pipeline/virus_ref_db

 ----------------------------------------------------------------------------------------
*/


def helpmsg() {

    log.info""""
    ____________________________________________
     Virus Genome Mapping Pipeline :  v${version}
    ____________________________________________
    
	Pipeline Usage:

    To run the pipeline, enter the following in the command line:

        nextflow run Virus_Genome_Mapping_Pipeline/main.nf --reads PATH_TO_FASTQ --virus_fasta .PATH_TO_VIR_FASTA --virus_index PATH_TO_VIR_INDEX --outdir ./output


    Valid CLI Arguments:
      --reads                       Path to input fastq.gz folder).
      --virus_fasta                 Path to fasta reference sequences (concatenated)
      --virus_index                 Path to indexed virus reference databases
      --singleEnd                   Specifies that the input fastq files are single end reads
      --notrim                      Specifying --notrim will skip the adapter trimming step
      --saveTrimmed                 Save the trimmed Fastq files in the the Results directory
      --trimmomatic_adapters_file   Adapters index for adapter removal
      --trimmomatic_mininum_length  Minimum length of reads
	  --withFastQC					Runs a quality control check on fastq files
      --outdir                      The output directory where the results will be saved
	  --helpmsg						Displays help message in terminal

    """.stripIndent()
}
// Check Nextflow version
nextflow_req_v = '20.10.0'
try {
    if( ! nextflow.version.matches(">= $nextflow_req_v") ){
        throw GroovyException('> ERROR: The version of Nextflow running on your machine is out dated.\n>Please update to Version '$nextflow_req_v)
    }
} catch (all) {
	log.error"ERROR: This version of Nextflow is out of date.\nPlease update to the latest version of Nextflow."
}
/*
 * Configuration Setup
 */
params.help = True

// Pipeline version
version = '1.0'

// Show help msg
if (params.help){
    helpmsg()
    exit 0
}
// Check for virus genome reference indexes
params.virus_fasta = false
if( params.virus_fasta ){
    virus_fasta_file = file(params.virus_fasta)
    if( !virus_fasta_file.exists() ) exit 1, "> Virus fasta file not found: ${params.virus_fasta}.\n> Please specify a valid file path!"
}
// Channel for input fastq files
Channel
    .fromFilePairs( params.reads, size: params.singleEnd ? 1 : 2 )
    .ifEmpty { exit 1, "> Invalid sequence read type.\n> Please retry with --singleEnd" }
    .into { raw_reads_fastqc; raw_reads_trimming }

if( params.virus_index ){
// Channel for virus genome reference indexes
	Channel
        .fromPath(params.virus_index)
        .ifEmpty { exit 1, "> Error: Virus index not found: ${params.virus_index}.\n> Please specify a valid file path!"}
        .into { virus_index_files; virus_index_files_ivar; virus_index_files_variant_calling }
}
// Check for fastq
params.reads = false
if (! params.reads ) exit 1, "> Error: Fastq files not found: $params.reads. Please specify a valid path with --reads"
// Single-end read option
params.singleEnd = false
// Trimming parameters
params.notrim = false
// Output files options
params.saveTrimmed = false
// Default trimming options
params.trimmomatic_adapters_file = "\$TRIMMOMATIC_PATH/adapters/TruSeq3-PE.fa"
params.trimmomatic_adapters_parameters = "2:30:10"
params.trimmomatic_window_length = "4"
params.trimmomatic_window_value = "20"
params.trimmomatic_mininum_length = "50"

// log files header
log.info "____________________________________________"
log.info " Virus Genome Mapping Pipeline :  v${version}"
log.info "____________________________________________"
def summary = [:]
summary['Fastq Files:']               = params.reads
summary['Read type:']           = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Virus Reference:']           = params.virus_fasta
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current directory path:']        = "$PWD"
summary['Working directory path:']         = workflow.workDir
summary['Output directory path:']          = params.outdir
summary['Pipeline directory path:']          = workflow.projectDir
if( params.notrim ){
    summary['Trimmomatic Options: '] = 'Skipped trimming step'
} else {
    summary['Trimmomatic adapters:'] = params.trimmomatic_adapters_file
    summary['Trimmomatic adapter parameters:'] = params.trimmomatic_adapters_parameters
    summary["Trimmomatic read length (minimum):"] = params.trimmomatic_mininum_length
}
summary['Configuration Profile:'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "____________________________________________"

if (! params.withFastQC ) {
/*
 * Fastq File Processing
 * 
 * Fastqc
 */
process fastqc {
	label "small"
	tag "$prefix"
	publishDir "${params.outdir}/01-fastQC", mode: 'copy',
		saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

	input:
	set val(name), file(reads) from raw_reads_fastqc

	output:
	file '*_fastqc.{zip,html}' into fastqc_results
	file '.command.out' into fastqc_stdout

	script:

	prefix = name - ~/(_S[0-9]{2})?(_L00[1-9])?(.R1)?(_1)?(_R1)?(_trimmed)?(_val_1)?(_00*)?(\.fq)?(\.fastq)?(\.gz)?$/
	"""
	mkdir tmp
	fastqc -t ${task.cpus} -dir tmp $reads
	rm -rf tmp
	"""
}

/*
 * Trimm sequence reads
 * 
 * Trimmomatic
 */
process trimming {
	label "small"
	tag "$prefix"
	publishDir "${params.outdir}/02-preprocessing", mode: 'copy',
		saveAs: {filename ->
			if (filename.indexOf("_fastqc") > 0) "../03-preprocQC/$filename"
			else if (filename.indexOf(".log") > 0) "logs/$filename"
      else if (params.saveTrimmed && filename.indexOf(".fastq.gz")) "trimmed/$filename"
			else null
	}

	input:
	set val(name), file(reads) from raw_reads_trimming

	output:
	file '*_paired_*.fastq.gz' into trimmed_paired_reads,trimmed_paired_reads_bwa,virustrimmed_paired_reads_irus
	file '*_unpaired_*.fastq.gz' into trimmed_unpaired_reads
	file '*_fastqc.{zip,html}' into trimmomatic_fastqc_reports
	file '*.log' into trimmomatic_results

	script:
	prefix = name - ~/(_S[0-9]{2})?(_L00[1-9])?(.R1)?(_1)?(_R1)?(_trimmed)?(_val_1)?(_00*)?(\.fq)?(\.fastq)?(\.gz)?$/
	"""
	trimmomatic PE -threads ${task.cpus} -phred33 $reads $prefix"_paired_R1.fastq" $prefix"_unpaired_R1.fastq" $prefix"_paired_R2.fastq" $prefix"_unpaired_R2.fastq" ILLUMINACLIP:${params.trimmomatic_adapters_file}:${params.trimmomatic_adapters_parameters} SLIDINGWINDOW:${params.trimmomatic_window_length}:${params.trimmomatic_window_value} MINLEN:${params.trimmomatic_mininum_length} 2> ${name}.log

	gzip *.fastq
	mkdir tmp
	fastqc -t ${task.cpus} -q *_paired_*.fastq.gz
	rm -rf tmp

	"""
}

/*
 * STEPS 2.2 Mapping virus
 */
process map_virus {
	tag "$prefix"
	publishDir "${params.outdir}/map_virus", mode: 'copy',
		saveAs: {filename ->
			if (filename.indexOf(".bam") > 0) "mapping/$filename"
			else if (filename.indexOf(".bai") > 0) "mapping/$filename"
      else if (filename.indexOf(".txt") > 0) "stats/$filename"
      else if (filename.indexOf(".stats") > 0) "stats/$filename"
	}

	input:
	set file(readsR1),file(readsR2) from virustrimmed_paired_reads_irus
    file refvirus from virus_fasta_file
    file index from virus_index_files.collect()

	output:
	file '*_sorted.bam' into mapping_virus_sorted_bam,mapping_virus_sorted_bam_consensus
    file '*_consensus_masked.fasta' into masked_fasta

	script:
  prefix = readsR1.toString() - '_paired_R1.fastq.gz'
	"""
  bowtie2 -p ${task.cpus} --local -x $refvirus -1 $readsR1 -2 $readsR2 --very-sensitive-local -S $prefix".sam"
  samtools sort -o $prefix"_sorted.bam" -O bam -T $prefix $prefix".sam"
  samtools index $prefix"_sorted.bam"
  samtools flagstat $prefix"_sorted.bam" > $prefix"_flagstat.txt"
	"""
}

/*
 * STEPS 3.3 Consensus Genome
 */
process genome_consensus {
  tag "$prefix"
  publishDir "${params.outdir}/map_consensus", mode: 'copy',
		saveAs: {filename ->
			if (filename.indexOf("_consensus.fasta") > 0) "consensus/$filename"
			else if (filename.indexOf("_consensus_masked.fasta") > 0) "masked/$filename"
	}

  input:
  file refvirus from virus_fasta_file
  file sorted_bam from sorted_bam_consensus
  file sorted_bai from bai_consensus

  output:
  file '*_consensus.fasta' into consensus_fasta

  script:
  refname = refvirus.baseName - ~/(\.2)?(\.fasta)?$/
  """
  bcftools index "$refname".vcf.gz"
  cat $refvirus | bcftools consensus "$refname".vcf.gz" > "$refname"_consensus.fasta"
  bedtools genomecov -bga -ibam $sorted_bam -g $refvirus | awk '\$4 < 20' | bedtools merge > "$refname"_bed4mask.bed"
  bedtools maskfasta -fi "$refname"_consensus.fasta" -bed "$refname"_bed4mask.bed" -fo "$refname"_consensus_masked.fasta"
  sed -i 's/$refname/g' "$refname"_consensus_masked.fasta"
  """
}

