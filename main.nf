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
 		-Bowtie2 - align to reference Virus genome.
 		-Samtools - SAM and BAM file processing.
 - 3. : Consensus generation (*.fasta):
  		-Bcftools - consensus formatting.

Optional:
	-withFastQC
	-Add Blast consensus to NCBI? (or Vapid functionalities?)

Dependencies:
  
  trimmomatic
  bowtie2
  samtools
  bedtools
  Bcftools

	PATHS FOR EASTLAKE KC iMac - TESTING & DEBUGGING
	PATH to Virus Reference Fasta:
	/Users/Kurtisc/Downloads/CURRENT/Virus_Genome_Mapping_Pipeline/Virus_Genome_Mapping_Pipeline/virus_ref_db/rhv_abc_sars2.fasta 
	
	PATH to indexed (indexed by bowtie2) Reference Fasta Database Files:
	/Users/Kurtisc/Downloads/CURRENT/Virus_Genome_Mapping_Pipeline/Virus_Genome_Mapping_Pipeline/virus_ref_db

	PATH to fastq files from 031221 shotgun run (down-sized to 1m reads) for testing & debugging purposes:
	/Users/Kurtisc/Downloads/CURRENT/test_fastq  ---> These are single-end reads

	CLI Commands

	Run Pipeline --helpMsg
	nextflow run /Users/Kurtisc/Downloads/CURRENT/Virus_Genome_Mapping_Pipeline/Virus_Genome_Mapping_Pipeline/main.nf --helpMsg helpMsg

	Run Pipeline on test fastqc:
	
    Single end:

    nextflow run /Users/kurtiscruz/Downloads/CURRENT/Virus_Genome_Mapping_Pipeline/main.nf --reads '/Users/kurtiscruz/Downloads/CURRENT/test_fastq_se/' --virus_fasta /Users/kurtiscruz/Downloads/CURRENT/Virus_Genome_Mapping_Pipeline/virus_ref_db/rhv_abc_sars2.fasta --virus_index /Users/kurtiscruz/Downloads/CURRENT/Virus_Genome_Mapping_Pipeline/virus_ref_db/virus_DB1' --outdir '/Users/kurtiscruz/Downloads/CURRENT/test_output/' --singleEnd singleEnd

    Paired end:

    nextflow run /Users/kurtiscruz/Downloads/CURRENT/Virus_Genome_Mapping_Pipeline/main.nf --reads '/Users/kurtiscruz/Downloads/CURRENT/test_fastq_pe' --virus_fasta /Users/kurtiscruz/Downloads/CURRENT/Virus_Genome_Mapping_Pipeline/virus_ref_db/rhv_abc_sars2.fasta --virus_index /Users/kurtiscruz/Downloads/CURRENT/Virus_Genome_Mapping_Pipeline/virus_ref_db/virus_DB1' --outdir '/Users/kurtiscruz/Downloads/CURRENT/test_output/'

 ----------------------------------------------------------------------------------------
*/

// Pipeline version
version = '1.0'
def helpMsg() {
    log.info"""
	 __________________________________________________
     Virus Genome Mapping Pipeline :  Version ${version}
	__________________________________________________
    
	Pipeline Usage:

    To run the pipeline, enter the following in the command line:

        nextflow run Virus_Genome_Mapping_Pipeline/main.nf --reads PATH_TO_FASTQ --virus_fasta PATH_TO_VIR_FASTA --virus_index PATH_TO_VIR_INDEX --outdir PATH_TO_OUTPUT_DIR


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
	  --helpMsg						Displays help message in terminal

    """.stripIndent()
}
// Initialize parameters
params.helpMsg = false
ADAPTERS = file("${baseDir}/All_adapters.fa")
MIN_LEN = 75
REFERENCE_FASTA = file("${baseDir}/virus_ref_db/rhv_abc_sars2.fasta")
REFERENCE_FASTA_FAI = file("${baseDir}/virus_ref_db/rhv_abc_sars2.fasta.fai")
// Bowtie2 index name: virus_ref_db01

// ${baseDir}/virus_ref_db/virus_ref_db01.1.bt2 
// ${baseDir}/virus_ref_db/virus_ref_db01.2.bt2 
// ${baseDir}/virus_ref_db/virus_ref_db01.3.bt2 
// ${baseDir}/virus_ref_db/virus_ref_db01.4.bt2 
// ${baseDir}/virus_ref_db/virus_ref_db01.rev.1.bt2 
// ${baseDir}/virus_ref_db/virus_ref_db01.rev.2.bt2
BOWTIE2_DB_PREFIX = file("${baseDir}/virus_ref_db/virus_ref_db01")
REF_BT2_INDEX1 = file("${baseDir}/virus_ref_db/virus_ref_db01.1.bt2")
REF_BT2_INDEX2 = file("${baseDir}/virus_ref_db/virus_ref_db01.2.bt2")
REF_BT2_INDEX3 = file("${baseDir}/virus_ref_db/virus_ref_db01.3.bt2")
REF_BT2_INDEX4 = file("${baseDir}/virus_ref_db/virus_ref_db01.4.bt2")
REF_BT2_INDEX5 = file("${baseDir}/virus_ref_db/virus_ref_db01.rev.1.bt2")
REF_BT2_INDEX6 = file("${baseDir}/virus_ref_db/virus_ref_db01.rev.2.bt2")
// Show help msg
if (params.helpMsg){
    helpMsg()
    exit 0
}
params.withFastQC = false
// Check Nextflow version
nextflow_req_v = '20.10.0'
try {
    if( ! nextflow.version.matches(">= $nextflow_req_v") ){
        throw GroovyException("> ERROR: The version of Nextflow running on your machine is out dated.\n>Please update to Version $nextflow_req_v")
    }
} catch (all) {
	log.error"ERROR: This version of Nextflow is out of date.\nPlease update to the latest version of Nextflow."
}

// Check for virus genome reference indexes
params.virus_fasta = false
if( params.virus_fasta ){
    virus_fasta_file = file(params.virus_fasta)
    if( !virus_fasta_file.exists() ) exit 1, "> Virus fasta file not found: ${params.virus_fasta}.\n> Please specify a valid file path!"
}
// Check for fastq
params.reads = false
if (! params.reads ) exit 1, "> Error: Fastq files not found. Please specify a valid path with --reads"
// Single-end read option
params.singleEnd = false
// Trimming parameters
params.notrim = false
// Output files options
params.saveTrimmed = false
// Default trimming options
params.trimmomatic_adapters_file_PE = "/Users/kurtiscruz/opt/anaconda3/pkgs/trimmomatic-0.39-0/share/trimmomatic-0.39-0/adapters/TruSeq2-PE.fa"
params.trimmomatic_adapters_file_SE = "/Users/kurtiscruz/opt/anaconda3/pkgs/trimmomatic-0.39-0/share/trimmomatic-0.39-0/adapters/TruSeq3-SE.fa"
params.trimmomatic_adapters_parameters = "2:30:10:1"
params.trimmomatic_window_length = "4"
params.trimmomatic_window_value = "20"
params.trimmomatic_mininum_length = "75"
// log files header
log.info "____________________________________________"
log.info " Virus Genome Mapping Pipeline :  v${version}"
log.info "____________________________________________"
def summary = [:]
summary['Fastq Files:']               = params.reads
summary['Read type:']           	  = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Virus Reference:']           = params.virus_fasta
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current directory path:']        = "$PWD"
summary['Working directory path:']         = workflow.workDir
summary['Output directory path:']          = params.outdir
summary['Pipeline directory path:']          = workflow.projectDir
if( params.notrim ){
    summary['Trimmomatic Options: '] = 'Skipped trimming step'
} else {
    summary['Trimmomatic adapters:'] = params.trimmomatic_adapters_file_SE
	summary['Trimmomatic adapters:'] = params.trimmomatic_adapters_file_PE
    summary['Trimmomatic adapter parameters:'] = params.trimmomatic_adapters_parameters
    summary["Trimmomatic read length (minimum):"] = params.trimmomatic_mininum_length
}
summary['Configuration Profile:'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "____________________________________________"
// Create channel for input reads.
// Import reads depending on single end vs. paired end
if(params.singleEnd == false) {
    // Check for R1s and R2s in input directory
    input_read_ch = Channel
        .fromFilePairs("${params.reads}*_R{1,2}*.fastq.gz")
        .ifEmpty { error "> Cannot located paired-end reads in: ${params.reads}.\n> Please enter a valid file path." }
        .map { it -> [it[0], it[1][0], it[1][1]]}
} else {
    // Looks for gzipped files, assumes all separate samples
    input_read_ch = Channel
        .fromPath("${params.reads}*.gz")
        //.map { it -> [ file(it)]}
        .map { it -> file(it)}
}
if(params.virus_index) {
// Channel for virus genome reference indexes
Channel
    .fromPath(params.virus_index)
    .ifEmpty { exit 1, "> Error: Virus index not found: ${params.virus_index}.\n> Please specify a valid file path!"}
    .set { virus_index_files }
}
/*
 * Processing: Trim fastq sequence reads
 * 
 * Trimmomatic
 */
if (params.singleEnd) {
	process Trim_Reads_SE {
    errorStrategy 'retry'
    maxRetries 3

    input:
        file R1 from input_read_ch
        file ADAPTERS
        val MIN_LEN
    output: 
        tuple env(base),file("*.trimmed.fastq.gz"),file("*summary.csv") into Trim_out_ch_SE
        tuple env(base),file("*.trimmed.fastq.gz") into Trim_out_ch2_SE

    publishDir "${params.outdir}trimmed_fastqs", mode: 'copy',pattern:'*.trimmed.fastq*'

    script:
    """
    #!/bin/bash
    base=`basename ${R1} ".fastq.gz"`
    echo \$base

	printf "> Now Trimming: " $R1
	
	trimmomatic SE -threads ${task.cpus} ${R1} \$base.trimmed.fastq.gz \
	ILLUMINACLIP:${params.trimmomatic_adapters_file_SE}:${params.trimmomatic_adapters_parameters} SLIDINGWINDOW:${params.trimmomatic_window_length}:${params.trimmomatic_window_value} MINLEN:${params.trimmomatic_mininum_length} 2> ${R1}.log
	

	num_untrimmed=\$((\$(gunzip -c ${R1} | wc -l)/4))
    num_trimmed=\$((\$(gunzip -c \$base'.trimmed.fastq.gz' | wc -l)/4))
    
    percent_trimmed=\$((100-\$((100*num_trimmed/num_untrimmed))))
    echo "> \$base,\$num_untrimmed,\$num_trimmed,\$percent_trimmed" >> \$base'_summary.csv'
	echo "> Trimming completed succesfully. See log for details."

    """
} 
} else {
	process Trim_Reads_PE {
    errorStrategy 'retry'
    maxRetries 3

   input:
        tuple val(base), file(R1), file(R2) from input_read_ch
        file ADAPTERS
        val MIN_LEN
    output: 
        tuple val(base), file("${base}.trimmed.fastq.gz"),file("${base}_summary.csv") into Trim_out_ch
        tuple val(base), file(R1),file(R2),file("${base}.R1.paired.fastq.gz"), file("${base}.R2.paired.fastq.gz"),file("${base}.R1.unpaired.fastq.gz"), file("${base}.R2.unpaired.fastq.gz") into Trim_out_ch2
        tuple val(base), file("${base}.trimmed.fastq.gz") into Trim_out_ch3

    publishDir "${params.outdir}trimmed_fastqs", mode: 'copy',pattern:'*.trimmed.fastq*'
    
    script:
    """
    #!/bin/bash

    trimmomatic PE -threads ${task.cpus} ${R1} ${R2} ${base}.R1.paired.fastq.gz ${base}.R1.unpaired.fastq.gz ${base}.R2.paired.fastq.gz ${base}.R2.unpaired.fastq.gz \
	ILLUMINACLIP:${params.trimmomatic_adapters_file_SE}:${params.trimmomatic_adapters_parameters} SLIDINGWINDOW:${params.trimmomatic_window_length}:${params.trimmomatic_window_value} MINLEN:${params.trimmomatic_mininum_length} 2> ${R1}.log

    num_r1_untrimmed=\$(gunzip -c ${R1} | wc -l)
    num_r2_untrimmed=\$(gunzip -c ${R2} | wc -l)
    num_untrimmed=\$((\$((num_r1_untrimmed + num_r2_untrimmed))/4))
    num_r1_paired=\$(gunzip -c ${base}.R1.paired.fastq.gz | wc -l)
    num_r2_paired=\$(gunzip -c ${base}.R2.paired.fastq.gz | wc -l)
    num_paired=\$((\$((num_r1_paired + num_r2_paired))/4))
    num_r1_unpaired=\$(gunzip -c ${base}.R1.unpaired.fastq.gz | wc -l)
    num_r2_unpaired=\$(gunzip -c ${base}.R2.unpaired.fastq.gz | wc -l)
    num_unpaired=\$((\$((num_r1_unpaired + num_r2_unpaired))/4))
    num_trimmed=\$((num_paired + num_unpaired))
    percent_trimmed=\$((100-\$((100*num_trimmed/num_untrimmed))))
    cat *paired.fastq.gz > ${base}.trimmed.fastq.gz

    """
}
}

/*
 * Map sequence reads to local virus database
 */
process Genome_Mapping {
	errorStrategy 'retry'
    maxRetries 3

    input: 
        tuple val(base), file("${base}.trimmed.fastq.gz") from Trim_out_ch2_SE
        file REFERENCE_FASTA
		file REF_BT2_INDEX1
        file REF_BT2_INDEX2
        file REF_BT2_INDEX3
        file REF_BT2_INDEX4
        file REF_BT2_INDEX5
        file REF_BT2_INDEX6
    output:
        tuple val(base), file("${base}.sam") into Aligned_sam_ch
        tuple val(base), file("${base}.bam") into Aligned_bam_ch
        tuple val(base), file("${base}.sorted.bam") into Sorted_bam_ch
        tuple val (base), file("*") into Dump_ch

    publishDir "${params.outdir}sam files", mode: 'copy', pattern:'*.sam*'
    publishDir "${params.outdir}bam files", mode: 'copy', pattern:'*.bam*'  
    publishDir "${params.outdir}sorted bam files", mode: 'copy', pattern:'*.sorted.bam*'  
    
    script:

    """
    #!/bin/bash

    bowtie2 -p ${task.cpus} -x $BOWTIE2_DB_PREFIX ${base}.trimmed.fastq.gz --very-sensitive-local -S ${base}.sam

    samtools view -S -b ${base}.sam > ${base}.bam
    samtools view ${base}.bam | head
    samtools sort ${base}.bam -o ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    """
}
/*
 * Fastq File Processing
 * 
 * Fastqc
 */
// if (params.withFastQC) {
//   process FastQC {
// 	tag "$prefix"
// 	publishDir "${params.outdir}/fastQC", mode: 'copy',
// 		saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}
	
// 	input:
// 	set val(name), file(reads) from input_read_ch

// 	output:
// 	file '*_fastqc.{zip,html}' into fastqc_results

// 	script:

// 	prefix = name - ~/(_S[0-9]{2})?(_L00[1-9])?(.R1)?(_1)?(_R1)?(_trimmed)?(_val_1)?(_00*)?(\.fq)?(\.fastq)?(\.gz)?$/
// 	"""
// 	mkdir tmp
// 	fastqc -t ${task.cpus} -dir tmp $reads
// 	rm -rf tmp
// 	"""
// }
// }
/*
 * Next Steps (aside from fine tuning this thing and getting it running perfectly)
 *
 * CODE - Blast Consensus? Specify with --WithBlastConsensus
 *
 * CODE - Vapid-like genbank prep? Specify with --GenBankPrep
 *
 * To Do 
 * CODE: Add logic for --noTrim 
 * 
 */
