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
	/Users/kurtiscruz/Downloads/CURRENT/Virus_Genome_Mapping_Pipeline/Virus_Genome_Mapping_Pipeline/virus_ref_db/rhv_abc_sars2.fasta 
	
	PATH to indexed (indexed by bowtie2) Reference Fasta Database Files:
	/Users/kurtiscruz/Downloads/CURRENT/Virus_Genome_Mapping_Pipeline/Virus_Genome_Mapping_Pipeline/virus_ref_db

	PATH to fastq files from 031221 shotgun run (down-sized to 1m reads) for testing & debugging purposes:
	/Users/kurtiscruz/Downloads/CURRENT/test_fastq  ---> These are single-end reads

	CLI Commands

	Run Pipeline --helpMsg
	nextflow run /Users/Kurtisc/Downloads/CURRENT/Virus_Genome_Mapping_Pipeline/main.nf --helpMsg helpMsg

	Run Pipeline on test fastqc:
	
    Single end:

    nextflow run /Users/Kurtisc/Downloads/CURRENT/Virus_Genome_Mapping_Pipeline/main.nf --reads '/Users/kurtiscruz/Downloads/CURRENT/test_fastq_se/' --outdir '/Users/kurtiscruz/Downloads/CURRENT/test_output/' --singleEnd singleEnd

    Paired end:

    nextflow run /Users/kurtiscruz/Downloads/CURRENT/Virus_Genome_Mapping_Pipeline/main.nf --reads '/Users/kurtiscruz/Downloads/CURRENT/test_fastq_pe' --outdir '/Users/kurtiscruz/Downloads/CURRENT/test_output/'

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
      --withBlast                   Blasts the resulting consensus sequence and outputs results
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
params.virus_index = false
params.virus_fasta = false
params.withBlast =false
ADAPTERS = file("${baseDir}/All_adapters.fa")
MIN_LEN = 75
REFERENCE_FASTA = file("${baseDir}/virus_ref_db/rhv_ref_db01.fasta")
REFERENCE_FASTA_INDEX = file("${baseDir}/virus_ref_db/rhv_ref_db01.fasta.fai")
// Bowtie2 index name: rhv_ref_db01
BOWTIE2_DB_PREFIX = file("${baseDir}/virus_ref_db/rhv_ref_db01")
REF_BT2_INDEX1 = file("${baseDir}/virus_ref_db/rhv_ref_db01.1.bt2")
REF_BT2_INDEX2 = file("${baseDir}/virus_ref_db/rhv_ref_db01.2.bt2")
REF_BT2_INDEX3 = file("${baseDir}/virus_ref_db/rhv_ref_db01.3.bt2")
REF_BT2_INDEX4 = file("${baseDir}/virus_ref_db/rhv_ref_db01.4.bt2")
REF_BT2_INDEX5 = file("${baseDir}/virus_ref_db/rhv_ref_db01.rev.1.bt2")
REF_BT2_INDEX6 = file("${baseDir}/virus_ref_db/rhv_ref_db01.rev.2.bt2")
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
params.trimmomatic_adapters_file_PE = "/Users/Kurtisc/anaconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/TruSeq2-PE.fa"
params.trimmomatic_adapters_file_SE = "/Users/Kurtisc/anaconda3/pkgs/trimmomatic-0.39-1/share/trimmomatic-0.39-1/adapters/TruSeq2-SE.fa"
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
        tuple env(base),file("*.trimmed.fastq.gz") into Trim_out_ch2_SE

    publishDir "${params.outdir}trimmed_fastqs", mode: 'copy',pattern:'*.trimmed.fastq*'

    script:
    """
    #!/bin/bash
    base=`basename ${R1} ".fastq.gz"`
    echo \$base
	
	trimmomatic SE -threads ${task.cpus} ${R1} \$base.trimmed.fastq.gz \
	ILLUMINACLIP:${params.trimmomatic_adapters_file_SE}:${params.trimmomatic_adapters_parameters} SLIDINGWINDOW:${params.trimmomatic_window_length}:${params.trimmomatic_window_value} MINLEN:${params.trimmomatic_mininum_length} 2> ${R1}.log

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
	ILLUMINACLIP:${params.trimmomatic_adapters_file_PE}:${params.trimmomatic_adapters_parameters} SLIDINGWINDOW:${params.trimmomatic_window_length}:${params.trimmomatic_window_value} MINLEN:${params.trimmomatic_mininum_length} 2> ${R1}.log

    """
}
}
/*
 * Map sequence reads to local virus database
 */
process Mapping {
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
        tuple val(base), file("${base}.bowtie2.log") into Bowtie2_log_ch
        tuple val (base), file("*") into Dump_ch

    publishDir "${params.outdir}sam files", mode: 'copy', pattern:'*.sam*'
    publishDir "${params.outdir}bowtie2 logs", mode: 'copy', pattern:'*.bowtie2.log*'  

    script:

    """
    #!/bin/bash

    bowtie2 -p ${task.cpus} --local -x $BOWTIE2_DB_PREFIX -U ${base}.trimmed.fastq.gz --very-sensitive-local 2> ${base}.bowtie2.log -S ${base}.sam

    """
}

    // bowtie2 -p ${task.cpus} -x $BOWTIE2_DB_PREFIX ${base}.trimmed.fastq.gz 2> ${base}.bowtie2.log -S ${base}.sam

/*
 * STEP 5.2: Convert BAM to coordinate sorted BAM
 */
process Sort_Bam {
	errorStrategy 'retry'
    maxRetries 3

    input: 
    tuple val(base), file("${base}.sam") from Aligned_sam_ch

    output:
    tuple val(base), file("${base}.bam") into Aligned_bam_ch
    tuple val(base), file("${base}.sorted.bam") into Sorted_bam_ch
    tuple val(base), file("${base}.sorted.bam") into Sorted_bam_variant_ch
    tuple val(base), file("${base}.sorted.bam") into Sorted_Cons_Bam_ch
    tuple val(base), file("${base}.sorted.bam.bai") into Indexed_bam_ch
    tuple val(base), file("${base}.sorted.bam.bai") into Indexed_bam_variant_ch
    tuple val(base), file("${base}.sorted.bam.bai") into Indexed_Cons_Bam_ch
    tuple val(base), file("${base}.sorted.bam.flagstat") into Flagstats_ch
    tuple val(base), file("${base}.sorted.bam.idxstats") into Idxstats_ch
    tuple val(base), file("${base}.sorted.bam.stats") into Stats_ch

    publishDir "${params.outdir}bam files", mode: 'copy', pattern:'*.bam*'  
    publishDir "${params.outdir}sorted bam files", mode: 'copy', pattern:'*.sorted.bam*' 
    publishDir "${params.outdir}indexed bam files", mode: 'copy', pattern:'*.sorted.bam.bai*'  
    publishDir "${params.outdir}sorted bam flagstats", mode: 'copy', pattern:'*.sorted.bam.flagstat*' 
    publishDir "${params.outdir}sorted bam idxstats", mode: 'copy', pattern:'*.sorted.bam.idxstats*'  
    publishDir "${params.outdir}sorted bam stats", mode: 'copy', pattern:'*.sorted.bam.stats*'  

    script:
    """
    samtools view -S -b ${base}.sam > ${base}.bam
    samtools view ${base}.bam | head
    samtools sort ${base}.bam -o ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools flagstat ${base}.sorted.bam > ${base}.sorted.bam.flagstat
    samtools idxstats ${base}.sorted.bam > ${base}.sorted.bam.idxstats
    samtools stats ${base}.sorted.bam > ${base}.sorted.bam.stats

    """
}


/*
 * Variant Calling
 */
process Variant_Calling {
	errorStrategy 'retry'
    maxRetries 3

	input:
    tuple val(base), file("${base}.sorted.bam") from Sorted_bam_variant_ch
    tuple val(base), file("${base}.sorted.bam.bai") from Indexed_bam_variant_ch 
    file REFERENCE_FASTA
    file REFERENCE_FASTA_INDEX

	output:
    tuple val(base), file("${base}.pileup") into Variant_calling_pileup_ch
    tuple val(base), file("${base}_majority.vcf") into Mariant_calling_pileup_ch, Majority_allele_vcf_annotation, Majority_allele_vcf_consensus
    tuple val(base), file("${base}_lowfreq.vcf") into Lowfreq_variants_vcf, Lowfreq_variants_vcf_annotation, Lowfreq_variants_vcf_annotation2 
    
    publishDir "${params.outdir}variant calling pileup", mode: 'copy', pattern:'*.pileup*'  
    publishDir "${params.outdir}majority allele vcf", mode: 'copy', pattern:'*_majority.vcf*'  
    publishDir "${params.outdir}low freq variants vcf", mode: 'copy', pattern:'*_lowfreq.vcf*'  

	script:

	"""
    samtools mpileup -A -d 20000 -Q 0 -f $REFERENCE_FASTA ${base}.sorted.bam > ${base}.pileup
    varscan mpileup2cns ${base}.pileup --min-var-freq 0.02 --p-value 0.99 --variants --output-vcf 1 > ${base}_lowfreq.vcf
    varscan mpileup2cns ${base}.pileup --min-var-freq 0.8 --p-value 0.05 --variants --output-vcf 1 > ${base}_majority.vcf
	"""
}

/*
 * Variant Calling annotation
 */
process Variant_Calling_Annotation {
	errorStrategy 'retry'
    maxRetries 3

 	input:
    tuple val(base), file("${base}_lowfreq.vcf") from Lowfreq_variants_vcf_annotation 
    tuple val(base), file("${base}_majority.vcf") from Majority_allele_vcf_annotation

 	output:
    tuple val(base), file("${base}_majority.ann.vcf") into Majority_annotated_variants
    tuple val(base), file("${base}_majority_snpEff_genes.txt") into Majority_snpeff_summary
    tuple val(base), file("${base}_majority_snpEff_summary.html") into Lowfreq_annotated_variants 
    tuple val(base), file("${base}_lowfreq.ann.vcf") into Variant_calling_pileup_ch2
    tuple val(base), file("${base}_lowfreq_snpEff_genes.txt") into Lowfreq_snpeff_genes
    tuple val(base), file("${base}_lowfreq_snpEff_summary.html") into Lowfreq_snpeff_summary 
    tuple val(base), file("${base}_majority.ann.table.txt") into Snpsift_majority_table
    tuple val(base), file("${base}_lowfreq.ann.table.txt") into Snpsift_lowfreq_table

    publishDir "${params.outdir}majority variants", mode: 'copy', pattern:'*_majority.ann.vcf*'  
    publishDir "${params.outdir}majority snpEff genes", mode: 'copy', pattern:'*_majority_snpEff_genes.txt*'  
    publishDir "${params.outdir}majority snpEff summary", mode: 'copy', pattern:'*_majority_snpEff_summary.html*'  
    publishDir "${params.outdir}lowfreq ann", mode: 'copy', pattern:'*_lowfreq.ann.vcf*'  
    publishDir "${params.outdir}lowfreq snpEff genes", mode: 'copy', pattern:'*_lowfreq_snpEff_genes.txt*'  
    publishDir "${params.outdir}lowfreq snpEff summary", mode: 'copy', pattern:'*_lowfreq_snpEff_summary.html*'  
    publishDir "${params.outdir}majority ann table.", mode: 'copy', pattern:'*_majority.ann.table.txt*'  
    publishDir "${params.outdir}lowfreq ann table", mode: 'copy', pattern:'*_lowfreq.ann.table.txt*'  

 	script:

 	"""
    snpEff sars-cov-2 ${base}_majority.vcf > ${base}_majority.ann.vcf"
    mv snpEff_genes.txt ${base}_majority_snpEff_genes.txt"
    mv snpEff_summary.html ${base}_majority_snpEff_summary.html"
    snpEff sars-cov-2 ${base}_lowfreq.vcf > ${base}_lowfreq.ann.vcf"
    mv snpEff_genes.txt ${base}_lowfreq_snpEff_genes.txt"
    mv snpEff_summary.html ${base}_lowfreq_snpEff_summary.html"
    SnpSift extractFields -s "," -e "." ${base}_majority.ann.vcf" CHROM POS REF ALT "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].IMPACT" "ANN[*].EFFECT" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" > ${base}_majority.ann.table.txt
    SnpSift extractFields -s "," -e "." ${base}_lowfreq.ann.vcf" CHROM POS REF ALT "ANN[*].GENE" "ANN[*].GENEID" "ANN[*].IMPACT" "ANN[*].EFFECT" "ANN[*].FEATURE" "ANN[*].FEATUREID" "ANN[*].BIOTYPE" "ANN[*].RANK" "ANN[*].HGVS_C" "ANN[*].HGVS_P" "ANN[*].CDNA_POS" "ANN[*].CDNA_LEN" "ANN[*].CDS_POS" "ANN[*].CDS_LEN" "ANN[*].AA_POS" "ANN[*].AA_LEN" "ANN[*].DISTANCE" "EFF[*].EFFECT" "EFF[*].FUNCLASS" "EFF[*].CODON" "EFF[*].AA" "EFF[*].AA_LEN" > ${base}_lowfreq.ann.table.txt
 	
     """
}

process Consensus {
	errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(base), file("${base}_lowfreq.vcf") from Lowfreq_variants_vcf_annotation2 
    tuple val(base), file("${base}.sorted.bam") from Sorted_Cons_Bam_ch
    tuple val(base), file("${base}.sorted.bam.bai") from Indexed_Cons_Bam_ch
    file REFERENCE_FASTA
    file REFERENCE_FASTA_INDEX

    output:
    tuple val(base), file("${base}_consensus.fasta") into Consensus_fasta_ch
    tuple val(base), file("${base}_consensus_masked.fasta") into Consensus_masked_fasta_ch 

    publishDir "${params.outdir}consensus fasta", mode: 'copy', pattern:'*_consensus.fasta*'  
    publishDir "${params.outdir}consensus masked fasta", mode: 'copy', pattern:'*_consensus_masked.fasta*' 

    script:

    """
    bgzip -c $variants > ${base}_lowfreq.vcf
    bcftools index ${base}_lowfreq.vcf
    cat $REFERENCE_FASTA | bcftools consensus ${base}_lowfreq.vcf > ${base}_consensus.fasta"
    bedtools genomecov -bga -ibam ${base}.sorted.bam -g $REFERENCE_FASTA | awk '\$4 < 20' | bedtools merge > ${base}_bed4mask.bed"
    bedtools maskfasta -fi ${base}_consensus.fasta" -bed ${base}_bed4mask.bed" -fo ${base}_consensus_masked.fasta"
    sed -i 's/${base}/g' ${base}_consensus_masked.fasta"
    """
}

if (params.withFastQC) {
/*
 * STEP 3: FastQC on input reads after concatenating libraries from the same sample
 */
process FASTQC {
    tag "$sample"
    label 'process_medium'
    publishDir "${params.outdir}/preprocess/fastqc", mode: params.publish_dir_mode,
        saveAs: { filename ->
                      filename.endsWith(".zip") ? "zips/$filename" : filename
                }

    when:
    !params.skip_fastqc

    input:
    tuple val(sample), val(single_end), path(reads) from ch_cat_fastqc

    output:
    path "*.{zip,html}" into ch_fastqc_raw_reports_mqc

    script:
    """
    fastqc --quiet --threads $task.cpus *.fastq.gz
    """
}
}
