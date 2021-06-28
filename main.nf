#!/usr/bin/env nextflow

/*
========================================================================================
                  Human Respiratory Virus (HRV) Pipeline v1.2
========================================================================================
 Github Repo:
 Greninger Lab
 
 Author:
 Paul RK Cruz <kurtisc@uw.edu>
 Michelle J Lin <Mjlin@uw.edu>
 Alex L Greninger <agrening@uw.edu>
 UW Medicine | Virology
 Department of Laboratory Medicine and Pathology
 University of Washington
 Created: April, 2021
 LICENSE: GNU
----------------------------------------------------------------------------------------
This pipeline was designed to run either single-end or paired end Next-Generation Sequencing reads to identify Human Respiratory Virus Complete Genomes for analysis and Genbank submission. 
The pipeline outputs mapping statistics, sam, bam, and consensus files, as well as a final report summary at the end of the workflow for analysis.

PIPELINE OVERVIEW:
 - 1. : Trim Reads
 		-Trimmomatic - sequence read trimming of adaptors and low quality reads.
 - 2. : Genome Mapping
 		-BBMap - align to MultiFasta Reference Virus Genome.
 		-Samtools - SAM and BAM file processing.
 - 3. : Reference Fasta Generation
 		-Generate a fasta reference from the genome mapping results.  
 - 4. : Sort Bam
  		-Convert Sam to Bam
        -Sort Bam file by coordinates
        -Generate Statistics about Bam file  
 - 5. : Variant Calling
        -Calculate the read coverage of positions in the genome
        -Detect the single nucleotide polymorphisms (SNPs)
        -Filter and report the SNP variants in variant calling format (VCF)
        CLI Command to view results:   less -S ${base}_final_variants.vcf
 - 6. : Consensus
        -Consensus generation using variants VCF, mapped reference fasta, and
        sorted bam. 
 - 7. : Final Consensus
        -Creates the Final Consensus by editing the fasta header.       
 - 6. : FastQC
 		-Sequence read quality control analysis.
    Dependencies:
    
    trimmomatic
    samtools
    bbtools  
    bcftools
    seqkit
    bgzip
    bedtools
    fastqc
    
    PIPELINE SETUP

    Setup fasta Reference:
    1. Reference_Fasta (Multifasta containing concatenated full length RhV genome sequences (6-10K bp) formatted with accession numbers only)
        Current file: rhv_ref_db01_accession_only.fasta - 500~ Human Rhinovirus Complete Genome Sequences courtesy of NCBI Genbank, 2021.
            source: https://www.ncbi.nlm.nih.gov/nucleotide/
    This pipeline can also run using one fasta reference.

    Setup File Paths:
    1. BBMAP_PATH
        Path to your installation of BBTools --> bbmap.sh
    2. trimmomatic_adapters_file_SE
        Path to your Trimmomatic single-end file
    3. trimmomatic_adapters_file_PE
            Path to your Trimmomatic paired-end file
    Setup Trimmomatic Parameters:
    1. params.trimmomatic_adapters_parameters = "2:30:10:1"
    2. params.trimmomatic_window_length = "4"
    3. params.trimmomatic_window_value = "20"
    4. params.trimmomatic_mininum_length = "75"
    
    EXAMPLE USAGE:
        Run Pipeline Help Message:
        nextflow run /Users/Kurtisc/Downloads/CURRENT/Virus_Genome_Mapping_Pipeline/main.nf --helpMsg helpMsg
/Users/uwvirongs/Documents/KC/input
        Run Pipeline on Single-end sequence reads ((SAMPLE_NAME)_S1_L001_R1_001.fastq, ((SAMPLE_NAME)_S1_L002_R1_001.fastq))
        nextflow run /Users/uwvirongs/Documents/KC/HRV_Genome_Mapping_Pipeline/main.nf --reads '/Users/uwvirongs/Documents/KC/input/' --outdir '/Users/uwvirongs/Documents/KC/input/' --singleEnd
        Run Pipeline on Paired-end sequence reads ((SAMPLE_NAME)_S1_L001_R1_001.fastq, ((SAMPLE_NAME)_S1_L001_R2_001.fastq))
        nextflow run /Users/Kurtisc/Downloads/CURRENT/Virus_Genome_Mapping_Pipeline/Virus_Genome_Mapping_Pipeline/main.nf --reads '/Users/Kurtisc/Downloads/CURRENT/test_fastq_pe/' --outdir '/Users/Kurtisc/Downloads/CURRENT/test_output/'
 ----------------------------------------------------------------------------------------
*/
// Pipeline version
params.helpMsg = false
version = '1.2'
def helpMsg() {
    log.info"""
	 _______________________________________________________________________________
     Human Respiratory Virus (HRV) Pipeline :  Version ${version}
	________________________________________________________________________________
    
	Usage:
    To run the pipeline, enter the following in the command line:

        nextflow run FILE_PATH/HRV_Genome_Mapping_Pipeline/main.nf --reads PATH_TO_FASTQ --outdir PATH_TO_OUTPUT_DIR --Reference_Fasta PATH_TO_REF_FASTA 
    
    Valid command-line arguments:
    REQUIRED:
      --reads                       Path to input fastq.gz folder (files must be *.gz)
      --outdir                      Path to output directory to store pipeline output files
    OPTIONAL:
      --withSampleSheet             Add FASTQ sample information to final summary report *.csv
      --singleEnd                   Specifies that the input fastq files are single-end reads
      --withFastQC                  Runs a quality control check on fastq files
      --helpMsg                     Display help message in terminal      
    """.stripIndent()
}
// Show help msg
if (params.helpMsg){
    helpMsg()
    exit 0
}
// Initialize parameters
params.virus_index = false
params.virus_fasta = false
params.withFastQC = false
params.skipTrim = false
params.reads = false
params.singleEnd = false
params.ADAPTERS = false
params.withSampleSheet = false
// Reference Files

// Rhinovirus Reference Multi-Fasta
Reference_Fasta_hrv = file("${baseDir}/hrv-ref/hrv_ref_rhinovirus.fa")
// Respiratoery Panel Reference Multi-Fasta
Reference_Fasta_hcov = file("${baseDir}/hrv-ref/hrv_ref_db01_accession_only-resp-panel.fa")
// HCoV Reference Multi-Fasta
Reference_Fasta_hcov = file("${baseDir}/hrv-ref/hrv_ref_hcov.fa")
// Influenza B Reference Multi-Fasta
Reference_Fasta_hcov = file("${baseDir}/hrv-ref/hrv_ref_Influenza_B.fa")
// Script Files
if(params.withSampleSheet != false) {
    SAMPLE_SHEET = file(params.withSampleSheet)
}
TRIM_ENDS=file("${baseDir}/scripts/trim_ends.py")
VCFUTILS=file("${baseDir}/scripts/vcfutils.pl")
SPLITCHR=file("${baseDir}/scripts/splitchr.txt")
FIX_COVERAGE = file("${baseDir}/scripts/fix_coverage.py")
ADAPTERS_SE = file("${baseDir}/adapters/TruSeq2-SE.fa")
ADAPTERS_PE = file("${baseDir}/adapters/TruSeq2-PE.fa")
// BBMap Path
BBMAP_PATH="/Users/greningerlab/Documents/bbmap/"
// BBMAP_PATH="/Volumes/HD2/bbmap/"
params.MINLEN = "35"
MINLEN = "35"
// Check Nextflow version
nextflow_req_v = '20.10.0'
try {
    if( ! nextflow.version.matches(">= $nextflow_req_v") ){
        throw GroovyException("> ERROR: The version of Nextflow running on your machine is out dated.\n>Please update to Version $nextflow_req_v")
    }
} catch (all) {
	log.error"ERROR: This version of Nextflow is out of date.\nPlease update to the latest version of Nextflow."
}

if (! params.reads ) exit 1, "> Error: Fastq files not found. Please specify a valid path with --reads"
// log files header
log.info "_______________________________________________________________________________"
log.info " Human Rhinovirus Genome Mapping Pipeline :  v${version}"
log.info "_______________________________________________________________________________"
def summary = [:]
summary['Fastq Files:']               = params.reads
summary['Sample Sheet:']               = params.withSampleSheet
summary['Read type:']           	  = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Virus Reference:']           = Reference_Fasta
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current directory path:']        = "$PWD"
summary['Working directory path:']         = workflow.workDir
summary['Output directory path:']          = params.outdir
summary['Pipeline directory path:']          = workflow.projectDir
if (params.singleEnd) {
// summary['Trimmomatic adapters:'] = params.trimmomatic_adapters_file_SE
} else {
// summary['Trimmomatic adapters:'] = params.trimmomatic_adapters_file_PE
}
summary["Trimmomatic read length (minimum):"] = params.MINLEN
summary['Configuration Profile:'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "_______________________________________________________________________________"
// Create channel for input reads.
// Import reads depending on single-end or paired-end
if(params.singleEnd == false) {
    // Check for R1s and R2s in input directory
    input_read_ch = Channel
        .fromFilePairs("${params.reads}*_R{1,2}*.fastq.gz")
        // .ifEmpty { error "> Cannot locate paired-end reads in: ${params.reads}.\n> Please enter a valid file path." }
        .map { it -> [it[0], it[1][0], it[1][1]]}
} else {
    // input: *.gz
    input_read_ch = Channel
        .fromPath("${params.reads}*.gz")
        .ifEmpty { error "> Cannot locate single-end reads in: ${params.reads}.\n> Please enter a valid file path." }
        .map { it -> file(it)}
}
/*
 * Trim Reads
 * 
 * Processing: Trim adaptors and repetitive bases from sequence reads and remove low quality sequence reads.
 */
// if(!params.skipTrimming) {
if (params.singleEnd) {
	process Trim_Reads {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"
    // errorStrategy 'retry'
    // maxRetries 3


    input:
        file R1 from input_read_ch
        file ADAPTERS_SE
        val MINLEN

    output: 
        tuple env(base),file("*.trimmed.fastq.gz"), file("${R1}_num_trimmed.txt"),file("*summary.csv") into Trim_out_SE, Trim_out_SE_FQC

    publishDir "${params.outdir}trimmed_fastqs", mode: 'copy',pattern:'*.trimmed.fastq*'


    script:
    """
    #!/bin/bash
    base=`basename ${R1} ".fastq.gz"`
    echo \$base
    /usr/local/miniconda/bin/trimmomatic SE -threads ${task.cpus} ${R1} \$base.trimmed.fastq.gz \
    ILLUMINACLIP:${ADAPTERS_SE}:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:${MINLEN}
    num_untrimmed=\$((\$(gunzip -c ${R1} | wc -l)/4))
    num_trimmed=\$((\$(gunzip -c \$base'.trimmed.fastq.gz' | wc -l)/4))
    printf "\$num_trimmed" >> ${R1}_num_trimmed.txt
    percent_trimmed=\$((100-\$((100*num_trimmed/num_untrimmed))))
    echo Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Reference_Genome,Reference_Length,Mapped_Reads,Percent_Ref_Coverage,Min_Coverage,Mean_Coverage,Max_Coverage,Bam_Size,Consensus_Length,Percent_N,%_Reads_On_Target, PCR_CT,Method, NCBI_Name> \$base'_summary.csv'
    printf "\$base,\$num_untrimmed,\$num_trimmed,\$percent_trimmed" >> \$base'_summary.csv'
    ls -latr
    """
} 
} else {
	process Trim_Reads_PE {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"
    errorStrategy 'retry'
    maxRetries 3

   input:
        tuple val(base), file(R1), file(R2) from input_read_ch
        file ADAPTERS_PE
        val MINLEN
    output: 
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${R1}_num_trimmed.txt"),file("*summary.csv") into Trim_out_PE, Trim_out_PE_FQC
        tuple val(base), file(R1),file(R2),file("${base}.R1.paired.fastq.gz"), file("${base}.R2.paired.fastq.gz"),file("${base}.R1.unpaired.fastq.gz"), file("${base}.R2.unpaired.fastq.gz") into Trim_out_ch2
        tuple val(base), file("${base}.trimmed.fastq.gz") into Trim_out_ch3


    publishDir "${params.outdir}fastq_trimmed", mode: 'copy',pattern:'*.trimmed.fastq*'
    
    script:
    """
    #!/bin/bash


    /usr/local/miniconda/bin/trimmomatic PE -threads ${task.cpus} ${R1} ${R2} ${base}.R1.paired.fastq.gz ${base}.R1.unpaired.fastq.gz ${base}.R2.paired.fastq.gz ${base}.R2.unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTERS_PE}:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:${MINLEN}

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
    printf "\$num_trimmed" >> ${R1}_num_trimmed.txt
    
    echo Sample_Name,Raw_Reads,Trimmed_Paired_Reads,Trimmed_Unpaired_Reads,Total_Trimmed_Reads, Percent_Trimmed,Reference_Genome,Reference_Length,Mapped_Reads,Percent_Ref_Coverage,Min_Coverage,Mean_Coverage,Max_Coverage,Bam_Size,Consensus_Length,Percent_N,%_Reads_On_Target, PCR_CT,Method, NCBI_Name> \$base'_summary.csv'
    printf "${base},\$num_untrimmed,\$num_paired,\$num_unpaired,\$num_trimmed,\$percent_trimmed" >> \$base'_summary.csv'
    ls -latr

    cat *paired.fastq.gz > ${base}.trimmed.fastq.gz

    """
}
}
/*
 * Map sequence reads to HRV Genomes using BBMap.
 */
 if (params.singleEnd) {
process Mapping {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"
    errorStrategy 'retry'
    maxRetries 3

    input: 
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_summary.csv") from Trim_out_SE
        file Reference_Fasta_hrv

    output:
        tuple val(base), file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"), file("${base}_summary2.csv"), file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),env(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),env(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt") into Everything_ch
        tuple val(base), file("${base}_map1_histogram.txt"),file("${base}_map2_histogram.txt") into BBmap_map1_hist_ch
        tuple val (base), file("*") into Dump_map1_ch

    publishDir "${params.outdir}sam_map1", mode: 'copy', pattern:'*_map1.sam*'
    publishDir "${params.outdir}sam_map2", mode: 'copy', pattern:'*_map2.sam*'
    publishDir "${params.outdir}txt_bbmap_map1_stats", mode: 'copy', pattern:'*_map1_bbmap_out.txt*'  
    publishDir "${params.outdir}txt_bbmap_map1_hist", mode: 'copy', pattern:'*_map2_histogram.txt*' 
    publishDir "${params.outdir}txt_bbmap_map2_stats", mode: 'copy', pattern:'*_map2_bbmap_out.txt*'  
    publishDir "${params.outdir}txt_bbmap_map2_hist", mode: 'copy', pattern:'*_map2_histogram.txt*'
    publishDir "${params.outdir}txt_indxstats_mapped_refs", mode: 'copy', pattern:'*_idxstats.txt*'   
    publishDir "${params.outdir}ref_id", mode: 'copy', pattern:'*_most_mapped_ref.txt*'  
    publishDir "${params.outdir}ref_fasta", mode: 'copy', pattern:'*_mapped_ref_genome.fa*'
    publishDir "${params.outdir}ref_size", mode: 'copy', pattern:'*_most_mapped_ref_size.txt*'
    publishDir "${params.outdir}ref_fai_index", mode: 'copy', pattern:'*_mapped_ref_genome.fa.fai*'
    
    script:

    """
    #!/bin/bash
    
    bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${Reference_Fasta_hrv} threads=${task.cpus} covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map1_stats.txt 2>&1

    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${Reference_Fasta_hrv} \$id > ${base}_mapped_ref_genome.fa




    bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map2.sam ref=${base}_mapped_ref_genome.fa threads=${task.cpus} covstats=${base}_map2_bbmap_out.txt covhist=${base}_map2_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map2_stats.txt 2>&1
    head -n 1 ${base}_mapped_ref_genome.fa > ${base}_mapped_ref_genome_edited.fa
    grep -v ">" ${base}_mapped_ref_genome.fa | sed 's/U/T/g' >> ${base}_mapped_ref_genome_edited.fa
    mv ${base}_mapped_ref_genome_edited.fa ${base}_mapped_ref_genome.fa
    samtools faidx ${base}_mapped_ref_genome.fa
    awk 'NR == 2 || \$5 > max {number = \$3; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref_size_out.txt
    id_ref_size=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref_size_out.txt)
    echo \$id_ref_size >> ${base}_most_mapped_ref_size.txt
    reads_mapped=\$(cat ${base}_map2_stats.txt | grep "mapped:" | cut -d\$'\\t' -f3)
    printf "\$reads_mapped" >> ${base}_num_mapped.txt
    cp ${base}_summary.csv ${base}_summary2.csv
    printf ",\$id" >> ${base}_summary2.csv
    printf ",\$id_ref_size" >> ${base}_summary2.csv
    printf ",\$reads_mapped" >> ${base}_summary2.csv
    printf ",\$ref_coverage" >> ${base}_summary2.csv
    """
}
} else {
process Mapping_PE {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest" 
    errorStrategy 'retry'
    maxRetries 3

    input: 
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_summary.csv") from Trim_out_PE
        file Reference_Fasta_hrv

    output:
        tuple val(base), file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"), file("${base}_summary2.csv"), file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),env(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),env(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt") into Everything_PE_ch
        tuple val(base), file("${base}_map1_histogram.txt"),file("${base}_map2_histogram.txt") into BBmap_map1_hist_PE_ch
        tuple val (base), file("*") into Dump_map1_ch

    publishDir "${params.outdir}sam_map1", mode: 'copy', pattern:'*_map1.sam*'
    publishDir "${params.outdir}sam_map2", mode: 'copy', pattern:'*_map2.sam*'
    publishDir "${params.outdir}txt_bbmap_map1_stats", mode: 'copy', pattern:'*_map1_bbmap_out.txt*'  
    publishDir "${params.outdir}txt_bbmap_map1_hist", mode: 'copy', pattern:'*_map2_histogram.txt*' 
    publishDir "${params.outdir}txt_bbmap_map2_stats", mode: 'copy', pattern:'*_map2_bbmap_out.txt*'  
    publishDir "${params.outdir}txt_bbmap_map2_hist", mode: 'copy', pattern:'*_map2_histogram.txt*'
    publishDir "${params.outdir}txt_indxstats_mapped_refs-map1", mode: 'copy', pattern:'*_idxstats.txt*'   
    publishDir "${params.outdir}ref_id", mode: 'copy', pattern:'*_most_mapped_ref.txt*'  
    publishDir "${params.outdir}ref_fasta", mode: 'copy', pattern:'*_mapped_ref_genome.fa*'
    publishDir "${params.outdir}ref_size", mode: 'copy', pattern:'*_most_mapped_ref_size.txt*'
    publishDir "${params.outdir}ref_fai_index", mode: 'copy', pattern:'*_mapped_ref_genome.fa.fai*'
    
    script:

    """
    #!/bin/bash
    bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${Reference_Fasta_hrv} threads=${task.cpus} covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt maxindel=20 strictmaxindel local=true -Xmx6g > ${base}_map1_stats.txt 2>&1

    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${Reference_Fasta_hrv} \$id > ${base}_mapped_ref_genome.fa

    bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map2.sam ref=${base}_mapped_ref_genome.fa threads=${task.cpus} covstats=${base}_map2_bbmap_out.txt covhist=${base}_map2_histogram.txt local=true -Xmx6g > ${base}_map2_stats.txt 2>&1
    head -n 1 ${base}_mapped_ref_genome.fa > ${base}_mapped_ref_genome_edited.fa
    grep -v ">" ${base}_mapped_ref_genome.fa | sed 's/U/T/g' >> ${base}_mapped_ref_genome_edited.fa
    mv ${base}_mapped_ref_genome_edited.fa ${base}_mapped_ref_genome.fa
    samtools faidx ${base}_mapped_ref_genome.fa
    awk 'NR == 2 || \$5 > max {number = \$3; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref_size_out.txt
    id_ref_size=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref_size_out.txt)
    echo \$id_ref_size >> ${base}_most_mapped_ref_size.txt
    reads_mapped=\$(cat ${base}_map2_stats.txt | grep "mapped:" | cut -d\$'\\t' -f3)
    printf "\$reads_mapped" >> ${base}_num_mapped.txt
    cp ${base}_summary.csv ${base}_summary2.csv
    printf ",\$id" >> ${base}_summary2.csv
    printf ",\$id_ref_size" >> ${base}_summary2.csv
    printf ",\$reads_mapped" >> ${base}_summary2.csv
    printf ",\$ref_coverage" >> ${base}_summary2.csv
    """
}
}

/*
 * Convert BAM to coordinate sorted BAM
 */
 // Step 1. Convert Sam to Bam
 // Step 2. Sort Bam file by coordinates
 // Step 3. Generate Statistics about Bam file
if (params.singleEnd) {
process Sort_Bam {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"    
	errorStrategy 'retry'
    maxRetries 3

    input: 
    tuple val(base), file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"), file("${base}_summary2.csv"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt") from Everything_ch
    output:
    tuple val(base), file("${base}.bam") into Aligned_bam_ch, Bam_ch
    tuple val(base), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),env(bamsize),file("${base}.sorted.bam.bai"),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"), file("${base}_summary.csv"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt") into Consensus_ch

    publishDir "${params.outdir}bam-map2", mode: 'copy', pattern:'*.bam*'
    publishDir "${params.outdir}bam_sorted-map2", mode: 'copy', pattern:'*.sorted.bam*'  
    publishDir "${params.outdir}txt_bam_flagstats-map2", mode: 'copy', pattern:'*_flagstats.txt*' 

    script:
    """
    #!/bin/bash
    samtools view -S -b ${base}_map2.sam > ${base}.bam
    samtools sort -@ ${task.cpus} ${base}.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools flagstat ${base}.sorted.bam > ${base}_flagstats.txt
    bedtools genomecov -d -ibam ${base}.sorted.bam > ${base}_coverage.txt
    
    awk 'NR == 3 || \$3 > max {number = \$3; max = \$1} END {if (NR) print number, max}' < ${base}_coverage.txt > ${base}_min_coverage.txt
    awk 'NR == 2 || \$3 > min {number = \$1; min = \$3} END {if (NR) print number, min}' < ${base}_coverage.txt > ${base}_max_coverage.txt
    
    meancoverage=\$(cat ${base}_coverage.txt | awk '{sum+=\$3} END { print sum/NR}')
    mincoverage=\$(awk 'FNR==1{print val,\$1}' ${base}_min_coverage.txt)
    maxcoverage=\$(awk 'FNR==1{print val,\$2}' ${base}_max_coverage.txt)
    bamsize=\$((\$(wc -c ${base}.sorted.bam | awk '{print \$1'})+0))

    cp ${base}_summary2.csv ${base}_summary.csv
    printf ",\$mincoverage" >> ${base}_summary.csv
    printf ",\$meancoverage" >> ${base}_summary.csv
    printf ",\$maxcoverage" >> ${base}_summary.csv
    printf ",\$bamsize" >> ${base}_summary.csv

    """
}
} else {
process Sort_Bam_PE {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"  
	errorStrategy 'retry'
    maxRetries 3

    input: 
    tuple val(base), file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"), file("${base}_summary2.csv"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt") from Everything_PE_ch
    output:
    tuple val(base), file("${base}.bam") into Aligned_bam_PE_ch, Bam_PE_ch
    tuple val(base), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),env(bamsize),file("${base}.sorted.bam.bai"),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"), file("${base}_summary.csv"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt") into Consensus_PE_ch

    publishDir "${params.outdir}bam-map2", mode: 'copy', pattern:'*.bam*'
    publishDir "${params.outdir}bam_sorted-map2", mode: 'copy', pattern:'*.sorted.bam*'  
    publishDir "${params.outdir}txt_bam_flagstats-map2", mode: 'copy', pattern:'*_flagstats.txt*'  

    script:
    """
    #!/bin/bash
    samtools view -S -b ${base}_map2.sam > ${base}.bam
    samtools sort -@ ${task.cpus} ${base}.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools flagstat ${base}.sorted.bam > ${base}_flagstats.txt
    bedtools genomecov -d -ibam ${base}.sorted.bam > ${base}_coverage.txt
    
    awk 'NR == 3 || \$3 > max {number = \$3; max = \$1} END {if (NR) print number, max}' < ${base}_coverage.txt > ${base}_min_coverage.txt
    awk 'NR == 2 || \$3 > min {number = \$1; min = \$3} END {if (NR) print number, min}' < ${base}_coverage.txt > ${base}_max_coverage.txt
    
    meancoverage=\$(cat ${base}_coverage.txt | awk '{sum+=\$3} END { print sum/NR}')
    mincoverage=\$(awk 'FNR==1{print val,\$1}' ${base}_min_coverage.txt)
    maxcoverage=\$(awk 'FNR==1{print val,\$2}' ${base}_max_coverage.txt)
    bamsize=\$((\$(wc -c ${base}.sorted.bam | awk '{print \$1'})+0))

    cp ${base}_summary2.csv ${base}_summary.csv
    printf ",\$mincoverage" >> ${base}_summary.csv
    printf ",\$meancoverage" >> ${base}_summary.csv
    printf ",\$maxcoverage" >> ${base}_summary.csv
    printf ",\$bamsize" >> ${base}_summary.csv

    """
}
}
/*
 * Call variants & Generate Consensus
  */
if (params.singleEnd) {
 process Generate_Consensus {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(base), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),val(bamsize),file("${base}.sorted.bam.bai"),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"), file("${base}_summary.csv"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt") from Consensus_ch
    
    file VCFUTILS
    file SPLITCHR
    file TRIM_ENDS

    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}.consensus.fa"), file("${base}_summary.csv"), val(bamsize), val(id),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt") into Consensus_Fasta_ch
    tuple val(base), file("${base}_pre_bcftools.vcf"), file("${base}_bcftools.vcf") into Consensus_Vcf_ch

    publishDir "${params.outdir}consensus-unmasked", mode: 'copy', pattern:'*.consensus.fa*' 
    publishDir "${params.outdir}vcf-map2", mode: 'copy', pattern:'*_bcftools.vcf*' 
    publishDir "${params.outdir}vcf_pre-map2", mode: 'copy', pattern:'*_pre_bcftools.vcf*' 

    shell:
    '''
    #!/bin/bash
    ls -latr
    R1=!{base}
    echo "bamsize: !{bamsize}"
    #if [ -s !{} ]
    # More reliable way of checking bam size, because of aliases
    if (( !{bamsize} > 92 ))
    then
        # Parallelize pileup based on number of cores
        splitnum=$(($((!{id_ref_size}/!{task.cpus}))+1))
        perl !{VCFUTILS} splitchr -l $splitnum !{base}_mapped_ref_genome.fa.fai | \\
        #cat !{SPLITCHR} | \\
            xargs -I {} -n 1 -P !{task.cpus} sh -c \\
                "/usr/local/miniconda/bin/bcftools mpileup \\
                    -f !{base}_mapped_ref_genome.fa -r {} \\
                    --count-orphans \\
                    --no-BAQ \\
                    --max-depth 50000 \\
                    --max-idepth 500000 \\
                    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \\
                !{base}.sorted.bam | /usr/local/miniconda/bin/bcftools call -A -m -M -Oz - > tmp.{}.vcf.gz"
        
        # Concatenate parallelized vcfs back together
        gunzip tmp*vcf.gz
        mv tmp.*:1-* \${R1}_catted.vcf
        for file in tmp*.vcf; do grep -v "#" $file >> \${R1}_catted.vcf; done
        cat \${R1}_catted.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' | /usr/local/miniconda/bin/bcftools norm -m -any > \${R1}_pre_bcftools.vcf
        
        # Make sure variants are majority variants for consensus calling
        /usr/local/miniconda/bin/bcftools filter -i '(DP4[0]+DP4[1]) < (DP4[2]+DP4[3]) && ((DP4[2]+DP4[3]) > 0) | (IMF > 0.5)' --threads !{task.cpus} \${R1}_pre_bcftools.vcf -o \${R1}.vcf
       # /usr/local/miniconda/bin/bcftools filter -e 'IMF < 0.5' \${R1}_pre2.vcf -o \${R1}.vcf
        # Index and generate consensus from vcf with majority variants
        /usr/local/miniconda/bin/bgzip \${R1}.vcf
        /usr/local/miniconda/bin/tabix \${R1}.vcf.gz 
        cat !{base}_mapped_ref_genome.fa | /usr/local/miniconda/bin/bcftools consensus \${R1}.vcf.gz > \${R1}.consensus.fa
        
        # Create coverage file from bam for whole genome, then pipe anything that has less than 6 coverage to bed file,
        # to be masked later
        /usr/local/miniconda/bin/bedtools genomecov \\
            -bga \\
            -ibam !{base}.sorted.bam \\
            -g !{base}_mapped_ref_genome.fa  \\
            | awk '\$4 < 6' | /usr/local/miniconda/bin/bedtools merge > \${R1}.mask.bed
        # Get rid of anything outside of the genome we care about, to prevent some sgrnas from screwing with masking
        awk '{ if(\$3 > 200 && \$2 < 29742) {print}}' \${R1}.mask.bed > a.tmp && mv a.tmp \${R1}.mask.bed
        # Mask refseq fasta for low coverage areas based on bed file
        /usr/local/miniconda/bin/bedtools maskfasta \\
            -fi !{base}_mapped_ref_genome.fa  \\
            -bed \${R1}.mask.bed \\
            -fo ref.mask.fasta
        
        # Align to refseq and unwrap fasta
        cat ref.mask.fasta \${R1}.consensus.fa > align_input.fasta
        /usr/local/miniconda/bin/mafft --auto --thread !{task.cpus} align_input.fasta > repositioned.fasta
        awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' repositioned.fasta > repositioned_unwrap.fasta
        # Trim ends and aligns masking of refseq to our consensus
        python3 !{TRIM_ENDS} \${R1}
        gunzip \${R1}.vcf.gz
        mv \${R1}.vcf \${R1}_bcftools.vcf
        sed -i 's/>.*/>!{base}.consensus/' \${R1}.consensus.fa
    else
       echo "Empty bam detected. Generating empty consensus fasta file..."
       touch \${R1}_bcftools.vcf
       touch \${R1}_pre_bcftools.vcf
    fi
    '''
}
} else {
 process Generate_Consensus_PE {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"   
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(base), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),val(bamsize),file("${base}.sorted.bam.bai"),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"), file("${base}_summary.csv"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt") from Consensus_PE_ch
    
    file VCFUTILS
    file SPLITCHR
    file TRIM_ENDS

    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}.consensus.fa"), file("${base}_summary.csv"), val(bamsize), val(id),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt") into Consensus_Fasta_PE_ch
    tuple val(base), file("${base}_pre_bcftools.vcf"), file("${base}_bcftools.vcf") into Consensus_Vcf_PE_ch

    publishDir "${params.outdir}consensus-unmasked", mode: 'copy', pattern:'*.consensus.fa*' 
    publishDir "${params.outdir}vcf-map2", mode: 'copy', pattern:'*_bcftools.vcf*' 
    publishDir "${params.outdir}vcf_pre-map2", mode: 'copy', pattern:'*_pre_bcftools.vcf*' 

    shell:
    '''
    #!/bin/bash
    ls -latr
    R1=!{base}
    echo "bamsize: !{bamsize}"
    #if [ -s !{} ]
    # More reliable way of checking bam size, because of aliases
    if (( !{bamsize} > 92 ))
    then
        # Parallelize pileup based on number of cores
        splitnum=$(($((!{id_ref_size}/!{task.cpus}))+1))
        perl !{VCFUTILS} splitchr -l $splitnum !{base}_mapped_ref_genome.fa.fai | \\
        #cat !{SPLITCHR} | \\
            xargs -I {} -n 1 -P !{task.cpus} sh -c \\
                "/usr/local/miniconda/bin/bcftools mpileup \\
                    -f !{base}_mapped_ref_genome.fa -r {} \\
                    --count-orphans \\
                    --no-BAQ \\
                    --max-depth 50000 \\
                    --max-idepth 500000 \\
                    --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \\
                !{base}.sorted.bam | /usr/local/miniconda/bin/bcftools call -A -m -M -Oz - > tmp.{}.vcf.gz"
        
        # Concatenate parallelized vcfs back together
        gunzip tmp*vcf.gz
        mv tmp.*:1-* \${R1}_catted.vcf
        for file in tmp*.vcf; do grep -v "#" $file >> \${R1}_catted.vcf; done
        cat \${R1}_catted.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' | /usr/local/miniconda/bin/bcftools norm -m -any > \${R1}_pre_bcftools.vcf
        
        # Make sure variants are majority variants for consensus calling
        /usr/local/miniconda/bin/bcftools filter -i '(DP4[0]+DP4[1]) < (DP4[2]+DP4[3]) && ((DP4[2]+DP4[3]) > 0) | (IMF > 0.5)' --threads !{task.cpus} \${R1}_pre_bcftools.vcf -o \${R1}.vcf
       # /usr/local/miniconda/bin/bcftools filter -e 'IMF < 0.5' \${R1}_pre2.vcf -o \${R1}.vcf
        # Index and generate consensus from vcf with majority variants
        /usr/local/miniconda/bin/bgzip \${R1}.vcf
        /usr/local/miniconda/bin/tabix \${R1}.vcf.gz 
        cat !{base}_mapped_ref_genome.fa | /usr/local/miniconda/bin/bcftools consensus \${R1}.vcf.gz > \${R1}.consensus.fa
        
        # Create coverage file from bam for whole genome, then pipe anything that has less than 6 coverage to bed file,
        # to be masked later
        /usr/local/miniconda/bin/bedtools genomecov \\
            -bga \\
            -ibam !{base}.sorted.bam \\
            -g !{base}_mapped_ref_genome.fa  \\
            | awk '\$4 < 6' | /usr/local/miniconda/bin/bedtools merge > \${R1}.mask.bed
        # Get rid of anything outside of the genome we care about, to prevent some sgrnas from screwing with masking
        awk '{ if(\$3 > 200 && \$2 < 29742) {print}}' \${R1}.mask.bed > a.tmp && mv a.tmp \${R1}.mask.bed
        # Mask refseq fasta for low coverage areas based on bed file
        /usr/local/miniconda/bin/bedtools maskfasta \\
            -fi !{base}_mapped_ref_genome.fa  \\
            -bed \${R1}.mask.bed \\
            -fo ref.mask.fasta
        
        # Align to refseq and unwrap fasta
        cat ref.mask.fasta \${R1}.consensus.fa > align_input.fasta
        /usr/local/miniconda/bin/mafft --auto --thread !{task.cpus} align_input.fasta > repositioned.fasta
        awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' repositioned.fasta > repositioned_unwrap.fasta
        # Trim ends and aligns masking of refseq to our consensus
        python3 !{TRIM_ENDS} \${R1}
        gunzip \${R1}.vcf.gz
        mv \${R1}.vcf \${R1}_bcftools.vcf
        sed -i 's/>.*/>!{base}.consensus/' \${R1}.consensus.fa
    else
       echo "Empty bam detected. Generating empty consensus fasta file..."
       touch \${R1}_bcftools.vcf
       touch \${R1}_pre_bcftools.vcf
    fi
    '''
}
}
if (params.singleEnd) {
process Final_Mapping {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"    
	errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(base), file("${base}_mapped_ref_genome.fa"), file("${base}.consensus.fa"), file("${base}_summary.csv"), val(bamsize), val(id),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt") from Consensus_Fasta_ch
    
    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}.consensus-final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map3_stats.txt"), file("${base}.mpileup"), file("${base}_final_summary.csv"), file("${base}.trimmed.fastq.gz"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt") into Mapping_Final_ch

    publishDir "${params.outdir}mpileup_map3", mode: 'copy', pattern:'*.mpileup*'
    publishDir "${params.outdir}bam_map3", mode: 'copy', pattern:'*.map3.sorted.bam*'
    publishDir "${params.outdir}sam_map3", mode: 'copy', pattern:'*_map3.sam*'
    publishDir "${params.outdir}consensus-final", mode: 'copy', pattern:'*.consensus-final*'
    // publishDir "${params.outdir}consensus-ivar-masked", mode: 'copy', pattern:'*.consensus.masked.fa*'
    publishDir "${params.outdir}txt_bbmap_map3_stats", mode: 'copy', pattern:'*_map3_stats.txt*'
    publishDir "${params.outdir}summary", mode: 'copy', pattern:'*_final_summary.csv*'

    script:

    """
    #!/bin/bash
    
    bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus.fa threads=${task.cpus} local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map3_stats.txt 2>&1
    samtools view -S -b ${base}_map3.sam > ${base}_map3.bam
    samtools sort -@ 4 ${base}_map3.bam > ${base}.map3.sorted.bam
    samtools index ${base}.map3.sorted.bam

    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}.consensus.fa \\
        --min-BQ 15 \\
        --output ${base}.mpileup \\
        ${base}.map3.sorted.bam
    cat ${base}.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}.consensus_final

    bedtools genomecov \\
        -bga \\
        -ibam ${base}.map3.sorted.bam \\
        -g ${base}_mapped_ref_genome.fa \\
        | awk '\$4 < 10' | bedtools merge > ${base}.mask.bed
    
    bedtools maskfasta \\
        -fi ${base}.consensus_final.fa \\
        -bed ${base}.mask.bed \\
        -fo ${base}.consensus.masked.fa

    sed -i 's/>.*/>${base}.ivar.masked.consensus/' ${base}.consensus.masked.fa
    sed -i 's/>.*/>${base}.ivar.consensus/' ${base}.consensus_final.fa

    awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length(\$0)}END{print l}' ${base}.consensus_final.fa > bases.txt
    num_bases=\$(awk 'FNR==2{print val,\$1}' bases.txt)

    seqkit -is replace -p "^n+|n+\$" -r "" ${base}.consensus_final.fa > ${base}.consensus-final.fa

    grep -v "^>" ${base}.consensus-final.fa | tr -cd N | wc -c > N.txt
    num_ns=\$(awk 'FNR==1{print val,\$1}' N.txt)
    echo "\$num_ns/\$num_bases*100" | bc -l > n_percent.txt
    percent_n=\$(awk 'FNR==1{print val,\$1}' n_percent.txt)
    printf ",\$num_bases" >> ${base}_summary.csv
    printf ",\$percent_n" >> ${base}_summary.csv

    cp ${base}_summary.csv ${base}_final_summary.csv
    
    """  
}
} else {
process Final_Mapping_PE {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"   
	errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(base), file("${base}_mapped_ref_genome.fa"), file("${base}.consensus.fa"), file("${base}_summary.csv"), val(bamsize), val(id),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt") from Consensus_Fasta_PE_ch
    
    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}.consensus-final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map3_stats.txt"), file("${base}.mpileup"), file("${base}_final_summary.csv"), file("${base}.trimmed.fastq.gz"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt") into Mapping_Final_PE_ch

    publishDir "${params.outdir}mpileup_map3", mode: 'copy', pattern:'*.mpileup*'
    publishDir "${params.outdir}bam_map3", mode: 'copy', pattern:'*.map3.sorted.bam*'
    publishDir "${params.outdir}sam_map3", mode: 'copy', pattern:'*_map3.sam*'
    publishDir "${params.outdir}consensus-final", mode: 'copy', pattern:'*.consensus-final*'
    // publishDir "${params.outdir}consensus-ivar-masked", mode: 'copy', pattern:'*.consensus.masked.fa*'
    publishDir "${params.outdir}txt_bbmap_map3_stats", mode: 'copy', pattern:'*_map3_stats.txt*'
    publishDir "${params.outdir}summary", mode: 'copy', pattern:'*_final_summary.csv*'

    script:

    """
    #!/bin/bash
    
    
    bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus.fa threads=${task.cpus} local=true -Xmx6g > ${base}_map3_stats.txt 2>&1
    samtools view -S -b ${base}_map3.sam > ${base}_map3.bam
    samtools sort -@ 4 ${base}_map3.bam > ${base}.map3.sorted.bam
    samtools index ${base}.map3.sorted.bam

    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}.consensus.fa \\
        --min-BQ 15 \\
        --output ${base}.mpileup \\
        ${base}.map3.sorted.bam
    cat ${base}.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}.consensus_final

    bedtools genomecov \\
        -bga \\
        -ibam ${base}.map3.sorted.bam \\
        -g ${base}_mapped_ref_genome.fa \\
        | awk '\$4 < 10' | bedtools merge > ${base}.mask.bed
    
    bedtools maskfasta \\
        -fi ${base}.consensus_final.fa \\
        -bed ${base}.mask.bed \\
        -fo ${base}.consensus.masked.fa

    sed -i 's/>.*/>${base}.ivar.masked.consensus/' ${base}.consensus.masked.fa
    sed -i 's/>.*/>${base}.ivar.consensus/' ${base}.consensus_final.fa

    awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length(\$0)}END{print l}' ${base}.consensus_final.fa > bases.txt
    num_bases=\$(awk 'FNR==2{print val,\$1}' bases.txt)

    seqkit -is replace -p "^n+|n+\$" -r "" ${base}.consensus_final.fa > ${base}.consensus-final.fa

    grep -v "^>" ${base}.consensus-final.fa | tr -cd N | wc -c > N.txt
    num_ns=\$(awk 'FNR==1{print val,\$1}' N.txt)
    echo "\$num_ns/\$num_bases*100" | bc -l > n_percent.txt
    percent_n=\$(awk 'FNR==1{print val,\$1}' n_percent.txt)
    printf ",\$num_bases" >> ${base}_summary.csv
    printf ",\$percent_n" >> ${base}_summary.csv

    cp ${base}_summary.csv ${base}_final_summary.csv
    
    """  
}
}
if (params.withSampleSheet) {
    if (params.singleEnd) {
    process Final_Processing {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"    
    errorStrategy 'retry'
    maxRetries 3

    input:
    file SAMPLE_LIST from SAMPLE_SHEET
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}.consensus-final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map3_stats.txt"), file("${base}.mpileup"), file("${base}_final_summary.csv"),file("${base}.trimmed.fastq.gz"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt") from Mapping_Final_ch    

    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}.consensus-final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map3_stats.txt"), file("${base}.mpileup"), file("${base}_sample_stats.csv"), file("${base}_summary.csv"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_sample_id.txt"), file("${base}_pcr_ct.txt"), file("${base}_method.txt") into Final_Processing_final_ch  

    publishDir "${params.outdir}summary", mode: 'copy', pattern:'*_summary.csv*'

    script:

    """
    #!/bin/bash

    R1=${base}
    NCBI_Name=\${R1:4:6}
    csvgrep -c sample_id -r \$NCBI_Name ${SAMPLE_LIST} > ${base}_sample_stats.csv

    csvcut -c 1 ${base}_sample_stats.csv > ${base}_sample_id.txt
    csvcut -c 2 ${base}_sample_stats.csv > ${base}_pcr_ct.txt
    csvcut -c 3 ${base}_sample_stats.csv > ${base}_method.txt

    sample_id=\$(cat ${base}_sample_id.txt | sed -n '2 p')
    pcr_ct=\$(cat ${base}_pcr_ct.txt | sed -n '2 p')
    method=\$(cat ${base}_method.txt | sed -n '2 p')

    reads_mapped=\$(cat ${base}_num_mapped.txt | tr -d " \t\n\r" | sed -n '1 p')
    num_trimmed=\$(cat ${base}_num_trimmed.txt | tr -d " \t\n\r" | sed -n '1 p')

    echo "\$reads_mapped/\$num_trimmed*100" | bc -l > reads-on-t_percent.txt
    Reads_On_Target=\$(awk 'FNR==1{print val,\$1}' reads-on-t_percent.txt)

    printf ",\$Reads_On_Target" >> ${base}_final_summary.csv
    printf ",\$pcr_ct" >> ${base}_final_summary.csv
    printf ",\$method" >> ${base}_final_summary.csv
    printf ",\$NCBI_Name" >> ${base}_final_summary.csv
    cp ${base}_final_summary.csv ${base}_summary.csv

    """
    }
    } else {
    process Final_Processing_PE {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"   
    errorStrategy 'retry'
    maxRetries 3

    input:
    file SAMPLE_LIST from SAMPLE_SHEET
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}.consensus-final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map3_stats.txt"), file("${base}.mpileup"), file("${base}_final_summary.csv"),file("${base}.trimmed.fastq.gz"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt") from Mapping_Final_PE_ch    

    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}.consensus-final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map3_stats.txt"), file("${base}.mpileup"), file("${base}_sample_stats.csv"), file("${base}_summary.csv"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_sample_id.txt"), file("${base}_pcr_ct.txt"), file("${base}_method.txt") into Final_Processing_PE_ch  

    publishDir "${params.outdir}summary", mode: 'copy', pattern:'*_summary.csv*'

    script:

    """
    #!/bin/bash

    R1=${base}
    NCBI_Name=\${R1:4:6}
    csvgrep -c sample_id -r \$NCBI_Name ${SAMPLE_LIST} > ${base}_sample_stats.csv

    csvcut -c 1 ${base}_sample_stats.csv > ${base}_sample_id.txt
    csvcut -c 2 ${base}_sample_stats.csv > ${base}_pcr_ct.txt
    csvcut -c 3 ${base}_sample_stats.csv > ${base}_method.txt

    sample_id=\$(cat ${base}_sample_id.txt | sed -n '2 p')
    pcr_ct=\$(cat ${base}_pcr_ct.txt | sed -n '2 p')
    method=\$(cat ${base}_method.txt | sed -n '2 p')

    reads_mapped=\$(cat ${base}_num_mapped.txt | tr -d " \t\n\r" | sed -n '1 p')
    num_trimmed=\$(cat ${base}_num_trimmed.txt | tr -d " \t\n\r" | sed -n '1 p')

    echo "\$reads_mapped/\$num_trimmed*100" | bc -l > reads-on-t_percent.txt
    Reads_On_Target=\$(awk 'FNR==1{print val,\$1}' reads-on-t_percent.txt)

    printf ",\$Reads_On_Target" >> ${base}_final_summary.csv
    printf ",\$pcr_ct" >> ${base}_final_summary.csv
    printf ",\$method" >> ${base}_final_summary.csv
    printf ",\$NCBI_Name" >> ${base}_final_summary.csv
    cp ${base}_final_summary.csv ${base}_summary.csv

    """
    }
    }
    }

if (params.withFastQC) {

    if (params.singleEnd) {
 /* FastQC
 *
 * Sequence read quality control analysis.
 */
process FastQC_SE {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"   
	errorStrategy 'retry'
    maxRetries 3

    input:
        file R1 from Trim_out_fastqc_SE

    output:
	file '*_fastqc.{zip,html}' into fastqc_results

    publishDir "${params.outdir}fastqc_results", mode: 'copy', pattern:'*_fastqc.{zip,html}*'  

    script:
    """
    #!/bin/bash
    fastqc --quiet --threads $task.cpus *.fastq.gz
    """
    }
} else {
process FastQC_PE {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"   
	errorStrategy 'retry'
    maxRetries 3

    input:
        file R1 from Trim_out_fastqc_PE

    output:
	file '*_fastqc.{zip,html}' into fastqc_results

    publishDir "${params.outdir}fastqc_results", mode: 'copy', pattern:'*_fastqc.{zip,html}*'  

    script:
    """
    #!/bin/bash
    fastqc --quiet --threads $task.cpus *.fastq.gz
    """
    }
}
}
