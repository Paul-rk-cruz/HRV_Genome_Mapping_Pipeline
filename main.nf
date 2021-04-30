#!/usr/bin/env nextflow

/*
========================================================================================
                  Rhinovirus Genome Mapping Pipeline v1.0
========================================================================================
 Github Repo:
 Greninger Lab
 
 Author:
 Paul RK Cruz <kurtisc@uw.edu>
 UW Medicine | Virology
 Department of Laboratory Medicine and Pathology
 University of Washington
 Created: April, 2021
 LICENSE: GNU
----------------------------------------------------------------------------------------

This pipeline was designed to run either single-end or paired end Next-Generation Sequencing reads to identify Human Rhinovirus complete genomes for analysis and Genbank submission.

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

    Setup Multifasta Reference:
    1. REFERENCE_FASTA (must be a multifasta containing concatenated full length RhV genome sequences (6-10K bp) formatted with accession numbers only)
        Current file: rhv_ref_db01_accession_only.fasta - 327 Human Rhinovirus Complete Genome Sequences courtesy of NCBI Genbank, 2021.
            source: https://www.ncbi.nlm.nih.gov/nucleotide/

    2. REFERENCE_FASTA_INDEX
        run:
             samtools faidx <reference.fasta>
        to create a multifasta index file.

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
version = '1.0'
def helpMsg() {
    log.info"""
	 __________________________________________________
     Human Rhinovirus Genome Mapping Pipeline :  Version ${version}
	__________________________________________________
    
	Pipeline Usage:

    To run the pipeline, enter the following in the command line:

        nextflow run FILE_PATH/HRV_Genome_Mapping_Pipeline/main.nf --reads PATH_TO_FASTQ --outdir PATH_TO_OUTPUT_DIR


    Valid CLI Arguments:
    REQUIRED:
      --reads                       Path to input fastq.gz folder).
      --outdir                      The output directory where the results will be saved
    OPTIONAL:
	  --helpMsg						Displays help message in terminal
      --singleEnd                   Specifies that the input fastq files are single end reads
	  --withFastQC					Runs a quality control check on fastq files
      --skipTrimming                Skips the fastq trimmming process

    """.stripIndent()
}
// Initialize parameters
params.helpMsg = false
params.virus_index = false
params.virus_fasta = false
REFERENCE_FASTA = file("${baseDir}/hrv_ref/hrv_ref_db01_accession_only.fa")
REFERENCE_FASTA_INDEX = file("${baseDir}/hrv_ref/hrv_ref_db01.fa.fai")
BBMAP_PATH="/Users/uwvirongs/Documents/KC/bbmap/"
// Show help msg
if (params.helpMsg){
    helpMsg()
    exit 0
}
params.withFastQC = false
params.skipTrim = false
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
// Default trimming options
params.trimmomatic_adapters_file_PE = "/Users/uwvirongs/miniconda3/share/trimmomatic-0.39-2/adapters/TruSeq2-PE.fa"
params.trimmomatic_adapters_file_SE = "/Users/uwvirongs/miniconda3/share/trimmomatic-0.39-2/adapters/TruSeq2-SE.fa"
params.trimmomatic_adapters_parameters = "2:30:10:1"
params.trimmomatic_window_length = "4"
params.trimmomatic_window_value = "20"
params.trimmomatic_mininum_length = "75"
trimmomatic_mininum_length = "75"
// Script Files
TRIM_ENDS=file("${baseDir}/scripts/trim_ends.py")
VCFUTILS=file("${baseDir}/scripts/vcfutils.pl")
SPLITCHR=file("${baseDir}/scripts/splitchr.txt")
FIX_COVERAGE = file("${baseDir}/scripts/fix_coverage.py")
// log files header
log.info "____________________________________________"
log.info " Human Rhinovirus Genome Mapping Pipeline :  v${version}"
log.info "____________________________________________"
def summary = [:]
summary['Fastq Files:']               = params.reads
summary['Read type:']           	  = params.singleEnd ? 'Single-End' : 'Paired-End'
summary['Virus Reference:']           = REFERENCE_FASTA
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
summary['Current directory path:']        = "$PWD"
summary['Working directory path:']         = workflow.workDir
summary['Output directory path:']          = params.outdir
summary['Pipeline directory path:']          = workflow.projectDir
if (params.singleEnd) {
summary['Trimmomatic adapters:'] = params.trimmomatic_adapters_file_SE
} else {
summary['Trimmomatic adapters:'] = params.trimmomatic_adapters_file_PE
}
summary['Trimmomatic adapter parameters:'] = params.trimmomatic_adapters_parameters
summary["Trimmomatic read length (minimum):"] = params.trimmomatic_mininum_length
summary['Configuration Profile:'] = workflow.profile
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "____________________________________________"
// Create channel for input reads.
// Import reads depending on single-end or paired-end
if(params.singleEnd == false) {
    // Check for R1s and R2s in input directory
    input_read_ch = Channel
        .fromFilePairs("${params.reads}*_R{1,2}*.fastq.gz")
        .ifEmpty { error "> Cannot located paired-end reads in: ${params.reads}.\n> Please enter a valid file path." }
        .map { it -> [it[0], it[1][0], it[1][1]]}
} else {
    // input: *.gz
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
 * Trim Reads
 * 
 * Processing: Trim adaptors and repetitive bases from sequence reads and remove low quality sequence reads.
 */
// if(!params.skipTrimming) {
if (params.singleEnd) {
	process Trim_Reads_SE {
    errorStrategy 'retry'
    maxRetries 3

    input:
        file R1 from input_read_ch
        val trimmomatic_mininum_length

    output: 
        tuple env(base),file("*.trimmed.fastq.gz") into Trim_out_map1_ch, Trim_out_map2_ch, Trim_out_fastqc_SE

    publishDir "${params.outdir}fastq_trimmed", mode: 'copy',pattern:'*.trimmed.fastq*'


    script:
    """
    #!/bin/bash
    base=`basename ${R1} ".fastq.gz"`

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
        val trimmomatic_mininum_length
    output: 
        tuple env(base),file("*.trimmed.fastq.gz") into Trim_out_map1_ch, Trim_out_map2_ch, Trim_out_fastqc_PE
        // tuple val(base),file("${base}_results.csv") into Results_trimmed_ch

    publishDir "${params.outdir}fastq_trimmed", mode: 'copy',pattern:'*.trimmed.fastq*'
    
    script:
    """
    #!/bin/bash

    trimmomatic PE -threads ${task.cpus} ${R1} ${R2} ${base}.R1.paired.fastq.gz ${base}.R1.unpaired.fastq.gz ${base}.R2.paired.fastq.gz ${base}.R2.unpaired.fastq.gz \
	ILLUMINACLIP:${params.trimmomatic_adapters_file_PE}:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:${params.trimmomatic_mininum_length}


    """
}
}
// } else {
//    input_read_ch
//        .set {Trim_out_map1_ch}
// }
/*
 * Map sequence reads to HRV Genomes using BBMap.
 */
process Mapping {
    errorStrategy 'retry'
    maxRetries 3

    input: 
        tuple val(base), file("${base}.trimmed.fastq.gz") from Trim_out_map1_ch
        file REFERENCE_FASTA

    output:
        tuple val(base), file("${base}_map1.sam")into Sam_first_mapping_ch
        tuple val(base), file("${base}_map2.sam")into Aligned_sam_ch
        tuple val(base), file("${base}_most_mapped_ref.txt") into Mapped_Ref_Id_ch, Mapped_Ref_Cons_Id_ch
        tuple val(base), file("${base}_most_mapped_ref_size.txt") into Mapped_Ref_Id_Size_ch       
        tuple val(base), file("${base}_most_mapped_ref_size_out.txt"),env(id_ref_size) into Mapped_Ref_Id_perc_ch           
        tuple val(base), file("${base}_idxstats.txt") into Indx_stats_Ch
        tuple val(base), file("${base}_mapped_ref_genome.fa"),env(id) into Mapped_Ref_Gen_ch, Mapped_Ref_Gen_map2_ch, Mapped_Ref_Gen_Cons_ch
        tuple val(base), file("${base}_map1_bbmap_out.txt")into Bbmap_map1_bbmap_txt_ch
        tuple val(base), file("${base}_map2_bbmap_out.txt")into Bbmap_map2_bbmap_txt_ch
        tuple val(base), file("${base}_map1_stats.txt") into Bbmap_map1_stats_ch
        tuple val(base), file("${base}_map2_stats.txt") into Bbmap_map2_stats_ch
        tuple val(base), file("${base}_map1_histogram.txt") into BBmap_map1_hist_ch
        tuple val(base), file("${base}_map2_histogram.txt") into BBmap_map2_hist_ch
        tuple val(base), file("${base}_mapped_ref_genome.fa.fai") into Mapped_Ref_Gen_Index_fai_Cons_ch
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

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${REFERENCE_FASTA} threads=8 covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false -Xmx6g > ${base}_map1_stats.txt 2>&1
    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    samtools faidx ${REFERENCE_FASTA} \$id > ${base}_mapped_ref_genome.fa
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map2.sam ref=${base}_mapped_ref_genome.fa threads=8 covstats=${base}_map2_bbmap_out.txt covhist=${base}_map2_histogram.txt local=true interleaved=false -Xmx6g > ${base}_map2_stats.txt 2>&1
    samtools faidx ${base}_mapped_ref_genome.fa
    awk 'NR == 2 || \$5 > max {number = \$3; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref_size_out.txt
    id_ref_size=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref_size_out.txt)
    echo \$id_ref_size >> ${base}_most_mapped_ref_size.txt

    """
}
/*
 * Convert BAM to coordinate sorted BAM
 */
 // Step 1. Convert Sam to Bam
 // Step 2. Sort Bam file by coordinates
 // Step 3. Generate Statistics about Bam file
process Sort_Bam {
	errorStrategy 'retry'
    maxRetries 3

    input: 
    tuple val(base), file("${base}_map2.sam") from Aligned_sam_ch

    output:
    tuple val(base), file("${base}.bam") into Aligned_bam_ch, Bam_ch
    tuple val(base), file("${base}.sorted.bam"),env(bamsize) into Sorted_bam_ch, Sorted_Cons_Bam_ch
    tuple val(base), file("${base}_flagstats.txt") into Flagstats_ch
    tuple val(base), file("${base}_bam_size.txt") into Bam_Size_ch

    publishDir "${params.outdir}bam", mode: 'copy', pattern:'*.bam*'
    publishDir "${params.outdir}bam_sorted", mode: 'copy', pattern:'*.sorted.bam*'  
    publishDir "${params.outdir}txt_bam_flagstats", mode: 'copy', pattern:'*_flagstats.txt*'  
    publishDir "${params.outdir}bam_size", mode: 'copy', pattern:'*_bam_size.txt*'  

    script:
    """
    #!/bin/bash
    samtools view -S -b ${base}_map2.sam > ${base}.bam
    samtools sort -@ ${task.cpus} ${base}.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools flagstat ${base}.sorted.bam > ${base}_flagstats.txt
    bedtools genomecov -d -ibam ${base}.sorted.bam > ${base}_coverage.txt
    meancoverage=\$(cat ${base}_coverage.txt | awk '{sum+=\$3} END { print sum/NR}')
    bamsize=\$((\$(wc -c ${base}.sorted.bam | awk '{print \$1'})+0)) 
    echo \$bamsize >> ${base}_bam_size.txt
    """
}
/*
 * Call variants & Generate Consensus
  */
 process Generate_Consensus {
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(base), file("${base}.sorted.bam"),val(bamsize) from Sorted_Cons_Bam_ch
    tuple val(base), file("${base}_mapped_ref_genome.fa"),val(id) from Mapped_Ref_Gen_Cons_ch
    tuple val(base), file("${base}_mapped_ref_genome.fa.fai") from Mapped_Ref_Gen_Index_fai_Cons_ch
    tuple val(base), file("${base}_most_mapped_ref.txt") from Mapped_Ref_Cons_Id_ch
    tuple val(base), file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size) from Mapped_Ref_Id_perc_ch    
    file TRIM_ENDS
    file FIX_COVERAGE
    file VCFUTILS

    output:
    tuple val(base), file("${base}_consensus.fa") into Consensus_Fasta_ch
    tuple val(base), file("${base}_bcftools.vcf") into Vcf_ch
    tuple val(base), val(bamsize), file("${base}_pre_bcftools.vcf") into Vcf_Pre_ch
    file(INDEX_FILE)

    publishDir "${params.outdir}consensus", mode: 'copy', pattern:'*_consensus.fa*'  
    publishDir "${params.outdir}vcf", mode: 'copy', pattern:'*_bcftools.vcf*' 
    publishDir "${params.outdir}vcf-pre", mode: 'copy', pattern:'*_pre_bcftools.vcf*' 

	script:

    shell:
    '''
    #!/bin/bash
    ls -latr

    R1=!{base}
    id_ref=!{id}
    echo "bamsize: !{bamsize}"

    #if [ -s !{} ]
    # More reliable way of checking bam size, because of aliases
    if (( !{bamsize} > 92 ))
    then
        # Parallelize pileup based on number of cores
            splitnum=$(($((!{id_ref_size}/!{task.cpus}))+1))
            perl !{VCFUTILS} splitchr -l $splitnum !{base}_mapped_ref_genome.fa.fai | \\
                xargs -I {} -n 1 -P !{task.cpus} sh -c \\
                    "bcftools mpileup \\
                        -f !{base}_mapped_ref_genome.fa -r {} \\
                        --count-orphans \\
                        --no-BAQ \\
                        --max-depth 50000 \\
                        --max-idepth 500000 \\
                        --annotate FORMAT/AD,FORMAT/ADF,FORMAT/ADR,FORMAT/DP,FORMAT/SP,INFO/AD,INFO/ADF,INFO/ADR \\
                    !{base}.sorted.bam | bcftools call -A -m -Oz - > tmp.{}.vcf.gz"
            
            # Concatenate parallelized vcfs back together
            gunzip tmp*vcf.gz
            mv tmp.\${id_ref}\\:1-* \${R1}_catted.vcf
            for file in tmp*.vcf; do grep -v "#" $file >> \${R1}_catted.vcf; done

            cat \${R1}_catted.vcf | awk '$1 ~ /^#/ {print $0;next} {print $0 | "sort -k1,1 -k2,2n"}' | bcftools norm -m -any > \${R1}_pre_bcftools.vcf

            # Make sure variants are majority variants for consensus calling
            bcftools filter -i '(DP4[0]+DP4[1]) < (DP4[2]+DP4[3]) && ((DP4[2]+DP4[3]) > 0)' --threads !{task.cpus} \${R1}_pre_bcftools.vcf -o \${R1}_pre2.vcf
            bcftools filter -e 'IMF < 0.5' \${R1}_pre2.vcf -o \${R1}.vcf

            # Index and generate consensus from vcf with majority variants
            bgzip \${R1}.vcf
            tabix \${R1}.vcf.gz 
            cat !{base}_mapped_ref_genome.fa | bcftools consensus \${R1}.vcf.gz > \${R1}.consensus.fa
    fi
    '''
}
/*
 * Variant Calling
 */
 // Step 1: Calculate the read coverage of positions in the genome
 // Step 2: Detect the single nucleotide polymorphisms (SNPs)
 // Step 3: Filter and report the SNP variants in variant calling format (VCF)
 // VIEW RESULTS:   less -S ${base}_final_variants.vcf
// process Variant_Calling {
// 	errorStrategy 'retry'
//     maxRetries 3

// 	input:
//     tuple val(base), file("${base}.sorted.bam") from Sorted_bam_ch
//     tuple val(base), file("${base}_mapped_ref_genome.fasta") from Mapped_Ref_Gen_ch
//     file REFERENCE_FASTA
//     file REFERENCE_FASTA_INDEX

// 	output:
//     tuple val(base), file("${base}.pileup") into Mpileup_ch
//     tuple val(base), file("${base}_lowfreq.vcf") into Low_Freq_ch      
//     tuple val(base), file("${base}_majority.vcf") into Majority_Freq_ch  

//     publishDir "${params.outdir}vcf_majority_freq", mode: 'copy', pattern:'*_majority.vcf*'  
//     publishDir "${params.outdir}vcf_low_freq", mode: 'copy', pattern:'*_lowfreq.vcf*'  
//     publishDir "${params.outdir}mpileup", mode: 'copy', pattern:'*.pileup*' 

// 	script:

// 	"""
//     #!/bin/bash

//   samtools mpileup -A -B -d 20000 -Q 15 -f ${base}_mapped_ref_genome.fasta ${base}.sorted.bam > ${base}.pileup
//   varscan mpileup2cns ${base}.pileup --min-var-freq 0.02 --min-avg-qual 15 --p-value 0.99 --output-vcf 1 > ${base}_lowfreq.vcf
//   varscan mpileup2cns ${base}.pileup --min-var-freq 0.4 --min-avg-qual 15 --p-value 0.99 --min-coverage 10 --output-vcf 1 > ${base}_majority.vcf

// 	"""
// }
/*
 * Consensus
 *
 * Consensus generation using sorted Bam, final_variants.vcf, and mapped_ref_genome.
 */
// process Consensus {
// 	errorStrategy 'retry'
//     maxRetries 3

//     input:
//     tuple val(base), file("${base}_majority.vcf") from Majority_Freq_ch  
//     tuple val(base), file("${base}.sorted.bam") from Sorted_Cons_Bam_ch
//     tuple val(base), file("${base}_mapped_ref_genome.fasta") from Mapped_Ref_Gen_Cons_ch
//     tuple val(base), file("${base}_most_mapped_ref.txt") from Mapped_Ref_Final_Cons_Id_ch

//     output:
//     tuple val(base), file("${base}_consensus.fasta") into Consensus_Fasta_ch, Consensus_Fasta_Processing_ch
//     tuple val(base), file("${base}_consensus_masked.fasta") into Consensus_fasta_Masked_ch
//     tuple val(base), file("${base}_bed4mask.bed") into Consensus_bed4mask_ch
//     tuple val(base), file("${base}_consensus_final.fasta") into Consensus_fasta_Complete_ch
    
//     publishDir "${params.outdir}consensus_fasta_files", mode: 'copy', pattern:'*_consensus.fasta*'  
//     publishDir "${params.outdir}consensus_masked_fasta_files", mode: 'copy', pattern:'*_consensus_masked.fasta*'  
//     publishDir "${params.outdir}bed4mask_bed_files", mode: 'copy', pattern:'*_bed4mask.bed*'  
//     publishDir "${params.outdir}consensus_final", mode: 'copy', pattern:'*_consensus_final.fasta*' 

//     script:

//     """
//     #!/bin/bash

//     bgzip -c ${base}_majority.vcf > ${base}_majority.vcf.gz
//     bcftools index ${base}_majority.vcf.gz
//     cat ${base}_mapped_ref_genome.fasta | bcftools consensus ${base}_majority.vcf.gz > ${base}_consensus.fasta
//     bedtools genomecov -bga -ibam ${base}.sorted.bam -g ${base}_mapped_ref_genome.fasta | awk '\$4 < 20' | bedtools merge > ${base}_bed4mask.bed
//     bedtools maskfasta -fi ${base}_consensus.fasta -bed ${base}_bed4mask.bed -fo ${base}_consensus_masked.fasta
//     id=\$(awk '{print \$1}' ${base}_most_mapped_ref.txt)
//     seqkit replace -p "\$id" -r '${base}' ${base}_consensus.fasta > ${base}_consensus_final.fasta

//     """
// }
if (params.withFastQC) {

    if (params.singleEnd) {
 /* FastQC
 *
 * Sequence read quality control analysis.
 */
process FastQC_SE {
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

    

    """
    }
} else {
process FastQC_PE {
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
