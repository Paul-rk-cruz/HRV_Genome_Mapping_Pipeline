#!/usr/bin/env nextflow

/*
========================================================================================
                 Human Respiratory Virus Pipeline v1.4
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
 Updated: August 24, 2021
 LICENSE: GNU
----------------------------------------------------------------------------------------
Human Respiratory Virus Pipeline was designed to run either single-end or paired end Illumina Next-Generation-Sequencing (NGS) sequence for Human respiratory virus discovery, analysis, and Genbank submission.
PIPELINE OVERVIEW:
 - 1. : Trim Reads
 		-Trimmomatic - sequence read trimming of adaptors and low quality reads.
 - 2. : Genome Mapping
 		-BBMap - align to MultiFasta Reference Virus Genome.
 		-Samtools - SAM and BAM file processing.
 - 3. : Reference Fasta Selection
 		-Selects the closest reference genome.  
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
 - 8. : Summary Report Generation
        Generates a run report summary.               
 - 9. : FastQC
 		-Sequence read quality control analysis.

Dependencies:

HRV-Docker includes all dependencies. Currently (7/2021), Mapping step requires local dependencies. Please see docker for dependencies required for mapping and viral annotation steps.
    
Setup Multifasta References:

1. Multifasta references containing Viral genome sequences formatted with accession numbers only.
   Current supported virus`s: Rhinovirus, Human Coronavirus, Influenza B, HPIV3, and HPV.

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

    Single-end Illumina Reads

    nextflow run ./HRV_Pipeline/main.nf --reads './HRV_Pipeline/example/' --outdir './HRV_Pipeline/example/output/' --withMetadata './HRV_Pipeline/example/example_metadata.csv' --withSerotype --singleEnd -resume

    Paired-end Illumina Reads

    nextflow run ./HRV_Pipeline/main.nf --reads './HRV_Pipeline/example/' --outdir './HRV_Pipeline/example/output/' --withMetadata './HRV_Pipeline/example/example_metadata.csv' --withSerotype -resume


 ----------------------------------------------------------------------------------------
*/
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                DISPLAY HELP MSG                    */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
// Pipeline version
version = '1.4'
def helpMsg() {
    log.info"""
	 _______________________________________________________________________________
     Human Rhinovirus Genome Mapping Pipeline :  Version ${version}
	________________________________________________________________________________
    
	Pipeline Usage:
    To run the pipeline, enter the following in the command line:
        nextflow run FILE_PATH/HRV_Genome_Mapping_Pipeline/main.nf --reads PATH_TO_FASTQ --outdir PATH_TO_OUTPUT_DIR
    Valid CLI Arguments:
    REQUIRED:
      --reads                       Path to input fastq.gz folder).
      --outdir                      The output directory where the results will be saved
    OPTIONAL:
      --withMetadata                Adds Metadata information to Final Run Report Summary
      --withVapid                   Annotate the resulting consensus fasta for GenBank submission         
      --ref_rv                      Overwrite set multi-fasta Rhinovirus reference file
      --ref_hcov                    Overwrite set multi-fasta Human Coronavirus reference file
      --ref_hpv                     Overwrite set multi-fasta HPV reference file
      --ref_inflb                   Overwrite set multi-fasta Influenza B reference file
      --ref_hpiv3                   Overwrite set multi-fasta HPIV3 reference file
	  --helpMsg						Displays help message in terminal
      --singleEnd                   Specifies that the input fastq files are single end reads
	  --withFastQC					Runs a quality control check on fastq files
      --skipTrimming                Skips the fastq trimmming process   
    """.stripIndent()
}
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*              CONFIGURATION VARIABLES               */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
// Make sure outdir path ends with trailing slash
if (!params.outdir.endsWith("/")){
   params.outdir = "${params.outdir}/"
}
// Make sure reads path ends with trailing slash
if (!params.reads.endsWith("/")){
   params.reads = "${params.outdir}/"
}
params.helpMsg = false
params.virus_index = false
params.virus_fasta = false
params.withFastQC = false
params.skipTrim = false
params.reads = false
params.singleEnd = false
params.ADAPTERS = false
params.withMetadata = false
params.withSerotype = false
params.withVapid = false
// Trimmomatic Paths and variables
params.ADAPTERS_SE = file("${baseDir}/adapters/TruSeq2-SE.fa")
params.ADAPTERS_EE = file("${baseDir}/adapters/TruSeq2-PE.fa")
ADAPTERS_SE = file("${baseDir}/adapters/TruSeq2-SE.fa")
ADAPTERS_PE = file("${baseDir}/adapters/TruSeq2-PE.fa")
vapid_rhinovirus_sbt = file("${baseDir}/vapid/rhinovirus.sbt")
params.SETTING = "2:30:10:1:true"
SETTING = "2:30:10:1:true"
params.LEADING = "3"
LEADING = "3"
params.TRAILING = "3"
TRAILING = "3"
params.SWINDOW = "4:20"
SWINDOW = "4:20"
// Script Files
if(params.withMetadata != false) {
    METADATA = file(params.withMetadata)
}
if(params.withVapid != false) {
    vapid_python = file("${baseDir}/vapid/vapid.py")
    vapid_python3 = file("${baseDir}/vapid/vapid3.py")
    VAPID_DB_ALL_1 = file("${baseDir}/vapid/all_virus.fasta")
    VAPID_DB_ALL_2 = file("${baseDir}/vapid/all_virus.fasta.ndb")
    VAPID_DB_ALL_3 = file("${baseDir}/vapid/all_virus.fasta.nhr")
    VAPID_DB_ALL_4 = file("${baseDir}/vapid/all_virus.fasta.nin")
    VAPID_DB_ALL_5 = file("${baseDir}/vapid/all_virus.fasta.not")
    VAPID_DB_ALL_6 = file("${baseDir}/vapid/all_virus.fasta.nsq")
    VAPID_DB_ALL_7 = file("${baseDir}/vapid/all_virus.fasta.ntf")
    VAPID_DB_ALL_8 = file("${baseDir}/vapid/all_virus.fasta.nto")
    tbl2asn = file("${baseDir}/vapid/tbl2asn")
}
if(params.withSerotype != false) {
    // Rhinovirus VP1 database files
    BLAST_DB_VP1_1 = file("${baseDir}/blast_db/VP1_164_annotated_nospaces.fasta")
    BLAST_DB_VP1_2 = file("${baseDir}/blast_db/VP1_164_annotated_nospaces.fasta.ndb")
    BLAST_DB_VP1_3 = file("${baseDir}/blast_db/VP1_164_annotated_nospaces.fasta.nhr")
    BLAST_DB_VP1_4 = file("${baseDir}/blast_db/VP1_164_annotated_nospaces.fasta.nin")
    BLAST_DB_VP1_5 = file("${baseDir}/blast_db/VP1_164_annotated_nospaces.fasta.not")
    BLAST_DB_VP1_6 = file("${baseDir}/blast_db/VP1_164_annotated_nospaces.fasta.nsq")
    BLAST_DB_VP1_7 = file("${baseDir}/blast_db/VP1_164_annotated_nospaces.fasta.ntf")
    BLAST_DB_VP1_8 = file("${baseDir}/blast_db/VP1_164_annotated_nospaces.fasta.nto")
    // HPV database files
    BLAST_DB_ALL_1hpv = file("${baseDir}/blast_db/hpv.fasta")
    BLAST_DB_ALL_2hpv = file("${baseDir}/blast_db/hpv.fasta.ndb")
    BLAST_DB_ALL_3hpv = file("${baseDir}/blast_db/hpv.fasta.nhr")
    BLAST_DB_ALL_4hpv = file("${baseDir}/blast_db/hpv.fasta.nin")
    BLAST_DB_ALL_5hpv = file("${baseDir}/blast_db/hpv.fasta.nog")
    BLAST_DB_ALL_6hpv = file("${baseDir}/blast_db/hpv.fasta.nos")
    BLAST_DB_ALL_7hpv = file("${baseDir}/blast_db/hpv.fasta.not")
    BLAST_DB_ALL_8hpv = file("${baseDir}/blast_db/hpv.fasta.nsq")
    BLAST_DB_ALL_9hpv = file("${baseDir}/blast_db/hpv.fasta.ntf")    
    BLAST_DB_ALL_10hpv = file("${baseDir}/blast_db/hpv.fasta.nto")  
    // All BLAST db files for respiratory viruses recognized by this pipeline
    BLAST_DB_ALL_1 = file("${baseDir}/blast_db/all_ref.fasta")
    BLAST_DB_ALL_2 = file("${baseDir}/blast_db/all_ref.fasta.ndb")
    BLAST_DB_ALL_3 = file("${baseDir}/blast_db/all_ref.fasta.nhr")
    BLAST_DB_ALL_4 = file("${baseDir}/blast_db/all_ref.fasta.nin")
    BLAST_DB_ALL_5 = file("${baseDir}/blast_db/all_ref.fasta.nog")
    BLAST_DB_ALL_6 = file("${baseDir}/blast_db/all_ref.fasta.nos")
    BLAST_DB_ALL_7 = file("${baseDir}/blast_db/all_ref.fasta.not")
    BLAST_DB_ALL_8 = file("${baseDir}/blast_db/all_ref.fasta.nsq")
    BLAST_DB_ALL_9 = file("${baseDir}/blast_db/all_ref.fasta.ntf")    
    BLAST_DB_ALL_10 = file("${baseDir}/blast_db/all_ref.fasta.nto")  
}
// File paths
BLAST_DB_VP1 = file("${baseDir}/blast_db/VP1_164_annotated_nospaces.fasta")
BLAST_DB_ALL = file("${baseDir}/blast_db/allref.fasta")
BBMAP_PATH="/Users/greningerlab/Documents/bbmap/"
params.MINLEN = "35"
MINLEN = "35"
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*            VIRAL MULTI-FASTA REFERENCES            */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
// Setup MULTIFASTA Reference parameters.
params.ref_rv = false
params.ref_hcov = false
params.ref_hpv = false
params.ref_inflb = false
params.ref_hpiv3 = false
// Setup MULTIFASTA Reference file paths for OPTIONAL-override of set file paths.
Reference_hpv_14 = file("${baseDir}/hrv_ref/hrv_ref_hpv_14.fa")
// Rhinovirus
if(params.ref_rv != false) {
    Reference_rv = file(params.ref_rv)
} else {
Reference_rv=file("${baseDir}/hrv_ref/hrv_ref_rhinovirus.fa")
}
// Human Coronavirus
if(params.ref_hcov != false) {
    Reference_hcov = file(params.ref_hcov)
} else {
Reference_hcov=file("${baseDir}/hrv_ref/hrv_ref_hcov.fa")
}
// HPV
if(params.ref_hpv != false) {
    Reference_hpv = file(params.ref_hpv)
} else {
Reference_hpv=file("${baseDir}/hrv_ref/hrv_ref_hpv.fa")
}
// Influenza B
if(params.ref_inflb != false) {
    Reference_inflb = file(params.ref_inflb)
} else {
Reference_inflb=file("${baseDir}/hrv_ref/hrv_ref_Influenza_B.fa")
}
// HPIV3 - Human Parainfluenza Virus 3
if(params.ref_hpiv3 != false) {
    Reference_hpiv3 = file(params.ref_hpiv3)
} else {
Reference_hpiv3=file("${baseDir}/hrv_ref/hrv_ref_hpiv3.fa")
}
// Show help msg
if (params.helpMsg){
    helpMsg()
    exit 0
}
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                  SET UP CHANNELS                   */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
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
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*              WORKFLOW DISPLAY HEADER               */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
def hrvheader() {
    
    return """
    """.stripIndent()
}
// log files header
// log.info hrvheader()
log.info "_______________________________________________________________________________"
log.info " Human Respiratory Virus Pipeline :  v${version}"
log.info "_______________________________________________________________________________"
def summary = [:]
summary['Configuration Profile:'] = workflow.profile
summary['Current directory path:']        = "$PWD"
summary['HRV Pipeline directory path:']          = workflow.projectDir
summary['Input directory path:']               = params.reads
summary['Output directory path:']          = params.outdir
summary['Work directory path:']         = workflow.workDir
summary['Metadata:']               = params.withMetadata ? 'True' : 'False'
summary['Serotyping:']               = params.withSerotype ? 'True' : 'False'
summary['Fastq type:']           	  = params.singleEnd ? 'Single-End' : 'Paired-End'
if(workflow.revision) summary['Pipeline Release'] = workflow.revision
if (params.singleEnd) {
summary['Trimmomatic adapters:'] = params.ADAPTERS_SE
} else {
summary['Trimmomatic adapters:'] = params.ADAPTERS_PE
}
summary["Trimmomatic read length (minimum):"] = params.MINLEN
summary["Trimmomatic Setting:"] = params.SETTING
summary["Trimmomatic Sliding Window:"] = params.SWINDOW
summary["Trimmomatic Leading:"] = params.LEADING
summary["Trimmomatic Trailing:"] = params.TRAILING
log.info summary.collect { k,v -> "${k.padRight(21)}: $v" }.join("\n")
log.info "_______________________________________________________________________________"
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                WORKFLOW PROCESSES                  */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*
 * Trim Reads
 * 
 * Processing: Trim adaptors and repetitive bases from sequence reads and remove low quality sequence reads.
 */
// if(!params.skipTrimming) {
if (params.singleEnd) {
	process Trim_Reads {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"
    errorStrategy 'retry'
    maxRetries 3

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
    ILLUMINACLIP:${ADAPTERS_SE}:${SETTING} LEADING:${LEADING} TRAILING:${TRAILING} SLIDINGWINDOW:${SWINDOW} MINLEN:${MINLEN}
    num_untrimmed=\$((\$(gunzip -c ${R1} | wc -l)/4))
    num_trimmed=\$((\$(gunzip -c \$base'.trimmed.fastq.gz' | wc -l)/4))
    printf "\$num_trimmed" >> ${R1}_num_trimmed.txt
    percent_trimmed=\$((100-\$((100*num_trimmed/num_untrimmed))))
    echo Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Reference_Genome,Reference_Length,Mapped_Reads,Percent_Ref_Coverage,Min_Coverage,Mean_Coverage,Max_Coverage,Bam_Size,Consensus_Length,Percent_N,%_Reads_On_Target, PCR_CT,Method, NCBI_Name, Serotype, Nomenclature, Reference_Name, Reference_Genome, Biosample_name, Biosample_accession, SRA_Accession, Release_date, Bioproject> \$base'_summary.csv'
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
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${R1}_num_trimmed.txt"),file("*summary.csv") into Trim_out_PE
        tuple val(base), file(R1),file(R2),file("${base}.R1.paired.fastq.gz"), file("${base}.R2.paired.fastq.gz"),file("${base}.R1.unpaired.fastq.gz"), file("${base}.R2.unpaired.fastq.gz") into Trim_out_fastqc_PE
        tuple val(base), file("${base}.trimmed.fastq.gz") into Trim_out_ch3

    publishDir "${params.outdir}fastq_trimmed", mode: 'copy',pattern:'*.trimmed.fastq*'
    
    script:
    """
    #!/bin/bash
    /usr/local/miniconda/bin/trimmomatic PE -threads ${task.cpus} ${R1} ${R2} ${base}.R1.paired.fastq.gz ${base}.R1.unpaired.fastq.gz ${base}.R2.paired.fastq.gz ${base}.R2.unpaired.fastq.gz \
    ILLUMINACLIP:${ADAPTERS_PE}:${SETTING} LEADING:${LEADING} TRAILING:${TRAILING} SLIDINGWINDOW:${SWINDOW} MINLEN:${MINLEN}
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
    
    echo Sample_Name,Raw_Reads,Trimmed_Paired_Reads,Trimmed_Unpaired_Reads,Total_Trimmed_Reads, Percent_Trimmed,Reference_Genome,Reference_Length,Mapped_Reads,Percent_Ref_Coverage,Min_Coverage,Mean_Coverage,Max_Coverage,Bam_Size,Consensus_Length,Percent_N,%_Reads_On_Target, PCR_CT,Method, NCBI_Name, Serotype, Nomenclature, Reference_Name, Reference_Genome, Biosample_name, Biosample_accession, SRA_Accession, Release_date, Bioproject> \$base'_summary.csv'
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
    echo true

    input: 
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_summary.csv") from Trim_out_SE
        file Reference_rv
        file Reference_hcov
        file Reference_hpv
        file Reference_hpv_14        
        file Reference_inflb
        file Reference_hpiv3

    output:
        tuple val(base), file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"), file("${base}_summary2.csv"), file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),env(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),env(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") into Everything_ch
        tuple val(base), file("${base}_map1_histogram.txt"),file("${base}_map2_histogram.txt") into BBmap_map1_hist_ch
        tuple val (base), file("*") into Dump_map1_ch

    publishDir "${params.outdir}all_ref", mode: 'copy', pattern:'*all_ref.sam*'
    publishDir "${params.outdir}all_ref", mode: 'copy', pattern:'*all_ref*'    
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

    cat ${Reference_hpv} ${Reference_rv} ${Reference_inflb} ${Reference_hcov} ${Reference_hpiv3} > ${base}_all_ref.fa

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_all_ref.sam ref=${base}_all_ref.fa threads=${task.cpus} covstats=${base}_all_ref_bbmap_out.txt covhist=${base}_all_ref_histogram.txt local=true interleaved=false -Xmx10g > ${base}_all_ref_stats.txt 2>&1

    samtools view -S -b ${base}_all_ref.sam > ${base}_all_ref.bam
    samtools sort -@ 4 ${base}_all_ref.bam > ${base}_all_ref.sorted.bam
    samtools index ${base}_all_ref.sorted.bam
    samtools idxstats ${base}_all_ref.sorted.bam > ${base}_all_ref_idxstats.txt
    
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_all_ref_bbmap_out.txt > ${base}_all_ref_id.txt
    all_ref_id=\$(awk '{print \$1}' ${base}_all_ref_id.txt)

    grep -B 0 ">" ${Reference_rv} | tr -d ">" > ${base}_rv_ids.txt
    grep -B 0 ">" ${Reference_hpv} | tr -d ">" > ${base}_hpv_ids.txt
    grep -B 0 ">" ${Reference_inflb} | tr -d ">" > ${base}_inbflb_ids.txt
    grep -B 0 ">" ${Reference_hcov} | tr -d ">" > ${base}_hcov_ids.txt
    grep -B 0 ">" ${Reference_hpiv3} | tr -d ">" > ${base}_hpiv3.txt

    # Rhinovirus
    if grep -q \$all_ref_id "${base}_rv_ids.txt";
    then
    echo "< Accession found in Rhinovirus multifasta file. hrv_ref_rhinovirus.fa will be used for mapping."


    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${Reference_rv} threads=${task.cpus} covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map1_stats.txt 2>&1
    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${Reference_rv} \$id > ${base}_mapped_ref_genome.fa
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map2.sam ref=${base}_mapped_ref_genome.fa threads=${task.cpus} covstats=${base}_map2_bbmap_out.txt covhist=${base}_map2_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map2_stats.txt 2>&1
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


    # HPV
    elif grep -q \$all_ref_id "${base}_hpv_ids.txt";
    then
    echo "< Accession found in HPV multifasta file. hrv_ref_hpv.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_hpv_all_map1.sam ref=${Reference_hpv} threads=${task.cpus} covstats=${base}_hpv_all_map1_bbmap_out.txt covhist=${base}_hpv_all_map1_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_hpv_all_map1_stats.txt 2>&1

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_hpv_14_map1.sam ref=${Reference_hpv_14} threads=${task.cpus} covstats=${base}_hpv_14_map1_bbmap_out.txt covhist=${base}_hpv_14_map1_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_hpv_14_map1_stats.txt 2>&1

    awk '{print \$5, \$1}' ${base}_hpv_all_map1_bbmap_out.txt  | sort -rn | head -n 2 > ${base}_most_mapped_ref.txt
    sed -n '2p' < ${base}_most_mapped_ref.txt > ${base}_ref2_percent.txt
    awk 'FNR == 1 {print \$1}' ${base}_ref2_percent.txt > ${base}_ref2_percent_num_parse.txt
    ref_2_percent=\$(sed -n '1p' < ${base}_ref2_percent_num_parse.txt | xargs)
    mixed_inf_cov=40

    # Check if Ref #2 has percent genome coverage higher than 40%. If true, Map Ref2 as a mixed infection.
    if [ "\$ref_2_percent" > "\$mixed_inf_cov" ]; then

    echo 'MIXED INFECTION - Mapping TOP TWO Genome References.'

    sed -n '1p' < ${base}_most_mapped_ref.txt > ${base}_ref1_name.txt
    sed -n '2p' < ${base}_most_mapped_ref.txt > ${base}_ref2_name.txt

    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    id_ref_1=\$(awk 'FNR==1{print val,\$2}' ${base}_ref1_name.txt | xargs)
    id_ref_2=\$(awk 'FNR==1{print val,\$2}'  ${base}_ref2_name.txt | xargs)

    ref_coverage_ref1=\$(awk '{print \$1; exit}' ${base}_ref1_name.txt)
    ref_coverage_ref2=\$(awk '{print \$1; exit}' ${base}_ref2_name.txt)

    samtools faidx ${Reference_hpv} \$id_ref_1 > ${base}_mapped_ref1_genome.fa
    samtools faidx ${Reference_hpv} \$id_ref_2 > ${base}_mapped_ref2_genome.fa

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_ref1_map2.sam ref=${base}_mapped_ref1_genome.fa threads=${task.cpus} covstats=${base}_ref1_map2_bbmap_out.txt covhist=${base}_ref1_map2_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_ref1_map2_stats.txt 2>&1

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_ref2_map2.sam ref=${base}_mapped_ref2_genome.fa threads=${task.cpus} covstats=${base}_ref2_map2_bbmap_out.txt covhist=${base}_ref2_map2_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_ref2_map2_stats.txt 2>&1

    samtools view -S -b ${base}_ref1_map2.sam > ${base}_ref1_map2.bam
    samtools view -S -b ${base}_ref2_map2.sam > ${base}_ref2_map2.bam

    samtools sort -@ 4 ${base}_ref1_map2.bam > ${base}_ref1_map2.sorted.bam
    samtools sort -@ 4 ${base}_ref2_map2.bam > ${base}_ref2_map2.sorted.bam

    samtools index ${base}_ref1_map2.sorted.bam
    samtools index ${base}_ref2_map2.sorted.bam

    samtools idxstats ${base}_ref1_map2.sorted.bam > ${base}_ref1_map2_idxstats.txt
    samtools idxstats ${base}_ref2_map2.sorted.bam > ${base}_ref2_map2_idxstats.txt

    head -n 1 ${base}_mapped_ref1_genome.fa > ${base}_mapped_ref1_edited_genome.fa
    head -n 1 ${base}_mapped_ref2_genome.fa > ${base}_mapped_ref2_edited_genome.fa

    grep -v ">" ${base}_mapped_ref1_genome.fa | sed 's/U/T/g' >> ${base}_mapped_ref1_edited_genome.fa
    grep -v ">" ${base}_mapped_ref2_genome.fa | sed 's/U/T/g' >> ${base}_mapped_ref2_edited_genome.fa

    mv ${base}_mapped_ref1_edited_genome.fa ${base}_mapped_ref1_genome.fa
    mv ${base}_mapped_ref2_edited_genome.fa ${base}_mapped_ref2_genome.fa

    samtools faidx ${base}_mapped_ref1_genome.fa
    samtools faidx ${base}_mapped_ref2_genome.fa
    cp ${base}_mapped_ref1_genome.fa ${base}_mapped_ref_genome.fa
    samtools faidx ${base}_mapped_ref_genome.fa

    awk 'NR == 2 || \$5 > max {number = \$3; max = \$5} END {if (NR) print number, max}' < ${base}_ref1_map2_bbmap_out.txt >  ${base}_ref1_most_mapped_ref_size_out.txt
    awk 'NR == 2 || \$5 > max {number = \$3; max = \$5} END {if (NR) print number, max}' < ${base}_ref2_map2_bbmap_out.txt >  ${base}_ref2_most_mapped_ref_size_out.txt

    id_ref_size=\$(awk 'FNR==1{print val,\$1}' ${base}_ref1_most_mapped_ref_size_out.txt | xargs)
    id_ref2_size=\$(awk 'FNR==1{print val,\$1}' ${base}_ref2_most_mapped_ref_size_out.txt | xargs)    

    echo \$id_ref_size >> ${base}_most_mapped_ref_size.txt
    echo \$id_ref2_size >> ${base}_most_mapped_ref2_size.txt

    reads_mapped=\$(cat ${base}_ref1_map2_stats.txt | grep "mapped:" | cut -d\$'\\t' -f3)
    reads_mapped_ref2=\$(cat ${base}_ref2_map2_stats.txt | grep "mapped:" | cut -d\$'\\t' -f3)

    printf "\$reads_mapped" >> ${base}_num_mapped.txt

    cp ${base}_most_mapped_ref.txt ${base}_most_mapped_ref.txt
    cp ${base}_mapped_ref1_genome.fa ${base}_mapped_ref_genome.fa
    cp ${base}_ref1_map1.sam ${base}_map1.sam
    cp ${base}_ref1_map2.sam ${base}_map2.sam
    cp ${base}_ref1_map2.sorted.bam ${base}.sorted.bam
    cp ${base}_ref1_map2_idxstats.txt ${base}_idxstats.txt
    cp ${base}_ref1_most_mapped_ref_size_out.txt ${base}_most_mapped_ref_size_out.txt
    cp ${base}_hpv_all_map1_bbmap_out.txt ${base}_map1_bbmap_out.txt
    cp ${base}_ref1_map2_bbmap_out.txt ${base}_map2_bbmap_out.txt
    cp ${base}_hpv_all_map1_stats.txt ${base}_map1_stats.txt
    cp ${base}_ref1_map2_stats.txt ${base}_map2_stats.txt
    cp ${base}_hpv_all_map1_histogram.txt ${base}_map1_histogram.txt   
    cp ${base}_ref1_map2_histogram.txt ${base}_map2_histogram.txt

    mkdir ${params.outdir}/bam_map2_mixed_infection/
    mkdir ${params.outdir}/sam_map2_mixed_infection/
    mkdir ${params.outdir}/txt_bbmap_map1_stats_mixed_infection/
    mkdir ${params.outdir}/txt_bbmap_map2_stats_mixed_infection/
    mkdir ${params.outdir}/txt_indxstats_mapped_refs_mixed_infection/
    mkdir ${params.outdir}/ref_id_mixed_infection/
    mkdir ${params.outdir}/ref_fasta_mixed_infection/
    mkdir ${params.outdir}/ref_size_mixed_infection/
    mkdir ${params.outdir}/ref_fai_index_mixed_infection/
    mkdir ${params.outdir}/hpv_ref_all/
    mkdir ${params.outdir}/hpv_ref_14/

    mv ${base}_ref1_map2_bbmap_out.txt ${params.outdir}/txt_bbmap_map1_stats_mixed_infection/
    mv ${base}_most_mapped_ref2_size.txt ${params.outdir}/ref_size_mixed_infection/
    mv ${base}_ref2_map2.sam ${params.outdir}/sam_map2_mixed_infection/
    mv ${base}_mapped_ref2_genome.fa ${params.outdir}/ref_fasta_mixed_infection/
    mv ${base}_ref2_map2.sorted.bam ${params.outdir}/bam_map2_mixed_infection/
    mv ${base}_ref2_map2_idxstats.txt ${params.outdir}/txt_indxstats_mapped_refs_mixed_infection/
    mv ${base}_hpv_all_map1_stats.txt ${params.outdir}/hpv_ref_all/
    mv ${base}_hpv_14_map1_stats.txt ${params.outdir}/hpv_ref_14/
    mv ${base}_hpv_all_map2_stats.txt ${params.outdir}/hpv_ref_all/
    mv ${base}_hpv_14_map2_stats.txt ${params.outdir}/hpv_ref_14/
    mv ${base}_ref2_map2_bbmap_out.txt ${params.outdir}/txt_bbmap_map2_stats_mixed_infection/
    mv ${base}_ref2_map2_stats.txt ${params.outdir}/txt_bbmap_map2_stats_mixed_infection/
    mv ${base}_hpv_all_map1.sam ${params.outdir}/hpv_ref_all/
    mv ${base}_hpv_14_map1.sam ${params.outdir}/hpv_ref_14/
    cp ${base}_most_mapped_ref.txt ${params.outdir}/ref_id_mixed_infection/
    cp ${base}_mapped_ref2_genome.fai ${params.outdir}/ref_fai_index_mixed_infection/
    cp ${base}_ref2_percent_num_parse.txt ${params.outdir}/hpv_ref_all/

    # Summary Statistics
    cp ${base}_summary.csv ${base}_summary2.csv
    printf ",\$id_ref_1" >> ${base}_summary2.csv
    printf ",\$id_ref_2" >> ${base}_summary2.csv
    printf ",\$id_ref_size" >> ${base}_summary2.csv
    printf ",\$id_ref2_size" >> ${base}_summary2.csv
    printf ",\$reads_mapped" >> ${base}_summary2.csv
    printf ",\$reads_mapped_ref2" >> ${base}_summary2.csv
    printf ",\$ref_coverage_ref1" >> ${base}_summary2.csv
    printf ",\$ref_coverage_ref2" >> ${base}_summary2.csv

    else

    echo 'NOT A MIXED INFECTION - Mapping TOP Genome Reference.'

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${Reference_hpv} threads=${task.cpus} covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map1_stats.txt 2>&1
    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${Reference_hpv} \$id > ${base}_mapped_ref_genome.fa
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map2.sam ref=${base}_mapped_ref_genome.fa threads=${task.cpus} covstats=${base}_map2_bbmap_out.txt covhist=${base}_map2_histogram.txt local=true interleaved=false maxindel=20 strictmaxindel -Xmx6g > ${base}_map2_stats.txt 2>&1
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

    fi


    # Influenza B
    elif grep -q \$all_ref_id "${base}_inbflb_ids.txt";
    then
    echo "< Accession found in Influenza B multifasta file. hrv_ref_Influenza_b.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${Reference_inflb} threads=${task.cpus} covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map1_stats.txt 2>&1
    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${Reference_inflb} \$id > ${base}_mapped_ref_genome.fa
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map2.sam ref=${base}_mapped_ref_genome.fa threads=${task.cpus} covstats=${base}_map2_bbmap_out.txt covhist=${base}_map2_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map2_stats.txt 2>&1
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


    # Human Coronavirus
    elif grep -q \$all_ref_id "${base}_hcov_ids.txt";
    then
    echo "Accession found in HCoVs multifasta file. hrv_ref_hcov.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${Reference_hcov} threads=${task.cpus} covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false maxindel=20 strictmaxindel -Xmx6g > ${base}_map1_stats.txt 2>&1
    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${Reference_hcov} \$id > ${base}_mapped_ref_genome.fa
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map2.sam ref=${base}_mapped_ref_genome.fa threads=${task.cpus} covstats=${base}_map2_bbmap_out.txt covhist=${base}_map2_histogram.txt local=true interleaved=false maxindel=20 strictmaxindel -Xmx6g > ${base}_map2_stats.txt 2>&1
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

    # HPIV3 - Human parainfluenza virus 3
    elif grep -q \$all_ref_id "${base}_hpiv3.txt";
    then
    echo "Accession found in HPIV3 multifasta file. hrv_ref_hpiv3.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${Reference_hpiv3} threads=${task.cpus} covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map1_stats.txt 2>&1
    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${Reference_hpiv3} \$id > ${base}_mapped_ref_genome.fa
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map2.sam ref=${base}_mapped_ref_genome.fa threads=${task.cpus} covstats=${base}_map2_bbmap_out.txt covhist=${base}_map2_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map2_stats.txt 2>&1
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

    else

    # Not Rhinovirus, HPV, Inlfuenza B, OR Human Coronavirus - Use Regular Mapping Settings
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${base}_all_ref.fa threads=${task.cpus} covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false maxindel=20 strictmaxindel -Xmx6g > ${base}_map1_stats.txt 2>&1
    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${base}_all_ref.fa \$id > ${base}_mapped_ref_genome.fa
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map2.sam ref=${base}_mapped_ref_genome.fa threads=${task.cpus} covstats=${base}_map2_bbmap_out.txt covhist=${base}_map2_histogram.txt local=true interleaved=false maxindel=20 strictmaxindel -Xmx6g > ${base}_map2_stats.txt 2>&1
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

    fi

    """
}
} else {
process Mapping_PE {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"     
    errorStrategy 'retry'
    maxRetries 3

    input: 
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_summary.csv") from Trim_out_PE
        file Reference_rv
        file Reference_hcov
        file Reference_hpv
        file Reference_inflb
        file Reference_hpiv3

    output:
        tuple val(base), file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"), file("${base}_summary2.csv"), file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),env(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),env(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") into Everything_PE_ch
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

    cat ${Reference_hpv} ${Reference_rv} ${Reference_inflb} ${Reference_hcov} ${Reference_hpiv3} > ${base}_all_ref.fa

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_all_ref.sam ref=${base}_all_ref.fa threads=${task.cpus} covstats=${base}_all_ref_bbmap_out.txt covhist=${base}_all_ref_histogram.txt local=true interleaved=false -Xmx10g > ${base}_all_ref_stats.txt 2>&1

    samtools view -S -b ${base}_all_ref.sam > ${base}_all_ref.bam
    samtools sort -@ 4 ${base}_all_ref.bam > ${base}_all_ref.sorted.bam
    samtools index ${base}_all_ref.sorted.bam
    samtools idxstats ${base}_all_ref.sorted.bam > ${base}_all_ref_idxstats.txt
    
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_all_ref_bbmap_out.txt > ${base}_all_ref_id.txt
    all_ref_id=\$(awk '{print \$1}' ${base}_all_ref_id.txt)

    grep -B 0 ">" ${Reference_rv} | tr -d ">" > ${base}_rv_ids.txt
    grep -B 0 ">" ${Reference_hpv} | tr -d ">" > ${base}_hpv_ids.txt
    grep -B 0 ">" ${Reference_inflb} | tr -d ">" > ${base}_inbflb_ids.txt
    grep -B 0 ">" ${Reference_hcov} | tr -d ">" > ${base}_hcov_ids.txt
    grep -B 0 ">" ${Reference_hpiv3} | tr -d ">" > ${base}_hpiv3.txt

    # Rhinovirus
    if grep -q \$all_ref_id "${base}_rv_ids.txt";
    then
    echo "< Accession found in Rhinovirus multifasta file. hrv_ref_rhinovirus.fa will be used for mapping."


    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${Reference_rv} threads=${task.cpus} covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map1_stats.txt 2>&1
    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${Reference_rv} \$id > ${base}_mapped_ref_genome.fa
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map2.sam ref=${base}_mapped_ref_genome.fa threads=${task.cpus} covstats=${base}_map2_bbmap_out.txt covhist=${base}_map2_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map2_stats.txt 2>&1
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


    # HPV
    elif grep -q \$all_ref_id "${base}_hpv_ids.txt";
    then
    echo "< Accession found in respiratory virus multifasta file. hrv_ref_hpv.fa will be used for mapping."


    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${Reference_hpv} threads=${task.cpus} covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map1_stats.txt 2>&1
    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${Reference_hpv} \$id > ${base}_mapped_ref_genome.fa
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map2.sam ref=${base}_mapped_ref_genome.fa threads=${task.cpus} covstats=${base}_map2_bbmap_out.txt covhist=${base}_map2_histogram.txt local=true interleaved=false maxindel=20 strictmaxindel -Xmx6g > ${base}_map2_stats.txt 2>&1
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


    # Influenza B
    elif grep -q \$all_ref_id "${base}_inbflb_ids.txt";
    then
    echo "< Accession found in Influenza B multifasta file. hrv_ref_Influenza_b.fa will be used for mapping."


    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${Reference_inflb} threads=${task.cpus} covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map1_stats.txt 2>&1
    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${Reference_inflb} \$id > ${base}_mapped_ref_genome.fa
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map2.sam ref=${base}_mapped_ref_genome.fa threads=${task.cpus} covstats=${base}_map2_bbmap_out.txt covhist=${base}_map2_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map2_stats.txt 2>&1
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


    # Human Coronavirus
    elif grep -q \$all_ref_id "${base}_hcov_ids.txt";
    then
    echo "Accession found in HCoVs multifasta file. hrv_ref_hcov.fa will be used for mapping."


    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${Reference_hcov} threads=${task.cpus} covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false maxindel=20 strictmaxindel -Xmx6g > ${base}_map1_stats.txt 2>&1
    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${Reference_hcov} \$id > ${base}_mapped_ref_genome.fa
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map2.sam ref=${base}_mapped_ref_genome.fa threads=${task.cpus} covstats=${base}_map2_bbmap_out.txt covhist=${base}_map2_histogram.txt local=true interleaved=false maxindel=20 strictmaxindel -Xmx6g > ${base}_map2_stats.txt 2>&1
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

    # HPIV3 - Human parainfluenza virus 3
    elif grep -q \$all_ref_id "${base}_hpiv3.txt";
    then
    echo "Accession found in HPIV3 multifasta file. hrv_ref_hpiv3.fa will be used for mapping."


    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${Reference_hpiv3} threads=${task.cpus} covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map1_stats.txt 2>&1
    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${Reference_hpiv3} \$id > ${base}_mapped_ref_genome.fa
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map2.sam ref=${base}_mapped_ref_genome.fa threads=${task.cpus} covstats=${base}_map2_bbmap_out.txt covhist=${base}_map2_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map2_stats.txt 2>&1
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

    else

    # Not Rhinovirus, Inlfuenza B, OR Human Coronavirus - Use Regular Mapping Settings
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${base}_all_ref.fa threads=${task.cpus} covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false maxindel=20 strictmaxindel -Xmx6g > ${base}_map1_stats.txt 2>&1
    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${base}_all_ref.fa \$id > ${base}_mapped_ref_genome.fa
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map2.sam ref=${base}_mapped_ref_genome.fa threads=${task.cpus} covstats=${base}_map2_bbmap_out.txt covhist=${base}_map2_histogram.txt local=true interleaved=false maxindel=20 strictmaxindel -Xmx6g > ${base}_map2_stats.txt 2>&1
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


    fi

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
    container "docker.io/paulrkcruz/hrv-pipeline:latest"    
	errorStrategy 'retry'
    maxRetries 3

    input: 
    tuple val(base), file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"), file("${base}_summary2.csv"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") from Everything_ch
    output:
    tuple val(base), file("${base}.bam") into Aligned_bam_ch, Bam_ch
    tuple val(base), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),env(bamsize),file("${base}.sorted.bam.bai"),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"), file("${base}_summary.csv"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") into Consensus_ch

    publishDir "${params.outdir}bam_map2", mode: 'copy', pattern:'*.sorted.bam*'  
    publishDir "${params.outdir}txt_bam_flagstats-map2", mode: 'copy', pattern:'*_flagstats.txt*'  

    script:
    """
    #!/bin/bash
    /usr/local/miniconda/bin/samtools view -S -b ${base}_map2.sam > ${base}.bam
    /usr/local/miniconda/bin/samtools sort -@ ${task.cpus} ${base}.bam > ${base}.sorted.bam
    /usr/local/miniconda/bin/samtools index ${base}.sorted.bam
    /usr/local/miniconda/bin/samtools flagstat ${base}.sorted.bam > ${base}_flagstats.txt
    /usr/local/miniconda/bin/bedtools genomecov -d -ibam ${base}.sorted.bam > ${base}_coverage.txt
    
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
    container "docker.io/paulrkcruz/hrv-pipeline:latest"    
	errorStrategy 'retry'
    maxRetries 3

    input: 
    tuple val(base), file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"), file("${base}_summary2.csv"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") from Everything_PE_ch
    output:
    tuple val(base), file("${base}.bam") into Aligned_bam_PE_ch, Bam_PE_ch
    tuple val(base), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),env(bamsize),file("${base}.sorted.bam.bai"),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"), file("${base}_summary.csv"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") into Consensus_PE_ch

    publishDir "${params.outdir}bam_map2", mode: 'copy', pattern:'*.sorted.bam*'  
    publishDir "${params.outdir}txt_bam_flagstats-map2", mode: 'copy', pattern:'*_flagstats.txt*' 

    script:
    """
    #!/bin/bash
    /usr/local/miniconda/bin/samtools view -S -b ${base}_map2.sam > ${base}.bam
    /usr/local/miniconda/bin/samtools sort -@ ${task.cpus} ${base}.bam > ${base}.sorted.bam
    /usr/local/miniconda/bin/samtools index ${base}.sorted.bam
    /usr/local/miniconda/bin/samtools flagstat ${base}.sorted.bam > ${base}_flagstats.txt
    /usr/local/miniconda/bin/bedtools genomecov -d -ibam ${base}.sorted.bam > ${base}_coverage.txt
    
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

if (params.singleEnd) {
process Generate_Consensus {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"     
	errorStrategy 'retry'
    // maxRetries 3
    input:
    tuple val(base), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),val(bamsize),file("${base}.sorted.bam.bai"),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"), file("${base}_summary.csv"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") from Consensus_ch
    
    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}_final_summary.csv"), val(bamsize), val(id),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt"),file("${base}.consensus_final.fa") into Consensus_Fasta_ch

    publishDir "${params.outdir}consensus-final", mode: 'copy', pattern:'*.consensus_final.fa*' 
    publishDir "${params.outdir}consensus_mpileup", mode: 'copy', pattern:'*.mpileup*' 
    
    script:
    """
    #!/bin/bash

    all_ref_id=\$(awk '{print \$1}' ${base}_all_ref_id.txt)

    # Rhinovirus
    if grep -q \$all_ref_id "${base}_rv_ids.txt"; 
    then
    echo "< Accession found in Rhinovirus multifasta file. hrv_ref_rhinovirus.fa will be used for mapping."

    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}.consensus.fa \\
        --min-BQ 15 \\
        --output ${base}.mpileup \\
        ${base}.sorted.bam
    cat ${base}.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}.consensus_final
    bedtools genomecov \\
        -bga \\
        -ibam ${base}.sorted.bam \\
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
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}.consensus_final.fa > ${base}.consensusfinal.fa

    sed 's/>.*/>${base}/' ${base}.consensusfinal.fa > ${base}.consensusfinal-renamed-header.fa
    grep -v "^>" ${base}.consensusfinal-renamed-header.fa | tr -cd N | wc -c > N.txt
    cp ${base}.consensusfinal-renamed-header.fa ${base}.consensus_final.fa
    
    num_ns=\$(awk 'FNR==1{print val,\$1}' N.txt)
    echo "\$num_ns/\$num_bases*100" | bc -l > n_percent.txt
    percent_n=\$(awk 'FNR==1{print val,\$1}' n_percent.txt)
    printf ",\$num_bases" >> ${base}_summary.csv
    printf ",\$percent_n" >> ${base}_summary.csv
    cp ${base}_summary.csv ${base}_final_summary.csv


    # HPV
    elif grep -q \$all_ref_id "${base}_hpv_ids.txt";
    then
    echo "< Accession found in HPV multifasta file. hrv_ref_hpv.fa will be used for mapping."


    cp ${base}_ref2_percent_num_parse.txt ${params.outdir}/hpv_ref_all/
    ref_2_percent=\$(sed -n '1p' < ${base}_ref2_percent_num_parse.txt | xargs)
    mixed_inf_cov=40

    # Check if Ref #2 has percent genome coverage higher than 40%. If true, Map Ref2 as a mixed infection.
    if [ "\$ref_2_percent" > "\$mixed_inf_cov" ]; then

    echo 'MIXED INFECTION - Mapping 2 Genomes'
    
    else

    echo 'NOT A MIXED INFECTION - Mapping 1 Genome'

    fi


cp ${params.outdir}/consensus_final_mixed_infection/${base}.consensus_final_mi.fa ${base}.consensus_final_ref2.fa

    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}.consensus.fa \\
        --min-BQ 15 \\
        --output ${base}.mpileup \\
        ${base}.sorted.bam
    cat ${base}.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}.consensus_final
    bedtools genomecov \\
        -bga \\
        -ibam ${base}.sorted.bam \\
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
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}.consensus_final.fa > ${base}.consensusfinal.fa

    sed 's/>.*/>${base}/' ${base}.consensusfinal.fa > ${base}.consensusfinal-renamed-header.fa
    grep -v "^>" ${base}.consensusfinal-renamed-header.fa | tr -cd N | wc -c > N.txt
    cp ${base}.consensusfinal-renamed-header.fa ${base}.consensus_final.fa
    
    num_ns=\$(awk 'FNR==1{print val,\$1}' N.txt)
    echo "\$num_ns/\$num_bases*100" | bc -l > n_percent.txt
    percent_n=\$(awk 'FNR==1{print val,\$1}' n_percent.txt)
    printf ",\$num_bases" >> ${base}_summary.csv
    printf ",\$percent_n" >> ${base}_summary.csv
    cp ${base}_summary.csv ${base}_final_summary.csv


    # Influenza B
    elif grep -q \$all_ref_id "${base}_inbflb_ids.txt";
    then
    echo "< Accession found in Influenza B multifasta file. hrv_ref_Influenza_b.fa will be used for mapping."

    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}.consensus.fa \\
        --min-BQ 15 \\
        --output ${base}.mpileup \\
        ${base}.sorted.bam
    cat ${base}.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}.consensus_final
    bedtools genomecov \\
        -bga \\
        -ibam ${base}.sorted.bam \\
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
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}.consensus_final.fa > ${base}.consensusfinal.fa

    sed 's/>.*/>${base}/' ${base}.consensusfinal.fa > ${base}.consensusfinal-renamed-header.fa
    grep -v "^>" ${base}.consensusfinal-renamed-header.fa | tr -cd N | wc -c > N.txt
    cp ${base}.consensusfinal-renamed-header.fa ${base}.consensus_final.fa
    
    num_ns=\$(awk 'FNR==1{print val,\$1}' N.txt)
    echo "\$num_ns/\$num_bases*100" | bc -l > n_percent.txt
    percent_n=\$(awk 'FNR==1{print val,\$1}' n_percent.txt)
    printf ",\$num_bases" >> ${base}_summary.csv
    printf ",\$percent_n" >> ${base}_summary.csv
    cp ${base}_summary.csv ${base}_final_summary.csv


    # Human Coronavirus
    elif grep -q \$all_ref_id "${base}_hcov_ids.txt";
    then
    echo "Accession found in HCoVs multifasta file. hrv_ref_hcov.fa will be used for mapping."

    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}.consensus.fa \\
        --min-BQ 15 \\
        --output ${base}.mpileup \\
        ${base}.sorted.bam
    cat ${base}.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}.consensus_final
    bedtools genomecov \\
        -bga \\
        -ibam ${base}.sorted.bam \\
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
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}.consensus_final.fa > ${base}.consensusfinal.fa

    sed 's/>.*/>${base}/' ${base}.consensusfinal.fa > ${base}.consensusfinal-renamed-header.fa
    grep -v "^>" ${base}.consensusfinal-renamed-header.fa | tr -cd N | wc -c > N.txt
    cp ${base}.consensusfinal-renamed-header.fa ${base}.consensus_final.fa
    
    num_ns=\$(awk 'FNR==1{print val,\$1}' N.txt)
    echo "\$num_ns/\$num_bases*100" | bc -l > n_percent.txt
    percent_n=\$(awk 'FNR==1{print val,\$1}' n_percent.txt)
    printf ",\$num_bases" >> ${base}_summary.csv
    printf ",\$percent_n" >> ${base}_summary.csv
    cp ${base}_summary.csv ${base}_final_summary.csv


    # HPIV3 - Human parainfluenza virus 3
    elif grep -q \$all_ref_id "${base}_hpiv3.txt";
    then
    echo "Accession found in HPIV3 multifasta file. hrv_ref_hpiv3.fa will be used for mapping."

    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}.consensus.fa \\
        --min-BQ 15 \\
        --output ${base}.mpileup \\
        ${base}.sorted.bam
    cat ${base}.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}.consensus_final
    bedtools genomecov \\
        -bga \\
        -ibam ${base}.sorted.bam \\
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
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}.consensus_final.fa > ${base}.consensusfinal.fa

    sed 's/>.*/>${base}/' ${base}.consensusfinal.fa > ${base}.consensusfinal-renamed-header.fa
    grep -v "^>" ${base}.consensusfinal-renamed-header.fa | tr -cd N | wc -c > N.txt
    cp ${base}.consensusfinal-renamed-header.fa ${base}.consensus_final.fa
    
    num_ns=\$(awk 'FNR==1{print val,\$1}' N.txt)
    echo "\$num_ns/\$num_bases*100" | bc -l > n_percent.txt
    percent_n=\$(awk 'FNR==1{print val,\$1}' n_percent.txt)
    printf ",\$num_bases" >> ${base}_summary.csv
    printf ",\$percent_n" >> ${base}_summary.csv
    cp ${base}_summary.csv ${base}_final_summary.csv

    else

    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}.consensus.fa \\
        --min-BQ 15 \\
        --output ${base}.mpileup \\
        ${base}.sorted.bam
    cat ${base}.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}.consensus_final
    bedtools genomecov \\
        -bga \\
        -ibam ${base}.sorted.bam \\
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
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}.consensus_final.fa > ${base}.consensusfinal.fa

    sed 's/>.*/>${base}/' ${base}.consensusfinal.fa > ${base}.consensusfinal-renamed-header.fa
    grep -v "^>" ${base}.consensusfinal-renamed-header.fa | tr -cd N | wc -c > N.txt
    cp ${base}.consensusfinal-renamed-header.fa ${base}.consensus_final.fa
    
    num_ns=\$(awk 'FNR==1{print val,\$1}' N.txt)
    echo "\$num_ns/\$num_bases*100" | bc -l > n_percent.txt
    percent_n=\$(awk 'FNR==1{print val,\$1}' n_percent.txt)
    printf ",\$num_bases" >> ${base}_summary.csv
    printf ",\$percent_n" >> ${base}_summary.csv
    cp ${base}_summary.csv ${base}_final_summary.csv

    fi

    """

}
} else {
process Generate_Consensus_PE {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"     
	errorStrategy 'retry'
    // maxRetries 3
    input:
    tuple val(base), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),val(bamsize),file("${base}.sorted.bam.bai"),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"), file("${base}_summary.csv"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") from Consensus_ch_PE
    
    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}_final_summary.csv"), val(bamsize), val(id),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt"),file("${base}.consensus_final.fa") into Consensus_Fasta_ch_PE

    publishDir "${params.outdir}consensus-final", mode: 'copy', pattern:'*.consensus_final.fa*' 
    publishDir "${params.outdir}consensus_mpileup", mode: 'copy', pattern:'*.mpileup*'
    
    script:
    """
    #!/bin/bash

    all_ref_id=\$(awk '{print \$1}' ${base}_all_ref_id.txt)

    # Rhinovirus
    if grep -q \$all_ref_id "${base}_rv_ids.txt"; 
    then
    echo "< Accession found in Rhinovirus multifasta file. hrv_ref_rhinovirus.fa will be used for mapping."

    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}.consensus.fa \\
        --min-BQ 15 \\
        --output ${base}.mpileup \\
        ${base}.sorted.bam
    cat ${base}.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}.consensus_final
    bedtools genomecov \\
        -bga \\
        -ibam ${base}.sorted.bam \\
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
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}.consensus_final.fa > ${base}.consensusfinal.fa

    sed 's/>.*/>${base}/' ${base}.consensusfinal.fa > ${base}.consensusfinal-renamed-header.fa
    grep -v "^>" ${base}.consensusfinal-renamed-header.fa | tr -cd N | wc -c > N.txt
    cp ${base}.consensusfinal-renamed-header.fa ${base}.consensus_final.fa
    
    num_ns=\$(awk 'FNR==1{print val,\$1}' N.txt)
    echo "\$num_ns/\$num_bases*100" | bc -l > n_percent.txt
    percent_n=\$(awk 'FNR==1{print val,\$1}' n_percent.txt)
    printf ",\$num_bases" >> ${base}_summary.csv
    printf ",\$percent_n" >> ${base}_summary.csv
    cp ${base}_summary.csv ${base}_final_summary.csv


    # HPV
    elif grep -q \$all_ref_id "${base}_hpv_ids.txt";
    then
    echo "< Accession found in HPV multifasta file. hrv_ref_hpv.fa will be used for mapping."

    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}.consensus.fa \\
        --min-BQ 15 \\
        --output ${base}.mpileup \\
        ${base}.sorted.bam
    cat ${base}.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}.consensus_final
    bedtools genomecov \\
        -bga \\
        -ibam ${base}.sorted.bam \\
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
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}.consensus_final.fa > ${base}.consensusfinal.fa

    sed 's/>.*/>${base}/' ${base}.consensusfinal.fa > ${base}.consensusfinal-renamed-header.fa
    grep -v "^>" ${base}.consensusfinal-renamed-header.fa | tr -cd N | wc -c > N.txt
    cp ${base}.consensusfinal-renamed-header.fa ${base}.consensus_final.fa
    
    num_ns=\$(awk 'FNR==1{print val,\$1}' N.txt)
    echo "\$num_ns/\$num_bases*100" | bc -l > n_percent.txt
    percent_n=\$(awk 'FNR==1{print val,\$1}' n_percent.txt)
    printf ",\$num_bases" >> ${base}_summary.csv
    printf ",\$percent_n" >> ${base}_summary.csv
    cp ${base}_summary.csv ${base}_final_summary.csv


    # Influenza B
    elif grep -q \$all_ref_id "${base}_inbflb_ids.txt";
    then
    echo "< Accession found in Influenza B multifasta file. hrv_ref_Influenza_b.fa will be used for mapping."

    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}.consensus.fa \\
        --min-BQ 15 \\
        --output ${base}.mpileup \\
        ${base}.sorted.bam
    cat ${base}.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}.consensus_final
    bedtools genomecov \\
        -bga \\
        -ibam ${base}.sorted.bam \\
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
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}.consensus_final.fa > ${base}.consensusfinal.fa

    sed 's/>.*/>${base}/' ${base}.consensusfinal.fa > ${base}.consensusfinal-renamed-header.fa
    grep -v "^>" ${base}.consensusfinal-renamed-header.fa | tr -cd N | wc -c > N.txt
    cp ${base}.consensusfinal-renamed-header.fa ${base}.consensus_final.fa
    
    num_ns=\$(awk 'FNR==1{print val,\$1}' N.txt)
    echo "\$num_ns/\$num_bases*100" | bc -l > n_percent.txt
    percent_n=\$(awk 'FNR==1{print val,\$1}' n_percent.txt)
    printf ",\$num_bases" >> ${base}_summary.csv
    printf ",\$percent_n" >> ${base}_summary.csv
    cp ${base}_summary.csv ${base}_final_summary.csv


    # Human Coronavirus
    elif grep -q \$all_ref_id "${base}_hcov_ids.txt";
    then
    echo "Accession found in HCoVs multifasta file. hrv_ref_hcov.fa will be used for mapping."

    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}.consensus.fa \\
        --min-BQ 15 \\
        --output ${base}.mpileup \\
        ${base}.sorted.bam
    cat ${base}.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}.consensus_final
    bedtools genomecov \\
        -bga \\
        -ibam ${base}.sorted.bam \\
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
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}.consensus_final.fa > ${base}.consensusfinal.fa

    sed 's/>.*/>${base}/' ${base}.consensusfinal.fa > ${base}.consensusfinal-renamed-header.fa
    grep -v "^>" ${base}.consensusfinal-renamed-header.fa | tr -cd N | wc -c > N.txt
    cp ${base}.consensusfinal-renamed-header.fa ${base}.consensus_final.fa
    
    num_ns=\$(awk 'FNR==1{print val,\$1}' N.txt)
    echo "\$num_ns/\$num_bases*100" | bc -l > n_percent.txt
    percent_n=\$(awk 'FNR==1{print val,\$1}' n_percent.txt)
    printf ",\$num_bases" >> ${base}_summary.csv
    printf ",\$percent_n" >> ${base}_summary.csv
    cp ${base}_summary.csv ${base}_final_summary.csv


    # HPIV3 - Human parainfluenza virus 3
    elif grep -q \$all_ref_id "${base}_hpiv3.txt";
    then
    echo "Accession found in HPIV3 multifasta file. hrv_ref_hpiv3.fa will be used for mapping."

    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}.consensus.fa \\
        --min-BQ 15 \\
        --output ${base}.mpileup \\
        ${base}.sorted.bam
    cat ${base}.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}.consensus_final
    bedtools genomecov \\
        -bga \\
        -ibam ${base}.sorted.bam \\
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
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}.consensus_final.fa > ${base}.consensusfinal.fa

    sed 's/>.*/>${base}/' ${base}.consensusfinal.fa > ${base}.consensusfinal-renamed-header.fa
    grep -v "^>" ${base}.consensusfinal-renamed-header.fa | tr -cd N | wc -c > N.txt
    cp ${base}.consensusfinal-renamed-header.fa ${base}.consensus_final.fa
    
    num_ns=\$(awk 'FNR==1{print val,\$1}' N.txt)
    echo "\$num_ns/\$num_bases*100" | bc -l > n_percent.txt
    percent_n=\$(awk 'FNR==1{print val,\$1}' n_percent.txt)
    printf ",\$num_bases" >> ${base}_summary.csv
    printf ",\$percent_n" >> ${base}_summary.csv
    cp ${base}_summary.csv ${base}_final_summary.csv


    else

    samtools mpileup \\
        --count-orphans \\
        --no-BAQ \\
        --max-depth 50000 \\
        --fasta-ref ${base}.consensus.fa \\
        --min-BQ 15 \\
        --output ${base}.mpileup \\
        ${base}.sorted.bam
    cat ${base}.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}.consensus_final
    bedtools genomecov \\
        -bga \\
        -ibam ${base}.sorted.bam \\
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
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}.consensus_final.fa > ${base}.consensusfinal.fa

    sed 's/>.*/>${base}/' ${base}.consensusfinal.fa > ${base}.consensusfinal-renamed-header.fa
    grep -v "^>" ${base}.consensusfinal-renamed-header.fa | tr -cd N | wc -c > N.txt
    cp ${base}.consensusfinal-renamed-header.fa ${base}.consensus_final.fa
    
    num_ns=\$(awk 'FNR==1{print val,\$1}' N.txt)
    echo "\$num_ns/\$num_bases*100" | bc -l > n_percent.txt
    percent_n=\$(awk 'FNR==1{print val,\$1}' n_percent.txt)
    printf ",\$num_bases" >> ${base}_summary.csv
    printf ",\$percent_n" >> ${base}_summary.csv
    cp ${base}_summary.csv ${base}_final_summary.csv

    fi

    """
}
}

if (params.singleEnd) {
process Final_Mapping {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"     
	errorStrategy 'retry'
    // maxRetries 3

    input:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}_final_summary.csv"), val(bamsize), val(id),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt"),file("${base}.consensus_final.fa") from Consensus_Fasta_ch

    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}_summary.csv"), val(bamsize), val(id),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt"),file("${base}.consensus_final.fa"),file("${base}_map3.sam") from into Mapping_Final_ch

    publishDir "${params.outdir}mpileup_map3", mode: 'copy', pattern:'*.mpileup*'
    publishDir "${params.outdir}bam_map3", mode: 'copy', pattern:'*.mpileup*'

    script:

    """
    #!/bin/bash
    
    # FINAL MAPPING - Map3

    # Rhinovirus
    if grep -q \$all_ref_id "${base}_rv_ids.txt"; 
    then
    echo "< Accession found in Rhinovirus multifasta file. hrv_ref_rhinovirus.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats_map3.txt 2>&1

    samtools view -S -b ${base}_map3.sam > ${base}_map3.bam
    samtools sort -@ 4 ${base}_map3.bam > ${base}_map3.sorted.bam
    samtools index ${base}_map3.sorted.bam
    cp ${base}_final_summary.csv ${base}_summary.csv


    # HPV
    elif grep -q \$all_ref_id "${base}_hpv_ids.txt";
    then
    echo "< Accession found in HPV multifasta file. hrv_ref_hpv.fa will be used for mapping."

    cp ${params.outdir}/hpv_ref_all/${base}_ref2_percent_num_parse.txt ${base}_ref2_percent_num_parse.txt
    ref_2_percent=\$(sed -n '1p' < ${base}_ref2_percent_num_parse.txt | xargs)
    mixed_inf_cov=40

    # Check if Ref #2 has percent genome coverage higher than 40%. If true, Map Ref2 as a mixed infection.
    if [ "\$ref_2_percent" > "\$mixed_inf_cov" ]; then

    echo 'MIXED INFECTION - Mapping 2 Genomes'

    # REF 1
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats_map3.txt 2>&1
    
    cp ${params.outdir}/consensus_final_mixed_infection/${base}.consensus_final_mi.fa ${base}.consensus_final_ref2.fa
    # REF 2
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus_final_ref2.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats_map3.txt 2>&1




    samtools view -S -b ${base}_map3.sam > ${base}_map3.bam
    samtools sort -@ 4 ${base}_map3.bam > ${base}.map3.sorted.bam
    samtools index ${base}.map3.sorted.bam

    picard MarkDuplicates -I ${base}.map3.sorted.bam -O ${base}.map3.sorted.deduplicated.bam -M ${base}.map3.sorted.metrics.txt -REMOVE_DUPLICATES  TRUE -ASSUME_SORTED TRUE -VALIDATION_STRINGENCY  SILENT

    samtools view -F 0x40 ${base}.map3.sorted.deduplicated.bam | cut -f1 | sort | uniq | wc -l > ${base}.map3.sorted.deduplicated.txt
    samtools view -F 0x40 ${base}.map3.sorted.bam | cut -f1 | sort | uniq | wc -l > ${base}.map3.sorted.txt

    deduplicated_reads=\$(head -n 1 ${base}.map3.sorted.deduplicated.txt)
    non_deduplicated_reads=\$(head -n 1 ${base}.map3.sorted.txt)

    num_trimmed=\$(cat ${base}_num_trimmed.txt | tr -d " \t\n\r" | sed -n '1 p')
    echo "\$deduplicated_reads/\$num_trimmed*100" | bc -l > reads-on-t_percent_dedup.txt
    reads_on_target_dedup=\$(awk 'FNR==1{print val,\$1}' reads-on-t_percent_dedup.txt)
    
    mkdir ${params.outdir}/map3_bam_deduplicated/
    mv ${base}.map3.sorted.deduplicated.bam ${params.outdir}/map3_bam_deduplicated/

    num_ns=\$(awk 'FNR==1{print val,\$1}' N.txt)
    echo "\$num_ns/\$num_bases*100" | bc -l > n_percent.txt
    percent_n=\$(awk 'FNR==1{print val,\$1}' n_percent.txt)
    printf ",\$num_bases" >> ${base}_final_summary.csv
    printf ",\$percent_n" >> ${base}_final_summary.csv
    printf ",\$non_deduplicated_reads" >> ${base}_final_summary.csv
    printf ",\$deduplicated_reads" >> ${base}_final_summary.csv
    printf ",\$reads_on_target_dedup" >> ${base}_final_summary.csv
    cp ${base}_final_summary.csv ${base}_summary.csv





    else

    echo 'NOT A MIXED INFECTION - Mapping 1 Genome'

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats_map3.txt 2>&1

    samtools view -S -b ${base}_map3.sam > ${base}_map3.bam
    samtools sort -@ 4 ${base}_map3.bam > ${base}.map3.sorted.bam
    samtools index ${base}.map3.sorted.bam

    picard MarkDuplicates -I ${base}.map3.sorted.bam -O ${base}.map3.sorted.deduplicated.bam -M ${base}.map3.sorted.metrics.txt -REMOVE_DUPLICATES  TRUE -ASSUME_SORTED TRUE -VALIDATION_STRINGENCY  SILENT

    samtools view -F 0x40 ${base}.map3.sorted.deduplicated.bam | cut -f1 | sort | uniq | wc -l > ${base}.map3.sorted.deduplicated.txt
    samtools view -F 0x40 ${base}.map3.sorted.bam | cut -f1 | sort | uniq | wc -l > ${base}.map3.sorted.txt

    deduplicated_reads=\$(head -n 1 ${base}.map3.sorted.deduplicated.txt)
    non_deduplicated_reads=\$(head -n 1 ${base}.map3.sorted.txt)

    num_trimmed=\$(cat ${base}_num_trimmed.txt | tr -d " \t\n\r" | sed -n '1 p')
    echo "\$deduplicated_reads/\$num_trimmed*100" | bc -l > reads-on-t_percent_dedup.txt
    reads_on_target_dedup=\$(awk 'FNR==1{print val,\$1}' reads-on-t_percent_dedup.txt)
    
    mkdir ${params.outdir}/map3_bam_deduplicated/
    mv ${base}.map3.sorted.deduplicated.bam ${params.outdir}/map3_bam_deduplicated/

    num_ns=\$(awk 'FNR==1{print val,\$1}' N.txt)
    echo "\$num_ns/\$num_bases*100" | bc -l > n_percent.txt
    percent_n=\$(awk 'FNR==1{print val,\$1}' n_percent.txt)
    printf ",\$num_bases" >> ${base}_final_summary.csv
    printf ",\$percent_n" >> ${base}_final_summary.csv
    printf ",\$non_deduplicated_reads" >> ${base}_final_summary.csv
    printf ",\$deduplicated_reads" >> ${base}_final_summary.csv
    printf ",\$reads_on_target_dedup" >> ${base}_final_summary.csv
    cp ${base}_final_summary.csv ${base}_summary.csv

    fi


    # Influenza B
    elif grep -q \$all_ref_id "${base}_inbflb_ids.txt";
    then
    echo "< Accession found in Influenza B multifasta file. hrv_ref_Influenza_b.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats_map3.txt 2>&1

    samtools view -S -b ${base}_map3.sam > ${base}_map3.bam
    samtools sort -@ 4 ${base}_map3.bam > ${base}_map3.sorted.bam
    samtools index ${base}_map3.sorted.bam
    cp ${base}_final_summary.csv ${base}_summary.csv


    # Human Coronavirus
    elif grep -q \$all_ref_id "${base}_hcov_ids.txt";
    then
    echo "Accession found in HCoVs multifasta file. hrv_ref_hcov.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats_map3.txt 2>&1

    samtools view -S -b ${base}_map3.sam > ${base}_map3.bam
    samtools sort -@ 4 ${base}_map3.bam > ${base}_map3.sorted.bam
    samtools index ${base}_map3.sorted.bam
    cp ${base}_final_summary.csv ${base}_summary.csv


    # HPIV3 - Human parainfluenza virus 3
    elif grep -q \$all_ref_id "${base}_hpiv3.txt";
    then
    echo "Accession found in HPIV3 multifasta file. hrv_ref_hpiv3.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats_map3.txt 2>&1

    samtools view -S -b ${base}_map3.sam > ${base}_map3.bam
    samtools sort -@ 4 ${base}_map3.bam > ${base}_map3.sorted.bam
    samtools index ${base}_map3.sorted.bam
    cp ${base}_final_summary.csv ${base}_summary.csv

    else

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats_map3.txt 2>&1

    samtools view -S -b ${base}_map3.sam > ${base}_map3.bam
    samtools sort -@ 4 ${base}_map3.bam > ${base}_map3.sorted.bam
    samtools index ${base}_map3.sorted.bam
    cp ${base}_final_summary.csv ${base}_summary.csv

    fi


    """  
}
} else {
process Final_Mapping_PE {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"     
	errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(base), file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus.fa"), file("${base}_final_summary.csv"), val(bamsize), val(id),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") from Consensus_Fasta_PE_ch

    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus_final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map4.sam"), file("${base}_map4.bam"), file("${base}.map4.sorted.bam"), file("${base}.map4.sorted.bam.bai"), file("${base}_final_mapping_stats.txt"), file("${base}_final_mapping_stats_map4.txt"), file("${base}.mpileup"), file("${base}_summary.csv"), file("${base}.trimmed.fastq.gz"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") into Mapping_Final_PE_ch

    publishDir "${params.outdir}mpileup_map3", mode: 'copy', pattern:'*.mpileup*'
    publishDir "${params.outdir}bam_map3", mode: 'copy', pattern:'*.map3.sorted.bam*'
    publishDir "${params.outdir}sam_map3", mode: 'copy', pattern:'*_map3.sam*'
    publishDir "${params.outdir}mpileup_map4", mode: 'copy', pattern:'*_map4.mpileup*'
    publishDir "${params.outdir}bam_map4", mode: 'copy', pattern:'*.map4.sorted.bam*'
    publishDir "${params.outdir}sam_map4", mode: 'copy', pattern:'*_map4.sam*'
    publishDir "${params.outdir}consensus-final", mode: 'copy', pattern:'*.consensus_final.fa*'
    // publishDir "${params.outdir}consensus-ivar-masked", mode: 'copy', pattern:'*.consensus.masked.fa*'
    publishDir "${params.outdir}txt_bbmap_final_mapping_stats_map3", mode: 'copy', pattern:'*_final_mapping_stats.txt*'
    publishDir "${params.outdir}txt_bbmap_final_mapping_stats_map4", mode: 'copy', pattern:'*_final_mapping_stats_map4.txt*'
    publishDir "${params.outdir}summary", mode: 'copy', pattern:'*_final_summary.csv*'

    script:

    """
    #!/bin/bash

    all_ref_id=\$(awk '{print \$1}' ${base}_all_ref_id.txt)

    # Rhinovirus
    if grep -q \$all_ref_id "${base}_rv_ids.txt"; 
    then
    echo "< Accession found in Rhinovirus multifasta file. hrv_ref_rhinovirus.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats.txt 2>&1

    # HPV
    elif grep -q \$all_ref_id "${base}_hpv_ids.txt";
    then
    echo "< Accession found in HPV multifasta file. hrv_ref_hpv.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats.txt 2>&1

    # Influenza B
    elif grep -q \$all_ref_id "${base}_inbflb_ids.txt";
    then
    echo "< Accession found in Influenza B multifasta file. hrv_ref_Influenza_b.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats.txt 2>&1

    # Human Coronavirus
    elif grep -q \$all_ref_id "${base}_hcov_ids.txt";
    then
    echo "Accession found in HCoVs multifasta file. hrv_ref_hcov.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus.fa threads=${task.cpus} local=true interleaved=false maxindel=20 -Xmx6g > ${base}_final_mapping_stats.txt 2>&1

    # HPIV3 - Human parainfluenza virus 3
    elif grep -q \$all_ref_id "${base}_hpiv3.txt";
    then
    echo "Accession found in HPIV3 multifasta file. hrv_ref_hpiv3.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus.fa threads=${task.cpus} local=true interleaved=false maxindel=20 -Xmx6g > ${base}_final_mapping_stats.txt 2>&1

    else

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats.txt 2>&1

    fi


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
    seqkit -is replace -p "^n+|n+\$" -r "" ${base}.consensus_final.fa > ${base}.consensusfinal.fa

    awk '/^>/{print ">${base}" ++i; next}{print}' < ${base}.consensusfinal.fa > ${base}.consensusfinal-renamed-header.fa
    grep -v "^>" ${base}.consensusfinal-renamed-header.fa | tr -cd N | wc -c > N.txt
    cp ${base}.consensusfinal-renamed-header.fa ${base}.consensus_final.fa

    num_ns=\$(awk 'FNR==1{print val,\$1}' N.txt)
    echo "\$num_ns/\$num_bases*100" | bc -l > n_percent.txt
    percent_n=\$(awk 'FNR==1{print val,\$1}' n_percent.txt)
    printf ",\$num_bases" >> ${base}_final_summary.csv
    printf ",\$percent_n" >> ${base}_final_summary.csv
    cp ${base}_final_summary.csv ${base}_summary.csv
    
    # Rhinovirus
    if grep -q \$all_ref_id "${base}_rv_ids.txt"; 
    then
    echo "< Accession found in Rhinovirus multifasta file. hrv_ref_rhinovirus.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map4.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats_map4.txt 2>&1

    # Respiratory Panel
    elif grep -q \$all_ref_id "${base}_hpv_ids.txt";
    then
    echo "< Accession found in respiratory virus multifasta file. hrv_ref_hpv.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map4.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=20 -Xmx6g > ${base}_final_mapping_stats_map4.txt 2>&1

    # Influenza B
    elif grep -q \$all_ref_id "${base}_inbflb_ids.txt";
    then
    echo "< Accession found in Influenza B multifasta file. hrv_ref_Influenza_b.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map4.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats_map4.txt 2>&1

    # Human Coronavirus
    elif grep -q \$all_ref_id "${base}_hcov_ids.txt";
    then
    echo "Accession found in HCoVs multifasta file. hrv_ref_hcov.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map4.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=20 -Xmx6g > ${base}_final_mapping_stats_map4.txt 2>&1

    # HPIV3 - Human parainfluenza virus 3
    elif grep -q \$all_ref_id "${base}_hpiv3.txt";
    then
    echo "Accession found in HPIV3 multifasta file. hrv_ref_hpiv3.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map4.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=20 -Xmx6g > ${base}_final_mapping_stats_map4.txt 2>&1

    else

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map4.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats_map4.txt 2>&1

    fi

    samtools view -S -b ${base}_map4.sam > ${base}_map4.bam
    samtools sort -@ 4 ${base}_map4.bam > ${base}.map4.sorted.bam
    samtools index ${base}.map4.sorted.bam

    """  
}
}
if (params.withMetadata) {
    if (params.singleEnd) {
process Summary_Generation {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"        
    errorStrategy 'retry'
    // maxRetries 3

    input:
    file SAMPLE_LIST from METADATA

    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus_final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map4.sam"), file("${base}_map4.bam"), file("${base}.map4.sorted.bam"), file("${base}.map4.sorted.bam.bai"), file("${base}_final_mapping_stats.txt"), file("${base}_final_mapping_stats_map4.txt"), file("${base}.mpileup"), file("${base}_summary.csv"), file("${base}.trimmed.fastq.gz"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") from Mapping_Final_ch  

    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus_final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map4.sam"), file("${base}_map4.bam"), file("${base}.map4.sorted.bam"), file("${base}.map4.sorted.bam.bai"), file("${base}_final_mapping_stats.txt"), file("${base}_final_mapping_stats_map4.txt"), file("${base}.mpileup"), file("${base}_final_summary.csv"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_sample_id.txt"), file("${base}_pcr_ct.txt"), file("${base}_method.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") into Final_Processing_final_ch  

    publishDir "${params.outdir}summary_withMetadata", mode: 'copy', pattern:'*_final_summary.csv*'

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
    printf ",\$Reads_On_Target" >> ${base}_summary.csv
    printf ",\$pcr_ct" >> ${base}_summary.csv
    printf ",\$method" >> ${base}_summary.csv
    printf ",\$NCBI_Name" >> ${base}_summary.csv 
    cp ${base}_summary.csv ${base}_final_summary.csv
    """
    }
    } else {
process Summary_Generation_PE {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"        
    errorStrategy 'retry'
    maxRetries 3

    input:
    file SAMPLE_LIST from SAMPLE_SHEET
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus_final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map4.sam"), file("${base}_map4.bam"), file("${base}.map4.sorted.bam"), file("${base}.map4.sorted.bam.bai"), file("${base}_final_mapping_stats.txt"), file("${base}_final_mapping_stats_map4.txt"), file("${base}.mpileup"), file("${base}_summary.csv"), file("${base}.trimmed.fastq.gz"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") from Mapping_Final_PE_ch    

    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus_final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map4.sam"), file("${base}_map4.bam"), file("${base}.map4.sorted.bam"), file("${base}.map4.sorted.bam.bai"), file("${base}_final_mapping_stats.txt"), file("${base}_final_mapping_stats_map4.txt"), file("${base}.mpileup"), file("${base}_final_summary.csv"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_sample_id.txt"), file("${base}_pcr_ct.txt"), file("${base}_method.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") into Final_Processing_PE_ch  

    publishDir "${params.outdir}summary_withMetadata", mode: 'copy', pattern:'*_final_summary.csv*'

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
    serotype="serotype"
    printf ",\$Reads_On_Target" >> ${base}_summary.csv
    printf ",\$pcr_ct" >> ${base}_summary.csv
    printf ",\$method" >> ${base}_summary.csv
    printf ",\$NCBI_Name" >> ${base}_summary.csv 
    cp ${base}_summary.csv ${base}_final_summary.csv
    """
    }
    }
    }
if (params.withSerotype) {
    if (params.singleEnd) {
process Serotyping {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"        
    errorStrategy 'retry'
    // maxRetries 3
    // errorStrategy 'ignore'

    input:
    file METADATA_INFO from METADATA
    file BLASTDB_VP1_1 from BLAST_DB_VP1_1
    file BLASTDB_VP1_2 from BLAST_DB_VP1_2
    file BLASTDB_VP1_3 from BLAST_DB_VP1_3
    file BLASTDB_VP1_4 from BLAST_DB_VP1_4
    file BLASTDB_VP1_5 from BLAST_DB_VP1_5
    file BLASTDB_VP1_6 from BLAST_DB_VP1_6
    file BLASTDB_VP1_7 from BLAST_DB_VP1_7
    file BLASTDB_VP1_8 from BLAST_DB_VP1_8

    file hpv_db_1 from BLAST_DB_ALL_1hpv
    file hpv_db_2 from BLAST_DB_ALL_2hpv
    file hpv_db_3 from BLAST_DB_ALL_3hpv
    file hpv_db_4 from BLAST_DB_ALL_4hpv
    file hpv_db_5 from BLAST_DB_ALL_5hpv
    file hpv_db_6 from BLAST_DB_ALL_6hpv
    file hpv_db_7 from BLAST_DB_ALL_7hpv
    file hpv_db_8 from BLAST_DB_ALL_8hpv
    file hpv_db_9 from BLAST_DB_ALL_9hpv
    file hpv_db_10 from BLAST_DB_ALL_10hpv

    file BLASTDB_ALL_1 from BLAST_DB_ALL_1
    file BLASTDB_ALL_2 from BLAST_DB_ALL_2
    file BLASTDB_ALL_3 from BLAST_DB_ALL_3
    file BLASTDB_ALL_4 from BLAST_DB_ALL_4
    file BLASTDB_ALL_5 from BLAST_DB_ALL_5
    file BLASTDB_ALL_6 from BLAST_DB_ALL_6
    file BLASTDB_ALL_7 from BLAST_DB_ALL_7
    file BLASTDB_ALL_8 from BLAST_DB_ALL_8
    file BLASTDB_ALL_9 from BLAST_DB_ALL_9
    file BLASTDB_ALL_10 from BLAST_DB_ALL_10

    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus_final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map4.sam"), file("${base}_map4.bam"), file("${base}.map4.sorted.bam"), file("${base}.map4.sorted.bam.bai"), file("${base}_final_mapping_stats.txt"), file("${base}_final_mapping_stats_map4.txt"), file("${base}.mpileup"), file("${base}_final_summary.csv"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_sample_id.txt"), file("${base}_pcr_ct.txt"), file("${base}_method.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") from Final_Processing_final_ch   

    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus_final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map4.sam"), file("${base}_map4.bam"), file("${base}.map4.sorted.bam"), file("${base}.map4.sorted.bam.bai"), file("${base}_final_mapping_stats.txt"), file("${base}_final_mapping_stats_map4.txt"), file("${base}.mpileup"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_sample_id.txt"), file("${base}_pcr_ct.txt"), file("${base}_method.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_collection_year.txt"), file("${base}_country_collected.txt"), file("${base}_blast_db_vp1.txt"), file("${base}_blast_db_all_ref.txt"), file("${base}_sample_stats.csv"), file("${base}_all_ref_id.txt"), file("${base}_nomen.txt") into Serotype_ch, Summary_file_ch
    tuple val(base), file("${base}_summary_final.csv") into Summary_cat_ch

    publishDir "${params.outdir}blast_serotype", mode: 'copy', pattern:'*_blast_db_vp1.txt*'
    publishDir "${params.outdir}blast_ref_genome", mode: 'copy', pattern:'*_blast_db_all_ref.txt*'    
    publishDir "${params.outdir}sample_Stats", mode: 'copy', pattern:'*_sample_stats.csv*'
    publishDir "${params.outdir}sample_Stats", mode: 'copy', pattern:'*.txt*'       
    publishDir "${params.outdir}sample_Stats", mode: 'copy', pattern:'*_collection_year.txt*'
    publishDir "${params.outdir}sample_Stats", mode: 'copy', pattern:'*_country_collected.txt*'
    publishDir "${params.outdir}sample_Stats", mode: 'copy', pattern:'*_nomen.txt*'
    publishDir "${params.outdir}summary_withserotype", mode: 'copy', pattern:'*_summary_final.csv*'

    script:

    """
    #!/bin/bash
    R1=${base}
    NCBI_Name=\${R1:4:6}
    SAMPLEName=\${R1:2:5}
    all_ref_id=\$(awk '{print \$1}' ${base}_all_ref_id.txt)

    # Rhinovirus
    if grep -q \$all_ref_id "${base}_rv_ids.txt";
    then
    echo "< Accession found in Rhinovirus multifasta file."

    csvgrep -c sample_id -r \$NCBI_Name ${METADATA_INFO} > ${base}_sample_stats.csv
    csvcut -c 1 ${base}_sample_stats.csv > ${base}_sample_id.txt
    csvcut -c 4 ${base}_sample_stats.csv > ${base}_collection_year.txt
    csvcut -c 5 ${base}_sample_stats.csv > ${base}_country_collected.txt
    csvcut -c 6 ${base}_sample_stats.csv > ${base}_biosample_name.txt
    csvcut -c 7 ${base}_sample_stats.csv > ${base}_biosample_accession.txt
    csvcut -c 8 ${base}_sample_stats.csv > ${base}_sra_accession.txt
    csvcut -c 10 ${base}_sample_stats.csv > ${base}_release_date.txt
    csvcut -c 9 ${base}_sample_stats.csv > ${base}_bioproject.txt

    sample_id=\$(cat ${base}_sample_id.txt | sed -n '2 p')
    collection_year=\$(cat ${base}_collection_year.txt | sed -n '2 p')
    country_collected=\$(cat ${base}_country_collected.txt | sed -n '2 p')

    blastn -out ${base}_blast_db_vp1.txt -query ${base}.consensus_final.fa -db ${BLASTDB_VP1_1} -outfmt 6 -task blastn -max_target_seqs 1 -evalue 1e-5

    serotype=\$(awk 'FNR==1{print val,\$2}' ${base}_blast_db_vp1.txt)
    cut -d "-" -f2- <<< "\$serotype" > ${base}_serotype-parse.txt
    serotype_parsed=\$(awk 'FNR==1{print val,\$1}' ${base}_serotype-parse.txt)
    rv='Rv'
    serotype_parse="\${rv} \${serotype_parsed}" 
    echo \$serotype_parse > ${base}_sero.txt
    cat ${base}_sero.txt | tr -d " \t\n\r" > ${base}_serot.txt
    serotype_parsed2=\$(awk 'FNR==1{print val,\$1}' ${base}_serot.txt)
    space='/'
    nomenclature="\${serotype_parsed2} \${space} \${country_collected} \${space} \${collection_year} \${space} \${NCBI_Name}" 
    echo \$nomenclature > ${base}_nomenclature.txt
    cat ${base}_nomenclature.txt | tr -d " \t\n\r" > ${base}_nomenclature_parsed.txt
    nomenclature_parsed=\$(awk 'FNR==1{print val,\$1}' ${base}_nomenclature_parsed.txt)	
    echo \$serotype_parsed2 | xargs > ${base}_serots.txt
    echo \$nomenclature_parsed | xargs > ${base}_nomen.txt
    serots=\$(awk 'FNR==1{print val,\$1}' ${base}_serots.txt)	
    nomen=\$(awk 'FNR==1{print val,\$1}' ${base}_nomen.txt)	

    blastn -out ${base}_blast_db_all_ref.txt -query ${base}.consensus_final.fa -db ${BLASTDB_ALL_1} -outfmt "5 std qlen" -task blastn -max_target_seqs 1 -evalue 1e-5

    awk 'NR==31' ${base}_blast_db_all_ref.txt > ${base}_strain.txt
    sed -i -e 's/<Hit_def>//g'  ${base}_strain.txt
    awk -F'</Hit_def>' '{print \$1}' ${base}_strain.txt | xargs > ${base}_strain-parsed.txt
    Reference_Name=\$(head -n 1 ${base}_strain-parsed.txt)
    biosample_name=\$(cat ${base}_biosample_name.txt | sed -n '2 p')
    biosample_accession=\$(cat ${base}_biosample_accession.txt | sed -n '2 p')
    sra_accession=\$(cat ${base}_sra_accession.txt | sed -n '2 p')
    release_date=\$(cat ${base}_release_date.txt | sed -n '2 p')
    bioproject=\$(cat ${base}_bioproject.txt | sed -n '2 p')
    
    printf ",\$serots" >> ${base}_final_summary.csv
    printf ",\$nomen" >> ${base}_final_summary.csv
    printf ",\$Reference_Name" >> ${base}_final_summary.csv
    printf ",\$biosample_name" >> ${base}_final_summary.csv
    printf ",\$biosample_accession" >> ${base}_final_summary.csv
    printf ",\$sra_accession" >> ${base}_final_summary.csv
    printf ",\$release_date" >> ${base}_final_summary.csv
    printf ",\$bioproject" >> ${base}_final_summary.csv
    cp ${base}_final_summary.csv ${base}_summary_final.csv
    
    
    # HPV
    elif grep -q \$all_ref_id "${base}_hpv_ids.txt";
    then
    echo "< Accession found in respiratory virus multifasta file."
    
    csvgrep -c sample_id -r \$SAMPLEName ${METADATA_INFO} > ${base}_sample_stats.csv
    csvcut -c 1 ${base}_sample_stats.csv > ${base}_sample_id.txt
    csvcut -c 4 ${base}_sample_stats.csv > ${base}_collection_year.txt
    csvcut -c 5 ${base}_sample_stats.csv > ${base}_country_collected.txt
    csvcut -c 6 ${base}_sample_stats.csv > ${base}_biosample_name.txt
    csvcut -c 7 ${base}_sample_stats.csv > ${base}_biosample_accession.txt
    csvcut -c 8 ${base}_sample_stats.csv > ${base}_sra_accession.txt
    csvcut -c 9 ${base}_sample_stats.csv > ${base}_release_date.txt
    csvcut -c 10 ${base}_sample_stats.csv > ${base}_bioproject.txt

    sample_id=\$(cat ${base}_sample_id.txt | sed -n '2 p')
    collection_year=\$(cat ${base}_collection_year.txt | sed -n '2 p')
    country_collected=\$(cat ${base}_country_collected.txt | sed -n '2 p')

    blastn -out ${base}_blast_db_vp1.txt -query ${base}.consensus_final.fa -db ${hpv_db_1} -outfmt "5 std qlen" -task blastn -max_target_seqs 1 -evalue 1e-5

    blastn -out ${base}_blast_db_all_ref.txt -query ${base}.consensus_final.fa -db ${BLASTDB_ALL_1} -outfmt "5 std qlen" -task blastn -max_target_seqs 1 -evalue 1e-5

    sed -n '31p' < ${base}_blast_db_vp1.txt > ${base}_genotype.txt
    sed -n '32p' < ${base}_blast_db_vp1.txt > ${base}_ref_id.txt

    ref_id_hpv=\$(awk -F'<Hit_accession>' '{print \$2}' ${base}_ref_id.txt > ${base}_ref_id_parsed_1.txt)
    ref_id_parse_1=\$(sed -n '1p' < ${base}_ref_id_parsed_1.txt) 
    echo "\$ref_id_parse_1" | cut -f1 -d"<" > ${base}_ref_id_parsed_2.txt
    ref_id_parsed=\$(sed -n '1p' < ${base}_ref_id_parsed_2.txt)

    genotype_hpv=\$(awk -F'<Hit_def>' '{print \$2}' ${base}_genotype.txt > ${base}_genotype_parsed_1.txt) 
    genotype_parse_1=\$(sed -n '1p' < ${base}_genotype_parsed_1.txt)
    echo "\$genotype_parse_1" | cut -f1 -d"<" > ${base}_genotype_parsed_2.txt
    genotype_parsed_2=\$(sed -n '1p' < ${base}_genotype_parsed_2.txt)

    cp ${base}_genotype_parsed_2.txt ${base}_serots.txt
    cp ${base}_ref_id_parsed_2.txt ${base}_nomen.txt

    biosample_name=\$(cat ${base}_biosample_name.txt | sed -n '2 p')
    biosample_accession=\$(cat ${base}_biosample_accession.txt | sed -n '2 p')
    sra_accession=\$(cat ${base}_sra_accession.txt | sed -n '2 p')
    release_date=\$(cat ${base}_release_date.txt | sed -n '2 p')
    bioproject=\$(cat ${base}_bioproject.txt | sed -n '2 p')

    printf ",\$ref_id_parsed" >> ${base}_final_summary.csv
    printf ",\$genotype_parsed_2" >> ${base}_final_summary.csv
    printf ",\$biosample_name" >> ${base}_final_summary.csv
    printf ",\$biosample_accession" >> ${base}_final_summary.csv
    printf ",\$sra_accession" >> ${base}_final_summary.csv
    printf ",\$release_date" >> ${base}_final_summary.csv
    printf ",\$bioproject" >> ${base}_final_summary.csv
    cp ${base}_final_summary.csv ${base}_summary_final.csv


    # Influenza B
    elif grep -q \$all_ref_id "${base}_inbflb_ids.txt";
    then
    echo "< Accession found in Influenza B multifasta file."


    printf ",\$serots" >> ${base}_final_summary.csv
    printf ",\$nomen" >> ${base}_final_summary.csv
    printf ",\$Reference_Name" >> ${base}_final_summary.csv
    printf ",\$biosample_name" >> ${base}_final_summary.csv
    printf ",\$biosample_accession" >> ${base}_final_summary.csv
    printf ",\$sra_accession" >> ${base}_final_summary.csv
    printf ",\$release_date" >> ${base}_final_summary.csv
    printf ",\$bioproject" >> ${base}_final_summary.csv
    cp ${base}_final_summary.csv ${base}_summary_final.csv

    # HPIV3 - Human parainfluenza virus 3
    elif grep -q \$all_ref_id "${base}_hpiv3.txt";
    then
    echo "Accession found in HPIV3 multifasta file."


    printf ",\$serots" >> ${base}_final_summary.csv
    printf ",\$nomen" >> ${base}_final_summary.csv
    printf ",\$Reference_Name" >> ${base}_final_summary.csv
    printf ",\$biosample_name" >> ${base}_final_summary.csv
    printf ",\$biosample_accession" >> ${base}_final_summary.csv
    printf ",\$sra_accession" >> ${base}_final_summary.csv
    printf ",\$release_date" >> ${base}_final_summary.csv
    printf ",\$bioproject" >> ${base}_final_summary.csv
    cp ${base}_final_summary.csv ${base}_summary_final.csv

    else

    printf ",\$serots" >> ${base}_final_summary.csv
    printf ",\$nomen" >> ${base}_final_summary.csv
    printf ",\$Reference_Name" >> ${base}_final_summary.csv
    printf ",\$biosample_name" >> ${base}_final_summary.csv
    printf ",\$biosample_accession" >> ${base}_final_summary.csv
    printf ",\$sra_accession" >> ${base}_final_summary.csv
    printf ",\$release_date" >> ${base}_final_summary.csv
    printf ",\$bioproject" >> ${base}_final_summary.csv
    cp ${base}_final_summary.csv ${base}_summary_final.csv

    fi



    """
    }
    } else {
process Serotyping_PE {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"        
    errorStrategy 'retry'
    maxRetries 3

    input:
    file SAMPLE_LIST from SAMPLE_SHEET

    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus_final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map4.sam"), file("${base}_map4.bam"), file("${base}.map4.sorted.bam"), file("${base}.map4.sorted.bam.bai"), file("${base}_final_mapping_stats.txt"), file("${base}_final_mapping_stats_map4.txt"), file("${base}.mpileup"), file("${base}_summary.csv"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_sample_id.txt"), file("${base}_pcr_ct.txt"), file("${base}_method.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt") from Final_Processing_PE_ch   

    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus_final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map4.sam"), file("${base}_map4.bam"), file("${base}.map4.sorted.bam"), file("${base}.map4.sorted.bam.bai"), file("${base}_final_mapping_stats.txt"), file("${base}_final_mapping_stats_map4.txt"), file("${base}.mpileup"), file("${base}_summary.csv"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_sample_id.txt"), file("${base}_pcr_ct.txt"), file("${base}_method.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt") into Serotyping_PE_ch

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
    serotype="serotype"
    printf ",\$Reads_On_Target" >> ${base}_final_summary.csv
    printf ",\$pcr_ct" >> ${base}_final_summary.csv
    printf ",\$method" >> ${base}_final_summary.csv
    printf ",\$NCBI_Name" >> ${base}_final_summary.csv
    printf ",\$serotype" >> ${base}_final_summary.csv    
    cp ${base}_final_summary.csv ${base}_summary.csv
    """
    }
    }
    }
    

if (params.withVapid) {
    if (params.singleEnd) {
process Vapid_Annotation {
    // container "docker.io/paulrkcruz/hrv-pl:latest"  
    // errorStrategy 'ignore'    
    // errorStrategy 'retry'
    // maxRetries 3
    
    input:
    file VAPID_DB_ALL_1
    file VAPID_DB_ALL_2
    file VAPID_DB_ALL_3
    file VAPID_DB_ALL_4
    file VAPID_DB_ALL_5
    file VAPID_DB_ALL_6
    file VAPID_DB_ALL_7
    file VAPID_DB_ALL_8
    file tbl2asn
    file vapid_python_main from vapid_python
    file vapid_python_main3 from vapid_python3
    file vapid_rhinovirus_sbt
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus_final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map4.sam"), file("${base}_map4.bam"), file("${base}.map4.sorted.bam"), file("${base}.map4.sorted.bam.bai"), file("${base}_final_mapping_stats.txt"), file("${base}_final_mapping_stats_map4.txt"), file("${base}.mpileup"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_sample_id.txt"), file("${base}_pcr_ct.txt"), file("${base}_method.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_collection_year.txt"), file("${base}_country_collected.txt"), file("${base}_blast_db_vp1.txt"), file("${base}_blast_db_all_ref.txt"), file("${base}_sample_stats.csv"), file("${base}_all_ref_id.txt"), file("${base}_nomen.txt") from Serotype_ch
    file("Run_Summary_cat.csv") from final_summary_out

    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus_final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map4.sam"), file("${base}_map4.bam"), file("${base}.map4.sorted.bam"), file("${base}.map4.sorted.bam.bai"), file("${base}_final_mapping_stats.txt"), file("${base}_final_mapping_stats_map4.txt"), file("${base}.mpileup"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_sample_id.txt"), file("${base}_pcr_ct.txt"), file("${base}_method.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_collection_year.txt"), file("${base}_country_collected.txt"), file("${base}_blast_db_vp1.txt"), file("${base}_blast_db_all_ref.txt"), file("${base}_sample_stats.csv"), file("${base}_all_ref_id.txt"), file ("${base}_vapid_metadata.csv") into All_files_ch
    // publishDir "${params.outdir}summary_vapid_annotation", mode: 'copy', pattern:'*_aligner.fasta*'
    // publishDir "${params.outdir}summary_vapid_annotation", mode: 'copy', pattern:'*_ref.fasta*'
    // publishDir "${params.outdir}summary_vapid_annotation", mode: 'copy', pattern:'*.ali*'
    // publishDir "${params.outdir}summary_vapid_annotation", mode: 'copy', pattern:'*.blastresults*'
    // publishDir "${params.outdir}summary_vapid_annotation", mode: 'copy', pattern:'*.fasta*'
    // publishDir "${params.outdir}summary_vapid_annotation", mode: 'copy', pattern:'*.fsa*'
    // publishDir "${params.outdir}summary_vapid_annotation", mode: 'copy', pattern:'*.gbf*'
    // publishDir "${params.outdir}summary_vapid_annotation", mode: 'copy', pattern:'*.sqn*'
    // publishDir "${params.outdir}summary_vapid_annotation", mode: 'copy', pattern:'*.tbl*'
    // publishDir "${params.outdir}summary_vapid_annotation", mode: 'copy', pattern:'*.val*'
    // publishDir "${params.outdir}summary_vapid_annotation", mode: 'copy', pattern:'*.cmt*'
    // publishDir "${params.outdir}summary_vapid_annotation", mode: 'copy', pattern:'*errorsummary.val*'
    // publishDir "${params.outdir}summary_vapid_annotation", mode: 'copy', pattern:'*_ref.gbk*'
    publishDir "${params.outdir}summary_vapid_metadata", mode: 'copy', pattern:'*_vapid_metadata.csv*'       

    script:

    """
    #!/bin/bash

    R1=${base}
    NCBI_Name=\${R1:4:6}
    SAMPLEName=\${R1:1:6}
    all_ref_id=\$(awk '{print \$1}' ${base}_all_ref_id.txt)


    # VIRAL ANNOTATION

    # Rhinovirus
    if grep -q \$all_ref_id "${base}_rv_ids.txt"; 
    then
    echo "< Accession found in Rhinovirus multifasta file. hrv_ref_rhinovirus.fa will be used for mapping."

    sample_id=\$(cat ${base}_sample_id.txt | sed -n '2 p')
    collection_year=\$(cat ${base}_collection_year.txt | sed -n '2 p')
    country_collected=\$(cat ${base}_country_collected.txt | sed -n '2 p')
    nomen=\$(awk 'FNR==1{print val,\$1}' ${base}_nomen.txt)	

    coverage=""
    echo strain, collection_date, country, coverage, full_name> ${base}_vapid_metadata.csv
    name=${base}
    printf "\$NCBI_Name, \$collection_year, \$country_collected, \$coverage, \$NCBI_Name" >> ${base}_vapid_metadata.csv

    cp ${base}.consensus_final.fa \${NCBI_Name}.final_consensus.fa
    
    seqkit replace -p '.+' -r \$NCBI_Name \${NCBI_Name}.final_consensus.fa > \${NCBI_Name}.fa

    reference_name=\$(grep -e ">" ${base}_mapped_ref_genome.fa > ${base}_ref_name.txt)
    ref_name_edit=\$(sed 's/>//g' ${base}_ref_name.txt > ${base}_ref_name_edit.txt)
    ref=\$(awk 'FNR==1{print val,\$1}' ${base}_ref_name_edit.txt)

    python3 ${vapid_python_main3} \${NCBI_Name}.fa ${vapid_rhinovirus_sbt} --metadata_loc ${base}_vapid_metadata.csv --output_location ${params.outdir}

	serotype=\$(awk 'FNR==1{print val,\$2}' ${base}_blast_db_vp1.txt)
    cut -d "-" -f2- <<< "\$serotype" > ${base}_serotype-parse.txt
    serotype_parsed=\$(awk 'FNR==1{print val,\$1}' ${base}_serotype-parse.txt)
    serotype_parse="\${serotype_parsed}" 
    echo \$serotype_parse > ${base}_sero.txt
    cat ${base}_sero.txt | tr -d " \t\n\r" > ${base}_serot.txt
    serotype_parsed2=\$(awk 'FNR==1{print val,\$1}' ${base}_serot.txt)
    space='/'
    nomenclature="\${serotype_parsed2} \${space} \${country_collected} \${space} \${collection_year} \${space} \${NCBI_Name}" 
    echo \$nomenclature > ${base}_nomenclature.txt
    cat ${base}_nomenclature.txt | tr -d " \t\n\r" > ${base}_nomenclature_parsed.txt
    nomenclature_parsed=\$(awk 'FNR==1{print val,\$1}' ${base}_nomenclature_parsed.txt)	
    echo \$serotype_parsed2 | xargs > ${base}_serots.txt
    echo \$nomenclature_parsed | xargs > ${base}_nomen.txt
    serots=\$(awk 'FNR==1{print val,\$1}' ${base}_serots.txt)	
    serots_adj=\$(awk 'FNR==1{print val,\$1}' ${base}_serots.txt  | xargs)	
    nomen=\$(awk 'FNR==1{print val,\$1}' ${base}_nomen.txt)	

    cp ${params.outdir}/summary_vapid_output/\${NCBI_Name}.sqn ${params.outdir}/summary_vapid_output/\${NCBI_Name}.txt

    # Edit Line:95 - taxname
    quote='"'
    rhino="rhinovirus \$serots_adj"
    taxname="              taxname \$quote"\$rhino"\$quote ,"
    sed "95s/.*/\$taxname/" ${params.outdir}/summary_vapid_output/\${NCBI_Name}.txt > ${params.outdir}/summary_vapid_output/\${NCBI_Name}_95.txt
    
    # Edit Line:100 - subname
    subname_1="\$quote\$nomen\$quote"
    subname_2="                   subname \$subname_1 } ,"
    sub=\$subname_2
    sed "100s/.*/"\$sub"/" 4N6KGN_95.txt > 4N6KGN_100.txt
    
    # Edit Line:188 - str - replace nomenclature with ncbi_name
    subname_1="\$quote\$nomen\$quote"
    subname_2="                   subname \$subname_1 } ,"
    sub=\$subname_2
    sed "100s/.*/"\$sub"/" 4N6KGN_95.txt > 4N6KGN_100.txt




    # HPV
    elif grep -q \$all_ref_id "${base}_hpv_ids.txt";
    then
    echo "< Accession found in HPV multifasta file. hrv_ref_hpv.fa will be used for mapping."




    # Influenza B
    elif grep -q \$all_ref_id "${base}_inbflb_ids.txt";
    then
    echo "< Accession found in Influenza B multifasta file. hrv_ref_Influenza_b.fa will be used for mapping."



    # Human Coronavirus
    elif grep -q \$all_ref_id "${base}_hcov_ids.txt";
    then
    echo "Accession found in HCoVs multifasta file. hrv_ref_hcov.fa will be used for mapping."




    # HPIV3 - Human parainfluenza virus 3
    elif grep -q \$all_ref_id "${base}_hpiv3.txt";
    then
    echo "Accession found in HPIV3 multifasta file. hrv_ref_hpiv3.fa will be used for mapping."



    else

    echo strain, collection_date, country, coverage, full_name> ${base}_vapid_metadata.csv

    fi



    """
    }
    } else {
process Vapid_Annotation_PE {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"        
    errorStrategy 'retry'
    maxRetries 3

    input:
    
    output:
    
    publishDir "${params.outdir}vapid_annotation_sqn", mode: 'copy', pattern:'*_final_summary.csv*'

    script:

    """
    #!/bin/bash

    """
    }
    }
    }
process Merge_run_summary {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"        
    errorStrategy 'retry'
    // maxRetries 3
    // errorStrategy 'ignore'

    input:
        file '*.csv' from Summary_cat_ch.collect()
        tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus_final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map4.sam"), file("${base}_map4.bam"), file("${base}.map4.sorted.bam"), file("${base}.map4.sorted.bam.bai"), file("${base}_final_mapping_stats.txt"), file("${base}_final_mapping_stats_map4.txt"), file("${base}.mpileup"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_sample_id.txt"), file("${base}_pcr_ct.txt"), file("${base}_method.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_collection_year.txt"), file("${base}_country_collected.txt"), file("${base}_blast_db_vp1.txt"), file("${base}_blast_db_all_ref.txt"), file("${base}_sample_stats.csv"), file("${base}_all_ref_id.txt"), file("${base}_nomen.txt") from Summary_file_ch

    output:
        file("Run_Summary_Final_cat.csv") into final_summary_out
        tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus_final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map4.sam"), file("${base}_map4.bam"), file("${base}.map4.sorted.bam"), file("${base}.map4.sorted.bam.bai"), file("${base}_final_mapping_stats.txt"), file("${base}_final_mapping_stats_map4.txt"), file("${base}.mpileup"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_sample_id.txt"), file("${base}_pcr_ct.txt"), file("${base}_method.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_collection_year.txt"), file("${base}_country_collected.txt"), file("${base}_blast_db_vp1.txt"), file("${base}_blast_db_all_ref.txt"), file("${base}_sample_stats.csv"), file("${base}_all_ref_id.txt"), file("${base}_nomen.txt") into Viral_annot_ch

    publishDir "${params.outdir}summary_withserotype_cat", mode: 'copy', pattern:'*Run_Summary_Final_cat.csv*'

    script:
    """
    #!/bin/bash
    R1=${base}
    NCBI_Name=\${R1:4:6}
    SAMPLEName=\${R1:2:5}
    all_ref_id=\$(awk '{print \$1}' ${base}_all_ref_id.txt)

    # Rhinovirus
    if grep -q \$all_ref_id "${base}_rv_ids.txt"; 
    then
    echo "< Accession found in Rhinovirus multifasta file. hrv_ref_rhinovirus.fa will be used for mapping."

    awk '(NR == 1) || (FNR > 1)' *.csv >  Run_Summary_cat.csv

    sed '1d' Run_Summary_cat.csv > Run_Summary_catted.csv

    echo -e "Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Reference_Genome,Reference_Length,Mapped_Reads,Percent_Ref_Coverage,Min_Coverage,Mean_Coverage,Max_Coverage,Bam_Size,Consensus_Length,Percent_N,%_Reads_On_Target, PCR_CT,Method, NCBI_Name, Serotype, Nomenclature, Reference_Name, Reference_Genome, Biosample_name, Biosample_accession, SRA_Accession, Release_date, Bioproject" | cat - Run_Summary_catted.csv > Run_Summary_Final_cat.csv


    # HPV
    elif grep -q \$all_ref_id "${base}_hpv_ids.txt";
    then
    echo "< Accession found in HPV multifasta file. hrv_ref_hpv.fa will be used for mapping."

    awk '(NR == 1) || (FNR > 1)' *.csv >  Run_Summary_cat.csv

    sed '1d' Run_Summary_cat.csv > Run_Summary_catted.csv

    echo -e "Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Reference_Genome,Reference_Length,Mapped_Reads,Percent_Ref_Coverage,Min_Coverage,Mean_Coverage,Max_Coverage,Bam_Size,Consensus_Length,Percent_N, Mapped_Reads_non-deduplicated, Mapped_Reads_Deduplicated, %_Reads_On_Target_deduplicated, %_Reads_On_Target_non-deduplicated, PCR_CT,Method, NCBI_Name, Reference_Name, Genotype, Genome, Biosample_name, Biosample_accession, SRA_Accession, Release_date, Bioproject" | cat - Run_Summary_catted.csv > Run_Summary_Final_cat.csv


    # Influenza B
    elif grep -q \$all_ref_id "${base}_inbflb_ids.txt";
    then
    echo "< Accession found in Influenza B multifasta file. hrv_ref_Influenza_b.fa will be used for mapping."

    awk '(NR == 1) || (FNR > 1)' *.csv >  Run_Summary_cat.csv

    sed '1d' Run_Summary_cat.csv > Run_Summary_catted.csv

    echo -e "Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Reference_Genome,Reference_Length,Mapped_Reads,Percent_Ref_Coverage,Min_Coverage,Mean_Coverage,Max_Coverage,Bam_Size,Consensus_Length,Percent_N,%_Reads_On_Target, PCR_CT,Method, NCBI_Name, Serotype, Nomenclature, Reference_Name, Reference_Genome, Biosample_name, Biosample_accession, SRA_Accession, Release_date, Bioproject" | cat - Run_Summary_catted.csv > Run_Summary_Final_cat.csv


    # Human Coronavirus
    elif grep -q \$all_ref_id "${base}_hcov_ids.txt";
    then
    echo "Accession found in HCoVs multifasta file. hrv_ref_hcov.fa will be used for mapping."

    awk '(NR == 1) || (FNR > 1)' *.csv >  Run_Summary_cat.csv

    sed '1d' Run_Summary_cat.csv > Run_Summary_catted.csv

    echo -e "Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Reference_Genome,Reference_Length,Mapped_Reads,Percent_Ref_Coverage,Min_Coverage,Mean_Coverage,Max_Coverage,Bam_Size,Consensus_Length,Percent_N,%_Reads_On_Target, PCR_CT,Method, NCBI_Name, Serotype, Nomenclature, Reference_Name, Reference_Genome, Biosample_name, Biosample_accession, SRA_Accession, Release_date, Bioproject" | cat - Run_Summary_catted.csv > Run_Summary_Final_cat.csv


    # HPIV3 - Human parainfluenza virus 3
    elif grep -q \$all_ref_id "${base}_hpiv3.txt";
    then
    echo "Accession found in HPIV3 multifasta file. hrv_ref_hpiv3.fa will be used for mapping."

    awk '(NR == 1) || (FNR > 1)' *.csv >  Run_Summary_cat.csv

    sed '1d' Run_Summary_cat.csv > Run_Summary_catted.csv

    echo -e "Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Reference_Genome,Reference_Length,Mapped_Reads,Percent_Ref_Coverage,Min_Coverage,Mean_Coverage,Max_Coverage,Bam_Size,Consensus_Length,Percent_N,%_Reads_On_Target, PCR_CT,Method, NCBI_Name, Serotype, Nomenclature, Reference_Name, Reference_Genome, Biosample_name, Biosample_accession, SRA_Accession, Release_date, Bioproject" | cat - Run_Summary_catted.csv > Run_Summary_Final_cat.csv

    else


    awk '(NR == 1) || (FNR > 1)' *.csv >  Run_Summary_cat.csv

    sed '1d' Run_Summary_cat.csv > Run_Summary_catted.csv

    echo -e "Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Reference_Genome,Reference_Length,Mapped_Reads,Percent_Ref_Coverage,Min_Coverage,Mean_Coverage,Max_Coverage,Bam_Size,Consensus_Length,Percent_N,%_Reads_On_Target, PCR_CT,Method, NCBI_Name, Serotype, Nomenclature, Reference_Name, Reference_Genome, Biosample_name, Biosample_accession, SRA_Accession, Release_date, Bioproject" | cat - Run_Summary_catted.csv > Run_Summary_Final_cat.csv

    fi


    """
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
	file '*_fastqc.{zip,html}' into Fastqc_results_SE

    publishDir "${params.outdir}fastqc_results", mode: 'copy', pattern:'*_fastqc.{zip,html}*'  

    script:
    """
    #!/bin/bash

    /usr/local/bin/fastqc --quiet --threads ${task.cpus} ${base}.trimmed.fastq.gz
    """
    }
} else {
process FastQC_PE {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"
	errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(base), file(R1),file(R2),file("${base}.R1.paired.fastq.gz"), file("${base}.R2.paired.fastq.gz"),file("${base}.R1.unpaired.fastq.gz"), file("${base}.R2.unpaired.fastq.gz") from Trim_out_fastqc_PE
    output: 
    file("*fastqc*") into Fastqc_results_PE 

    publishDir "${params.OUTDIR}fastqc", mode: 'copy'

    script:
    """
    #!/bin/bash

    /usr/local/bin/fastqc ${R1} ${R2} ${base}.R1.paired.fastq.gz ${base}.R2.paired.fastq.gz

    """
    }
}
}