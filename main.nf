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

HRV-Docker includes all dependencies. Currently (7/2021), Mapping step requires local dependencies. Please see docker for dependencies required.
    
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
version = '1.3'
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
    vapid_python = file("${baseDir}/vapid/vapid3.py")
    VAPID_DB_ALL_1 = file("${baseDir}/blast_db/all_virus.fasta")
    VAPID_DB_ALL_2 = file("${baseDir}/blast_db/all_virus.fasta.ndb")
    VAPID_DB_ALL_3 = file("${baseDir}/blast_db/all_virus.fasta.nhr")
    VAPID_DB_ALL_4 = file("${baseDir}/blast_db/all_virus.fasta.nin")
    VAPID_DB_ALL_5 = file("${baseDir}/blast_db/all_virus.fasta.nog")
    VAPID_DB_ALL_6 = file("${baseDir}/blast_db/all_virus.fasta.nos")
    VAPID_DB_ALL_7 = file("${baseDir}/blast_db/all_virus.fasta.not")
    VAPID_DB_ALL_8 = file("${baseDir}/blast_db/all_virus.fasta.nsq")

    tbl2asn = file("${baseDir}/blast_db/tbl2asn")  


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
// Script paths
TRIM_ENDS=file("${baseDir}/scripts/trim_ends.py")
VCFUTILS=file("${baseDir}/scripts/vcfutils.pl")
SPLITCHR=file("${baseDir}/scripts/splitchr.txt")
FIX_COVERAGE = file("${baseDir}/scripts/fix_coverage.py")
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
summary['Configuration Profile:'] = workflow.profile
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

    input: 
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_summary.csv") from Trim_out_SE
        file Reference_rv
        file Reference_hcov
        file Reference_hpv
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
/*
 * Call variants & Generate Consensus
  */
if (params.singleEnd) {
 process Generate_Consensus {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(base), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),val(bamsize),file("${base}.sorted.bam.bai"),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"), file("${base}_summary.csv"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") from Consensus_ch
    
    file VCFUTILS
    file SPLITCHR
    file TRIM_ENDS

    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus.fa"), file("${base}_final_summary.csv"), val(bamsize), val(id),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") into Consensus_Fasta_ch
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
    cp \${R1}_summary.csv \${R1}_final_summary.csv
    '''
}
} else {
 process Generate_Consensus_PE {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(base), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),val(bamsize),file("${base}.sorted.bam.bai"),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"), file("${base}_summary.csv"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") from Consensus_PE_ch
    
    file VCFUTILS
    file SPLITCHR
    file TRIM_ENDS

    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus.fa"), file("${base}_final_summary.csv"), val(bamsize), val(id),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") into Consensus_Fasta_PE_ch
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
    cp \${R1}_summary.csv \${R1}_final_summary.csv
    '''
}
}
if (params.singleEnd) {
process Final_Mapping {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"     
	errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(base), file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus.fa"), file("${base}_final_summary.csv"), val(bamsize), val(id),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") from Consensus_Fasta_ch

    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus_final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map4.sam"), file("${base}_map4.bam"), file("${base}.map4.sorted.bam"), file("${base}.map4.sorted.bam.bai"), file("${base}_final_mapping_stats.txt"), file("${base}_final_mapping_stats_map4.txt"), file("${base}.mpileup"), file("${base}_summary.csv"), file("${base}.trimmed.fastq.gz"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt") into Mapping_Final_ch

    publishDir "${params.outdir}mpileup_map3", mode: 'copy', pattern:'*.mpileup*'
    publishDir "${params.outdir}bam_map3", mode: 'copy', pattern:'*.map3.sorted.bam*'
    publishDir "${params.outdir}sam_map3", mode: 'copy', pattern:'*_map3.sam*'
    publishDir "${params.outdir}bam_map4", mode: 'copy', pattern:'*.map4.sorted.bam*'
    publishDir "${params.outdir}sam_map4", mode: 'copy', pattern:'*_map4.sam*'
    publishDir "${params.outdir}consensus-final", mode: 'copy', pattern:'*consensus_final.fa*'
    publishDir "${params.outdir}consensus-final", mode: 'copy', pattern:'*consensus_final.fa*'
    // publishDir "${params.outdir}consensus-ivar-masked", mode: 'copy', pattern:'*.consensus.masked.fa*'
    publishDir "${params.outdir}txt_bbmap_final_mapping_stats_map3", mode: 'copy', pattern:'*_final_mapping_stats.txt*'
    publishDir "${params.outdir}txt_bbmap_final_mapping_stats_map4", mode: 'copy', pattern:'*_final_mapping_stats_map4.txt*'
    publishDir "${params.outdir}summary", mode: 'copy', pattern:'*_summary.csv*'

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
    echo "< Accession found in respiratory virus multifasta file. hrv_ref_hpv.fa will be used for mapping."

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

    sed 's/>.*/>${base}/' ${base}.consensusfinal.fa > ${base}.consensusfinal-renamed-header.fa
    grep -v "^>" ${base}.consensusfinal-renamed-header.fa | tr -cd N | wc -c > N.txt
    cp ${base}.consensusfinal-renamed-header.fa ${base}.consensus_final.fa
    
    num_ns=\$(awk 'FNR==1{print val,\$1}' N.txt)
    echo "\$num_ns/\$num_bases*100" | bc -l > n_percent.txt
    percent_n=\$(awk 'FNR==1{print val,\$1}' n_percent.txt)
    printf ",\$num_bases" >> ${base}_final_summary.csv
    printf ",\$percent_n" >> ${base}_final_summary.csv
    cp ${base}_final_summary.csv ${base}_summary.csv

    # FINAL MAPPING: map4

    # Rhinovirus
    if grep -q \$all_ref_id "${base}_rv_ids.txt"; 
    then
    echo "< Accession found in Rhinovirus multifasta file. hrv_ref_rhinovirus.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map4.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats_map4.txt 2>&1

    # HPV
    elif grep -q \$all_ref_id "${base}_hpv_ids.txt";
    then
    echo "< Accession found in HPV multifasta file. hrv_ref_hpv.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map4.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats_map4.txt 2>&1

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
    maxRetries 3

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
    // errorStrategy 'retry'
    // maxRetries 3

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
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus_final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map4.sam"), file("${base}_map4.bam"), file("${base}.map4.sorted.bam"), file("${base}.map4.sorted.bam.bai"), file("${base}_final_mapping_stats.txt"), file("${base}_final_mapping_stats_map4.txt"), file("${base}.mpileup"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_sample_id.txt"), file("${base}_pcr_ct.txt"), file("${base}_method.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_collection_year.txt"), file("${base}_country_collected.txt"), file("${base}_blast_db_vp1.txt"), file("${base}_blast_db_all_ref.txt"), file("${base}_sample_stats.csv"), file("${base}_all_ref_id.txt") into Serotype_ch
    tuple val(base), file("${base}_vapid_metadata.csv") into Vapid_metadata_cat_ch
    tuple val(base), file("${base}_summary_final.csv") into Summary_cat_ch

    publishDir "${params.outdir}blast_serotype", mode: 'copy', pattern:'*_blast_db_vp1.txt*'
    publishDir "${params.outdir}blast_ref_genome", mode: 'copy', pattern:'*_blast_db_all_ref.txt*'    
    publishDir "${params.outdir}sample_Stats", mode: 'copy', pattern:'*_sample_stats.csv*'
    publishDir "${params.outdir}sample_Stats", mode: 'copy', pattern:'*.txt*'       
    publishDir "${params.outdir}sample_Stats", mode: 'copy', pattern:'*_collection_year.txt*'
    publishDir "${params.outdir}sample_Stats", mode: 'copy', pattern:'*_country_collected.txt*'
    publishDir "${params.outdir}sample_Stats", mode: 'copy', pattern:'*_nomenclature.txt*'
    publishDir "${params.outdir}summary_vapid_metadata", mode: 'copy', pattern:'*_vapid_metadata.csv*'
    publishDir "${params.outdir}summary_withserotype", mode: 'copy', pattern:'*_summary_final.csv*'

    script:

    """
    #!/bin/bash
    R1=${base}
    NCBI_Name=\${R1:4:6}
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
    cut -d "-" -f2- <<< "\$serotype" > serotype-parse.txt
    serotype_parsed=\$(awk 'FNR==1{print val,\$1}' serotype-parse.txt)
    rv='Rv'
    serotype_parse="\${rv} \${serotype_parsed}" 
    echo \$serotype_parse > sero.txt
    cat sero.txt | tr -d " \t\n\r" > serot.txt
    serotype_parsed2=\$(awk 'FNR==1{print val,\$1}' serot.txt)
    space='/'
    nomenclature="\${serotype_parsed2} \${space} \${country_collected} \${space} \${collection_year} \${space} \${NCBI_Name}" 
    echo \$nomenclature > nomenclature.txt
    cat nomenclature.txt | tr -d " \t\n\r" > nomenclature_parsed.txt
    nomenclature_parsed=\$(awk 'FNR==1{print val,\$1}' nomenclature_parsed.txt)	
    echo \$serotype_parsed2 | xargs > serots.txt
    echo \$nomenclature_parsed | xargs > nomen.txt
    serots=\$(awk 'FNR==1{print val,\$1}' serots.txt)	
    nomen=\$(awk 'FNR==1{print val,\$1}' nomen.txt)	

    blastn -out ${base}_blast_db_all_ref.txt -query ${base}.consensus_final.fa -db ${BLASTDB_ALL_1} -outfmt "5 std qlen" -task blastn -max_target_seqs 1 -evalue 1e-5

    awk 'NR==31' ${base}_blast_db_all_ref.txt > strain.txt
    sed -i -e 's/<Hit_def>//g'  strain.txt
    awk -F'</Hit_def>' '{print \$1}' strain.txt | xargs > strain-parsed.txt
    Reference_Name=\$(head -n 1 strain-parsed.txt)

    biosample_name=\$(cat ${base}_biosample_name.txt | sed -n '2 p')
    biosample_accession=\$(cat ${base}_biosample_accession.txt | sed -n '2 p')
    sra_accession=\$(cat ${base}_sra_accession.txt | sed -n '2 p')
    release_date=\$(cat ${base}_release_date.txt | sed -n '2 p')
    bioproject=\$(cat ${base}_bioproject.txt | sed -n '2 p')
    
    # Respiratory Panel
    elif grep -q \$all_ref_id "${base}_hpv_ids.txt";
    then
    echo "< Accession found in respiratory virus multifasta file."

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
    cut -d "-" -f2- <<< "\$serotype" > serotype-parse.txt
    serotype_parsed=\$(awk 'FNR==1{print val,\$1}' serotype-parse.txt)
    rv='HCoV_'
    serotype_parse="\${rv} \${serotype_parsed}" 
    echo \$serotype_parse > sero.txt
    cat sero.txt | tr -d " \t\n\r" > serot.txt
    serotype_parsed2=\$(awk 'FNR==1{print val,\$1}' serot.txt)
    space='/'
    nomenclature="\${serotype_parsed2} \${space} \${country_collected} \${space} \${collection_year} \${space} \${NCBI_Name}" 
    echo \$nomenclature > nomenclature.txt
    cat nomenclature.txt | tr -d " \t\n\r" > nomenclature_parsed.txt
    nomenclature_parsed=\$(awk 'FNR==1{print val,\$1}' nomenclature_parsed.txt)	
    echo \$serotype_parsed2 | xargs > serots.txt
    echo \$nomenclature_parsed | xargs > nomen.txt
    serots=\$(awk 'FNR==1{print val,\$1}' serots.txt)	
    nomen=\$(awk 'FNR==1{print val,\$1}' nomen.txt)	

    blastn -out ${base}_blast_db_all_ref.txt -query ${base}.consensus_final.fa -db ${BLASTDB_ALL_1} -outfmt "5 std qlen" -task blastn -max_target_seqs 1 -evalue 1e-5

    awk 'NR==31' ${base}_blast_db_all_ref.txt > strain.txt
    sed -i -e 's/<Hit_def>//g'  strain.txt
    awk -F'</Hit_def>' '{print \$1}' strain.txt | xargs > strain-parsed.txt
    Reference_Name=\$(head -n 1 strain-parsed.txt)

    biosample_name=\$(cat ${base}_biosample_name.txt | sed -n '2 p')
    biosample_accession=\$(cat ${base}_biosample_accession.txt | sed -n '2 p')
    sra_accession=\$(cat ${base}_sra_accession.txt | sed -n '2 p')
    release_date=\$(cat ${base}_release_date.txt | sed -n '2 p')
    bioproject=\$(cat ${base}_bioproject.txt | sed -n '2 p')

    # Influenza B
    elif grep -q \$all_ref_id "${base}_inbflb_ids.txt";
    then
    echo "< Accession found in Influenza B multifasta file."

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
    cut -d "-" -f2- <<< "\$serotype" > serotype-parse.txt
    serotype_parsed=\$(awk 'FNR==1{print val,\$1}' serotype-parse.txt)
    rv='Influenza_'
    serotype_parse="\${rv} \${serotype_parsed}" 
    echo \$serotype_parse > sero.txt
    cat sero.txt | tr -d " \t\n\r" > serot.txt
    serotype_parsed2=\$(awk 'FNR==1{print val,\$1}' serot.txt)
    space='/'
    nomenclature="\${serotype_parsed2} \${space} \${country_collected} \${space} \${collection_year} \${space} \${NCBI_Name}" 
    echo \$nomenclature > nomenclature.txt
    cat nomenclature.txt | tr -d " \t\n\r" > nomenclature_parsed.txt
    nomenclature_parsed=\$(awk 'FNR==1{print val,\$1}' nomenclature_parsed.txt)	
    echo \$serotype_parsed2 | xargs > serots.txt
    echo \$nomenclature_parsed | xargs > nomen.txt
    serots=\$(awk 'FNR==1{print val,\$1}' serots.txt)	
    nomen=\$(awk 'FNR==1{print val,\$1}' nomen.txt)	

    blastn -out ${base}_blast_db_all_ref.txt -query ${base}.consensus_final.fa -db ${BLASTDB_ALL_1} -outfmt "5 std qlen" -task blastn -max_target_seqs 1 -evalue 1e-5

    awk 'NR==31' ${base}_blast_db_all_ref.txt > strain.txt
    sed -i -e 's/<Hit_def>//g'  strain.txt
    awk -F'</Hit_def>' '{print \$1}' strain.txt | xargs > strain-parsed.txt
    Reference_Name=\$(head -n 1 strain-parsed.txt)

    biosample_name=\$(cat ${base}_biosample_name.txt | sed -n '2 p')
    biosample_accession=\$(cat ${base}_biosample_accession.txt | sed -n '2 p')
    sra_accession=\$(cat ${base}_sra_accession.txt | sed -n '2 p')
    release_date=\$(cat ${base}_release_date.txt | sed -n '2 p')
    bioproject=\$(cat ${base}_bioproject.txt | sed -n '2 p')

    # Human Coronavirus
    elif grep -q \$all_ref_id "${base}_hcov_ids.txt";
    then
    echo "Accession found in HCoVs multifasta file."

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
    cut -d "-" -f2- <<< "\$serotype" > serotype-parse.txt
    serotype_parsed=\$(awk 'FNR==1{print val,\$1}' serotype-parse.txt)
    rv='HCoV_'
    serotype_parse="\${rv} \${serotype_parsed}" 
    echo \$serotype_parse > sero.txt
    cat sero.txt | tr -d " \t\n\r" > serot.txt
    serotype_parsed2=\$(awk 'FNR==1{print val,\$1}' serot.txt)
    space='/'
    nomenclature="\${serotype_parsed2} \${space} \${country_collected} \${space} \${collection_year} \${space} \${NCBI_Name}" 
    echo \$nomenclature > nomenclature.txt
    cat nomenclature.txt | tr -d " \t\n\r" > nomenclature_parsed.txt
    nomenclature_parsed=\$(awk 'FNR==1{print val,\$1}' nomenclature_parsed.txt)	
    echo \$serotype_parsed2 | xargs > serots.txt
    echo \$nomenclature_parsed | xargs > nomen.txt
    serots=\$(awk 'FNR==1{print val,\$1}' serots.txt)	
    nomen=\$(awk 'FNR==1{print val,\$1}' nomen.txt)	

    blastn -out ${base}_blast_db_all_ref.txt -query ${base}.consensus_final.fa -db ${BLASTDB_ALL_1} -outfmt "5 std qlen" -task blastn -max_target_seqs 1 -evalue 1e-5

    awk 'NR==31' ${base}_blast_db_all_ref.txt > strain.txt
    sed -i -e 's/<Hit_def>//g'  strain.txt
    awk -F'</Hit_def>' '{print \$1}' strain.txt | xargs > strain-parsed.txt
    Reference_Name=\$(head -n 1 strain-parsed.txt)

    biosample_name=\$(cat ${base}_biosample_name.txt | sed -n '2 p')
    biosample_accession=\$(cat ${base}_biosample_accession.txt | sed -n '2 p')
    sra_accession=\$(cat ${base}_sra_accession.txt | sed -n '2 p')
    release_date=\$(cat ${base}_release_date.txt | sed -n '2 p')
    bioproject=\$(cat ${base}_bioproject.txt | sed -n '2 p')

    # HPIV3 - Human parainfluenza virus 3
    elif grep -q \$all_ref_id "${base}_hpiv3.txt";
    then
    echo "Accession found in HPIV3 multifasta file."

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
    cut -d "-" -f2- <<< "\$serotype" > serotype-parse.txt
    serotype_parsed=\$(awk 'FNR==1{print val,\$1}' serotype-parse.txt)
    rv='HCoV_'
    serotype_parse="\${rv} \${serotype_parsed}" 
    echo \$serotype_parse > sero.txt
    cat sero.txt | tr -d " \t\n\r" > serot.txt
    serotype_parsed2=\$(awk 'FNR==1{print val,\$1}' serot.txt)
    space='/'
    nomenclature="\${serotype_parsed2} \${space} \${country_collected} \${space} \${collection_year} \${space} \${NCBI_Name}" 
    echo \$nomenclature > nomenclature.txt
    cat nomenclature.txt | tr -d " \t\n\r" > nomenclature_parsed.txt
    nomenclature_parsed=\$(awk 'FNR==1{print val,\$1}' nomenclature_parsed.txt)	
    echo \$serotype_parsed2 | xargs > serots.txt
    echo \$nomenclature_parsed | xargs > nomen.txt
    serots=\$(awk 'FNR==1{print val,\$1}' serots.txt)	
    nomen=\$(awk 'FNR==1{print val,\$1}' nomen.txt)	

    blastn -out ${base}_blast_db_all_ref.txt -query ${base}.consensus_final.fa -db ${BLASTDB_ALL_1} -outfmt "5 std qlen" -task blastn -max_target_seqs 1 -evalue 1e-5

    awk 'NR==31' ${base}_blast_db_all_ref.txt > strain.txt
    sed -i -e 's/<Hit_def>//g'  strain.txt
    awk -F'</Hit_def>' '{print \$1}' strain.txt | xargs > strain-parsed.txt
    Reference_Name=\$(head -n 1 strain-parsed.txt)

    biosample_name=\$(cat ${base}_biosample_name.txt | sed -n '2 p')
    biosample_accession=\$(cat ${base}_biosample_accession.txt | sed -n '2 p')
    sra_accession=\$(cat ${base}_sra_accession.txt | sed -n '2 p')
    release_date=\$(cat ${base}_release_date.txt | sed -n '2 p')
    bioproject=\$(cat ${base}_bioproject.txt | sed -n '2 p')
    
    else

    printf ",\$serotype" >> ${base}_final_summary.csv   
    printf ",\$nomenclature" >> ${base}_final_summary.csv 

    fi
   
    coverage=""
    echo strain, collection_date, country, coverage, full_name> ${base}_vapid_metadata.csv
    name=${base}
    printf "\$name, \$collection_year, \$country_collected, \$coverage, \$nomen" >> ${base}_vapid_metadata.csv

    printf ",\$serots" >> ${base}_final_summary.csv
    printf ",\$nomen" >> ${base}_final_summary.csv
    printf ",\$Reference_Name" >> ${base}_final_summary.csv
    printf ",\$biosample_name" >> ${base}_final_summary.csv
    printf ",\$biosample_accession" >> ${base}_final_summary.csv
    printf ",\$sra_accession" >> ${base}_final_summary.csv
    printf ",\$release_date" >> ${base}_final_summary.csv
    printf ",\$bioproject" >> ${base}_final_summary.csv
    cp ${base}_final_summary.csv ${base}_summary_final.csv

    """
    }

process Merge_run_summary {
    echo true

    input:
        file '*.csv' from Summary_cat_ch.collect()

    output:
        file("Run_Summary_Final_cat.csv") into final_summary_out

    publishDir "${params.outdir}summary_withserotype_cat", mode: 'copy', pattern:'*Run_Summary_Final_cat.csv*'

    script:
    """
    
    awk '(NR == 1) || (FNR > 1)' *.csv >  Run_Summary_cat.csv

    sed '1d' Run_Summary_cat.csv > Run_Summary_catted.csv

    echo -e "Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Reference_Genome,Reference_Length,Mapped_Reads,Percent_Ref_Coverage,Min_Coverage,Mean_Coverage,Max_Coverage,Bam_Size,Consensus_Length,Percent_N,%_Reads_On_Target, PCR_CT,Method, NCBI_Name, Serotype, Nomenclature, Reference_Name, Reference_Genome, Biosample_name, Biosample_accession, SRA_Accession, Release_date, Bioproject" | cat - Run_Summary_catted.csv > Run_Summary_Final_cat.csv

    """
    }

// process Merge_vapid_metadata {
//     echo true

//     input:
//         file '*.csv' from Vapid_metadata_cat_ch.collect()

//     output:
//         file("vapid_metadata_cat.csv") into final_metadata_out

//     publishDir "${params.outdir}summary_vapid_metadata_cat", mode: 'copy', pattern:'*vapid_metadata_cat.csv*'    

//     script:
//     """

//     awk '(NR == 1) || (FNR > 1)' *.csv >  vapid_metadata_cat.csv


//     """
//     }

    // echo strain, collection_date, country, coverage, full_name> vapid_metadata_cat.csv
    // echo Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Reference_Genome,Reference_Length,Mapped_Reads,Percent_Ref_Coverage,Min_Coverage,Mean_Coverage,Max_Coverage,Bam_Size,Consensus_Length,Percent_N,%_Reads_On_Target, PCR_CT,Method, NCBI_Name, Serotype, Nomenclature, Reference_Name, Reference_Genome, Biosample_name, Biosample_accession, SRA_Accession, Release_date, Bioproject> Run_Summary_cat.csv

    // sed -i '1s/^/strain, collection_date, country, coverage, full_name\n/' vapid_metadata_cat.csv

    // sed -i '1s/^/Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Reference_Genome,Reference_Length,Mapped_Reads,Percent_Ref_Coverage,Min_Coverage,Mean_Coverage,Max_Coverage,Bam_Size,Consensus_Length,Percent_N,%_Reads_On_Target, PCR_CT,Method, NCBI_Name, Serotype, Nomenclature, Reference_Name, Reference_Genome, Biosample_name, Biosample_accession, SRA_Accession, Release_date, Bioproject\n/' Run_Summary_cat.csv

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
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"
    container "docker.io/confurious/blastn:latest"     
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
    file vapid_rhinovirus_sbt
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
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus_final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map4.sam"), file("${base}_map4.bam"), file("${base}.map4.sorted.bam"), file("${base}.map4.sorted.bam.bai"), file("${base}_final_mapping_stats.txt"), file("${base}_final_mapping_stats_map4.txt"), file("${base}.mpileup"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_sample_id.txt"), file("${base}_pcr_ct.txt"), file("${base}_method.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_collection_year.txt"), file("${base}_country_collected.txt"), file("${base}_blast_db_vp1.txt"), file("${base}_blast_db_all_ref.txt"), file("${base}_sample_stats.csv"), file("${base}_all_ref_id.txt") from Serotype_ch
    tuple val(base), file("${base}_vapid_metadata.csv") from Vapid_metadata_cat_ch
    // tuple val(base), file("${base}_summary_final.csv") from Summary_cat_ch
    file("Run_Summary_cat.csv") from final_summary_out
    // file("vapid_metadata_cat.csv") from final_metadata_out

    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}_most_mapped_ref.txt"), file("${base}.consensus_final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map4.sam"), file("${base}_map4.bam"), file("${base}.map4.sorted.bam"), file("${base}.map4.sorted.bam.bai"), file("${base}_final_mapping_stats.txt"), file("${base}_final_mapping_stats_map4.txt"), file("${base}.mpileup"), val(bamsize), val(id), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_sample_id.txt"), file("${base}_pcr_ct.txt"), file("${base}_method.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_collection_year.txt"), file("${base}_country_collected.txt"), file("${base}_blast_db_vp1.txt"), file("${base}_blast_db_all_ref.txt"), file("${base}_sample_stats.csv"), file("${base}_all_ref_id.txt") into All_files_ch
    tuple val(base), file("${base}_vapid_metadata.csv"), file("${base}.sqn") into Vapid_md_output_ch
    tuple val(base), file("${base}_summary_final.csv") into Summary_outputs_ch

    publishDir "${params.outdir}Summary_vapid_annotation_sqn", mode: 'copy', pattern:'*.sqn*'
    publishDir "${params.outdir}Summary_merged", mode: 'copy', pattern:'*Run_Summary_Final.csv*'

    script:

    """
    #!/bin/bash

    reference_name=\$(grep -e ">" ${base}_mapped_ref_genome.fa > ${base}_ref_name.txt)
    ref_name_edit=\$(sed 's/>//g' ${base}_ref_name.txt > ${base}_ref_name_edit.txt)
    ref=\$(awk 'FNR==1{print val,\$1}' ${base}_ref_name_edit.txt)

    python3 ${vapid_python_main} ${base}.consensus_final.fa ${vapid_rhinovirus_sbt} --r \${ref} --metadata_loc ${base}_vapid_metadata.csv


    

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