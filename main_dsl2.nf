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
Updated: October 27, 2021
Revision: v1.4
LICENSE: GNU
----------------------------------------------------------------------------------------
*/
// Nextflow dsl2
nextflow.enable.dsl=2
// Pipeline version
version = '1.4'
params.helpMsg = false
// Setup pipeline help msg
def helpMsg() {
    log.info"""
    _______________________________________________________________________________
    Human Respiratory Virus Pipeline :  Version ${version}
	________________________________________________________________________________
    
	Pipeline Usage:
    To run the pipeline, enter the following in the command line:
        nextflow run FILE_PATH/HRV_Genome_Mapping_Pipeline/main.nf --reads PATH_TO_FASTQ --outdir PATH_TO_OUTPUT_DIR
    Valid CLI Arguments:
    REQUIRED:
        --reads                       Path to input fastq.gz folder).
        --outdir                      The output directory where the results will be saved
        --singleEnd                   Specifies that the input fastq files are single end reads    
    OPTIONAL:
        --withMetadata                Adds Metadata information to Final Run Report Summary
        --withVapid                   Annotate the resulting consensus fasta for GenBank submission         
        --ref_rv                      Overwrite set multi-fasta Rhinovirus reference file
        --ref_hcov                    Overwrite set multi-fasta Human Coronavirus reference file
        --ref_hpv                     Overwrite set multi-fasta HPV reference file
        --ref_inflb                   Overwrite set multi-fasta Influenza B reference file
        --ref_hpiv3                   Overwrite set multi-fasta HPIV3 reference file
        --helpMsg					  Displays help message in terminal
        --singleEnd                   Specifies that the input fastq files are single end reads
        --withFastQC					Runs a quality control check on fastq files  
    """.stripIndent()
}
// Show help msg
if (params.helpMsg){
    helpMsg()
    exit 0
}
// Setup command-line parameters
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
params.ref_rv = false
params.ref_hcov = false
params.ref_hpv = false
params.ref_inflb = false
params.ref_hpiv3 = false
// Make sure INPUT ends with trailing slash
if (!params.reads.endsWith("/")){
   params.reads = "${params.reads}/"
}
// if OUTDIR not set
if (params.outdir == false) {
    println( "Must provide an output directory with --outdir") 
    exit(1)
}
// Make sure OUTDIR ends with trailing slash
if (!params.outdir.endsWith("/")){
   params.outdir = "${params.outdir}/"
}
// Check Nextflow version
nextflow_req_v = '20.10.0'
try {
    if( ! nextflow.version.matches(">= $nextflow_req_v") ){
        throw GroovyException("> ERROR: The version of Nextflow running on your machine is out dated.\n>Please update to Version $nextflow_req_v")
    }
} catch (all) {
	log.error"ERROR: This version of Nextflow is out of date.\nPlease update to the latest version of Nextflow."
}
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
// Setup metadata
if(params.withMetadata != false) {
    METADATA = file(params.withMetadata)
}
// File Paths
BLAST_DB_VP1 = file("${baseDir}/blast_db/VP1_164_annotated_nospaces.fasta")
BLAST_DB_ALL = file("${baseDir}/blast_db/allref.fasta")
BBMAP_PATH="/Users/greningerlab/Documents/bbmap/"
params.ADAPTERS_SE = file("${baseDir}/adapters/TruSeq2-SE.fa")
params.ADAPTERS_EE = file("${baseDir}/adapters/TruSeq2-PE.fa")
ADAPTERS_SE = file("${baseDir}/adapters/TruSeq2-SE.fa")
ADAPTERS_PE = file("${baseDir}/adapters/TruSeq2-PE.fa")
vapid_rhinovirus_sbt = file("${baseDir}/vapid/rhinovirus.sbt")
params.SETTING = "2:30:10:1:true"
params.LEADING = "3"
params.TRAILING = "3"
params.SWINDOW = "4:20"
params.MINLEN = "35"
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
// Setup hrv pipeline header
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
summary['Metadata:']               = params.withMetadata ? params.withMetadata : 'False'
summary['Serotyping:']               = params.withSerotype ? 'True' : 'False'
summary['Viral Annotation:']               = params.withVapid ? 'True' : 'False'
summary['Sequence type:']           	  = params.singleEnd ? 'Single-End' : 'Paired-End'
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

//
// Import processes
//
include { Trimming } from './modules_dsl2.nf'
include { Aligning } from './modules_dsl2.nf'
include { Bam_Sorting } from './modules_dsl2.nf'
include { Consensus_Generation } from './modules_dsl2.nf'
include { Aligning_Final } from './modules_dsl2.nf'
include { Final_Processing } from './modules_dsl2.nf'
include { Summary_Generation } from './modules_dsl2.nf'
include { Serotyping } from './modules_dsl2.nf'
include { Vapid_Annotation } from './modules_dsl2.nf'
include { FastQC } from './modules_dsl2.nf'

// Create channel for input reads: single-end or paired-end
if(params.singleEnd == false) {
    // Check for R1s and R2s in input directory
    input_read_ch = Channel
        .fromFilePairs("${params.reads}*_R{1,2}*.gz")
        .ifEmpty { error "Cannot find any FASTQ pairs in ${params.reads} ending with .gz" }
        .map { it -> [it[0], it[1][0], it[1][1]]}
} else {
    // Looks for gzipped files, assumes all separate samples
    input_read_ch = Channel
        .fromPath("${params.reads}*.gz")
        //.map { it -> [ file(it)]}
        .map { it -> file(it)}
}

////////////////////////////////////////////////////////
////////////////////////////////////////////////////////
/*                                                    */
/*                 RUN THE WORKFLOW                   */
/*                                                    */
////////////////////////////////////////////////////////
////////////////////////////////////////////////////////

workflow {
    if(params.singleEnd == true) {
    Trimming (
        input_read_ch, 
        ADAPTERS_SE,
        params.MINLEN,
        params.SETTING, 
        params.LEADING,
        params.TRAILING,
        params.SWINDOW,
    )
    Aligning (
        Trimming.out[0],
        Reference_rv,
        Reference_hcov,
        Reference_hpv,
        Reference_hpv_14,        
        Reference_inflb,
        Reference_hpiv3
    )
    Bam_Sorting (
        Aligning.out[0]
    )
    Consensus_Generation (
        Bam_Sorting.out[0]
    )
    Aligning_Final (
        Consensus_Generation.out[0]
    )
    if(params.withMetadata == true) {
    Summary_Generation (
        Aligning_Final.out[0],
        METADATA
    )
    }
    if(params.withSerotype == true) {
    Serotyping (
        Summary_Generation.out[0],
        Summary_Generation.out[1],
        METADATA,
        BLAST_DB_VP1_1,
        BLAST_DB_VP1_2,
        BLAST_DB_VP1_3,
        BLAST_DB_VP1_4,
        BLAST_DB_VP1_5,
        BLAST_DB_VP1_6,
        BLAST_DB_VP1_7,
        BLAST_DB_VP1_8,
        BLAST_DB_ALL_1,
        BLAST_DB_ALL_2,
        BLAST_DB_ALL_3,
        BLAST_DB_ALL_4,
        BLAST_DB_ALL_5,
        BLAST_DB_ALL_6,
        BLAST_DB_ALL_7,
        BLAST_DB_ALL_8,
        BLAST_DB_ALL_9,
        BLAST_DB_ALL_10
    )
    Final_Processing (
        Serotyping.out[0].collect(),
        Summary_Generation.out[0]
    )
    }
    if(params.withVapid == true) {
    Vapid_Annotation (
        Serotyping.out[0]
        VAPID_DB_ALL_1,
        VAPID_DB_ALL_2,
        VAPID_DB_ALL_3,
        VAPID_DB_ALL_4,
        VAPID_DB_ALL_5,
        VAPID_DB_ALL_6,
        VAPID_DB_ALL_7,
        VAPID_DB_ALL_8,
        tbl2asn,
        vapid_python,
        vapid_python3,
        vapid_rhinovirus_sbt
    )
    }   
    } else {
        // paired-end workflow
    }
}