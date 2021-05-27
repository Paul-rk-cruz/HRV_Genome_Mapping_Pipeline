#!/usr/bin/env nextflow

/*
========================================================================================
                  Rhinovirus Genome Mapping Pipeline v1.0
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
params.withFastQC = false
params.skipTrim = false
params.reads = false
params.singleEnd = false
params.ADAPTERS = false
// Script Files
TRIM_ENDS=file("${baseDir}/scripts/trim_ends.py")
VCFUTILS=file("${baseDir}/scripts/vcfutils.pl")
SPLITCHR=file("${baseDir}/scripts/splitchr.txt")
FIX_COVERAGE = file("${baseDir}/scripts/fix_coverage.py")
ADAPTERS_SE = file("${baseDir}/adapters/TruSeq2-SE.fa")
ADAPTERS_PE = file("${baseDir}/adapters/TruSeq2-PE.fa")
REFERENCE_FASTA = file("${baseDir}/hrv_ref/hrv_ref_db01_accession_only.fa")
REFERENCE_FASTA_INDEX = file("${baseDir}/hrv_ref/hrv_ref_db01.fa.fai")
BBMAP_PATH="/Users/uwvirongs/Documents/KC/bbmap/"
params.trimmomatic_adapters_file_PE = "/Users/uwvirongs/miniconda3/share/trimmomatic-0.39-2/adapters/TruSeq2-PE.fa"
params.trimmomatic_adapters_file_SE = "/Users/uwvirongs/miniconda3/share/trimmomatic-0.39-2/adapters/TruSeq2-SE.fa"
params.trimmomatic_adapters_parameters = "2:30:10:1"
params.trimmomatic_window_length = "4"
params.trimmomatic_window_value = "20"
params.MINLEN = "35"
MINLEN = "35"
// Show help msg
if (params.helpMsg){
    helpMsg()
    exit 0
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

if (! params.reads ) exit 1, "> Error: Fastq files not found. Please specify a valid path with --reads"
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
// summary['Trimmomatic adapters:'] = params.trimmomatic_adapters_file_SE
} else {
// summary['Trimmomatic adapters:'] = params.trimmomatic_adapters_file_PE
}
summary["Trimmomatic read length (minimum):"] = params.MINLEN
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

/*
 * Trim Reads
 * 
 * Processing: Trim adaptors and repetitive bases from sequence reads and remove low quality sequence reads.
 */
// if(!params.skipTrimming) {
if (params.singleEnd) {
	process Trim_Reads_SE {
    // container "quay.io/biocontainers/trimmomatic:0.35--6"
    errorStrategy 'retry'
    maxRetries 3

    input:
        file R1 from input_read_ch
        file ADAPTERS_SE
        val MINLEN

    output: 
        tuple env(base),file("*.trimmed.fastq.gz"),file("*summary.csv") into Trim_out_SE, Trim_out_SE_FQC

    publishDir "${params.outdir}trimmed_fastqs", mode: 'copy',pattern:'*.trimmed.fastq*'


    script:
    """
    #!/bin/bash
    base=`basename ${R1} ".fastq.gz"`
    echo \$base
    trimmomatic SE -threads ${task.cpus} ${R1} \$base.trimmed.fastq.gz \
    ILLUMINACLIP:${ADAPTERS_SE}:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:${MINLEN}
    num_untrimmed=\$((\$(gunzip -c ${R1} | wc -l)/4))
    num_trimmed=\$((\$(gunzip -c \$base'.trimmed.fastq.gz' | wc -l)/4))
    percent_trimmed=\$((100-\$((100*num_trimmed/num_untrimmed))))
    echo Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Reference_Genome,Reference_Length,Mapped_Reads,Percent_Ref_Coverage,Min_Coverage,Mean_Coverage,Max_Coverage,Bam_Size,Consensus_Length,Percent_N,PCR CT,NCBI_Name,Method,Genbank_Submission > \$base'_summary.csv'
    printf "\$base,\$num_untrimmed,\$num_trimmed,\$percent_trimmed" >> \$base'_summary.csv'
    ls -latr
    """
} 
} else {
	process Trim_Reads_PE {
    // container "quay.io/biocontainers/trimmomatic:0.35--6"
    errorStrategy 'retry'
    maxRetries 3

   input:
        tuple val(base), file(R1), file(R2) from input_read_ch
        file ADAPTERS_PE
        val MINLEN
    output: 
        tuple env(base),file("*.trimmed.fastq.gz"),file("*summary.csv") into Trim_out_PE, Trim_out_PE_FQC
        // tuple val(base),file("${base}_results.csv") into Results_trimmed_ch

    publishDir "${params.outdir}fastq_trimmed", mode: 'copy',pattern:'*.trimmed.fastq*'
    
    script:
    """
    #!/bin/bash
    trimmomatic PE -threads ${task.cpus} ${R1} ${R2} ${base}.R1.paired.fastq.gz ${base}.R1.unpaired.fastq.gz ${base}.R2.paired.fastq.gz ${base}.R2.unpaired.fastq.gz \
	ILLUMINACLIP:${ADAPTERS_PE}:2:30:10:1:true LEADING:3 TRAILING:3 SLIDINGWINDOW:4:20 MINLEN:${MINLEN}
    num_untrimmed=\$((\$(gunzip -c ${R1} | wc -l)/4))
    num_trimmed=\$((\$(gunzip -c \$base'.trimmed.fastq.gz' | wc -l)/4))
    percent_trimmed=\$((100-\$((100*num_trimmed/num_untrimmed))))
    echo Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Reference_Genome,Reference_Length,Mapped_Reads,Percent_Ref_Coverage,Min_Coverage,Mean_Coverage,Max_Coverage,Bam_Size,Consensus_Length,Percent_N,PCR CT,NCBI_Name,Method,Genbank_Submission > \$base'_summary.csv'
    printf "\$base,\$num_untrimmed,\$num_trimmed,\$percent_trimmed" >> \$base'_summary.csv'
    ls -latr

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
        tuple val(base), file("${base}.trimmed.fastq.gz"),file("${base}_summary.csv") from Trim_out_SE
        file REFERENCE_FASTA

    output:
        tuple val(base), file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"), file("${base}_summary2.csv"), file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),env(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),env(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz") into Everything_ch
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
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${REFERENCE_FASTA} threads=8 covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map1_stats.txt 2>&1
    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${REFERENCE_FASTA} \$id > ${base}_mapped_ref_genome.fa
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map2.sam ref=${base}_mapped_ref_genome.fa threads=8 covstats=${base}_map2_bbmap_out.txt covhist=${base}_map2_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map2_stats.txt 2>&1
    head -n 1 ${base}_mapped_ref_genome.fa > ${base}_mapped_ref_genome_edited.fa
    grep -v ">" ${base}_mapped_ref_genome.fa | sed 's/U/T/g' >> ${base}_mapped_ref_genome_edited.fa
    mv ${base}_mapped_ref_genome_edited.fa ${base}_mapped_ref_genome.fa
    samtools faidx ${base}_mapped_ref_genome.fa
    awk 'NR == 2 || \$5 > max {number = \$3; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref_size_out.txt
    id_ref_size=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref_size_out.txt)
    echo \$id_ref_size >> ${base}_most_mapped_ref_size.txt
    reads_mapped=\$(cat ${base}_map2_stats.txt | grep "mapped:" | cut -d\$'\\t' -f3)
    cp ${base}_summary.csv ${base}_summary2.csv
    printf ",\$id" >> ${base}_summary2.csv
    printf ",\$id_ref_size" >> ${base}_summary2.csv
    printf ",\$reads_mapped" >> ${base}_summary2.csv
    printf ",\$ref_coverage" >> ${base}_summary2.csv
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
    tuple val(base), file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"), file("${base}_summary2.csv"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz") from Everything_ch
    output:
    tuple val(base), file("${base}.bam") into Aligned_bam_ch, Bam_ch
    tuple val(base), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),env(bamsize),file("${base}.sorted.bam.bai"),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"), file("${base}_summary.csv"),file("${base}.trimmed.fastq.gz") into Consensus_ch

    publishDir "${params.outdir}bam", mode: 'copy', pattern:'*.bam*'
    publishDir "${params.outdir}bam_sorted", mode: 'copy', pattern:'*.sorted.bam*'  
    publishDir "${params.outdir}txt_bam_flagstats", mode: 'copy', pattern:'*_flagstats.txt*'  

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
/*
 * Call variants & Generate Consensus
  */
 process Generate_Consensus {
    container "quay.io/greninger-lab/swift-pipeline:latest"
    errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(base), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),val(bamsize),file("${base}.sorted.bam.bai"),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"), file("${base}_summary.csv"),file("${base}.trimmed.fastq.gz") from Consensus_ch
    
    file VCFUTILS
    file SPLITCHR
    file TRIM_ENDS

    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}.consensus.fa"), file("${base}_summary.csv"), val(bamsize), val(id),file("${base}.trimmed.fastq.gz") into Consensus_Fasta_ch
    tuple val(base), file("${base}_pre_bcftools.vcf"), file("${base}_bcftools.vcf") into Consensus_Vcf_ch

    publishDir "${params.outdir}consensus", mode: 'copy', pattern:'*.consensus.fa*' 
    publishDir "${params.outdir}vcf", mode: 'copy', pattern:'*_bcftools.vcf*' 
    publishDir "${params.outdir}vcf_pre", mode: 'copy', pattern:'*_pre_bcftools.vcf*' 

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
process Mapping_final {
	errorStrategy 'retry'
    maxRetries 3

    input:
    tuple val(base), file("${base}_mapped_ref_genome.fa"), file("${base}.consensus.fa"), file("${base}_summary.csv"), val(bamsize), val(id),file("${base}.trimmed.fastq.gz") from Consensus_Fasta_ch
    
    output:
    tuple val(base),file("${base}_mapped_ref_genome.fa"), file("${base}.consensus-final.fa"), file("${base}.consensus.masked.fa"), file("${base}_map3.sam"), file("${base}_map3.bam"), file("${base}.map3.sorted.bam"), file("${base}.map3.sorted.bam.bai"), file("${base}_map3_stats.txt"), file("${base}.mpileup"), file("${base}_final_summary.csv"), val(bamsize), val(id) into Mapping_Final_ch

    publishDir "${params.outdir}mpileup_map3", mode: 'copy', pattern:'*.mpileup*'
    publishDir "${params.outdir}bam_map3", mode: 'copy', pattern:'*.map3.sorted.bam*'
    publishDir "${params.outdir}sam_map3", mode: 'copy', pattern:'*_map3.sam*'
    publishDir "${params.outdir}consensus-ivar", mode: 'copy', pattern:'*.consensus-final*'
    publishDir "${params.outdir}consensus-ivar-masked", mode: 'copy', pattern:'*.consensus.masked.fa*'
    publishDir "${params.outdir}txt_bbmap_map3_stats", mode: 'copy', pattern:'*_map3_stats.txt*'
    publishDir "${params.outdir}summary", mode: 'copy', pattern:'*.csv*'

    script:

    """
    #!/bin/bash
    
    
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus.fa threads=8 local=true interleaved=false maxindel=9 strictmaxindel > ${base}_map3_stats.txt 2>&1
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
    cat ${base}.mpileup | ivar consensus -q 15 -t 0.6 -m 3 -n N -p ${base}.consensus-final

    bedtools genomecov \\
        -bga \\
        -ibam ${base}.map3.sorted.bam \\
        -g ${base}_mapped_ref_genome.fa \\
        | awk '\$4 < 10' | bedtools merge > ${base}.mask.bed
    
    bedtools maskfasta \\
        -fi ${base}.consensus-final.fa \\
        -bed ${base}.mask.bed \\
        -fo ${base}.consensus.masked.fa

    sed -i 's/>.*/>${base}.ivar.masked.consensus/' ${base}.consensus.masked.fa
    sed -i 's/>.*/>${base}.ivar.consensus/' ${base}.consensus-final.fa

    awk '/^>/{if (l!="") print l; print; l=0; next}{l+=length(\$0)}END{print l}' ${base}.consensus-final.fa > bases.txt
    num_bases=\$(awk 'FNR==2{print val,\$1}' bases.txt)

    grep -v "^>" ${base}.consensus-final.fa | tr -cd N | wc -c > N.txt
    num_ns=\$(awk 'FNR==1{print val,\$1}' N.txt)
    echo "\$num_ns/\$num_bases*100" | bc -l > n_percent.txt
    percent_n=\$(awk 'FNR==1{print val,\$1}' n_percent.txt)
    printf ",\$num_bases" >> ${base}_summary.csv
    printf ",\$percent_n" >> ${base}_summary.csv

    cp ${base}_summary.csv ${base}_final_summary.csv
    
    """  
// CALL in summary directory to combine all summary.csv files for all samples into one summary with one header
// cd <out_directory>
// awk '(NR == 1) || (FNR > 1)' *.csv > hrv-pl_<run_name>_summary.csv
}
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
