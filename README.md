## HRV Pipeline v1.3

Github Repo:
https://github.com/greninger-lab/HRV_Genome_Mapping_Pipeline

Authors:
Paul RK Cruz <kurtisc@uw.edu>
Michelle J Lin <Mjlin@uw.edu>
Alex L Greninger <agrening@uw.edu>

 
UW Medicine | Virology

Department of Laboratory Medicine and Pathology

University of Washington

Created: April, 2021

Last Update: August 16, 2021

License: MIT

## Usage

    To run the pipeline, enter the following in the command line:
        nextflow run FILE_PATH/HRV_Genome_Mapping_Pipeline/main.nf --reads PATH_TO_FASTQ --outdir PATH_TO_OUTPUT_DIR
    Valid CLI Arguments:
    REQUIRED:
      --reads                       Path to input fastq.gz folder).
      --outdir                      The output directory where the results will be saved
    OPTIONAL:
      --withSampleSheet             Adds Sample Sheet information to Final Report Summary
      --ref_rv                      Overwrite set multi-fasta Rhinovirus reference file
      --ref_hcov                    Overwrite set multi-fasta Human Coronavirus reference file
      --ref_respp                   Overwrite set multi-fasta Influenza B reference file
      --ref_inflb                   Overwrite set multi-fasta Respiratory Panel reference file
      --ref_hpiv3                   Overwrite set multi-fasta HPIV3 reference file
	  --helpMsg		    Displays help message in terminal
      --singleEnd                   Specifies that the input fastq files are single end reads. Pipeline assumes paired-end reads if this flag is not specified.
	  --withFastQC		    Runs a quality control check on fastq files
      --skipTrimming                Skips the fastq trimmming process (not available in v1.3; will be released in v1.4)

### Description:
Human Respiratory Virus Pipeline was designed to run either single-end or paired end Illumina Next-Generation-Sequencing (NGS) sequence reads for Human respiratory virus discovery, analysis, and Genbank submission.

### PIPELINE OVERVIEW:
1. Trim Reads
    1.1. Trimmomatic - sequence read trimming of adaptors and low quality reads.
    
 2. Genome Mapping & alignment
    2.1. BBMap - align to MultiFasta Reference Virus Genome.
    2.2. Samtools - SAM and BAM file processing.

 4. Sort Bam
    4.1. Convert Sam to Bam
    4.2. Sort Bam file by coordinates
    4.3. Generate bam statistics
    
 5. Variant Calling
    5.1. Calculate the read coverage of positions in the genome
    5.2. Detect the single nucleotide polymorphisms (SNPs)
    5.3. Filter and report the SNP variants in variant calling format (VCF)
    
 6. Consensus Generation
    6.1. Generates an un-masked consensus.
   
 7. Final Mapping
    7.1. Perform final mapping to unmasked consensus and generates a masked consensus.
    
 8. Final Processing (optional; --withSampleSheet)
    8.1. Sample Sheet information is added to the Final Summary Report.
  
 9. Final Summary Report Generation
    9.1. Generates summary statistics for each sample in csv format.
       
 10. FastQC (optional)
    10.1. Sequence read quality control analysis. Output in HTML format.

Dependencies:

HRV-Docker includes all dependencies. Currently in v1.3, mapping step requires local dependencies.

    trimmomatic         conda install -c bioconda trimmomatic
    bbtools             conda install -c bioconda bbmap    
    bgzip               conda install -c bioconda tabix
    samtools            conda install -c bioconda samtools
    varscan             conda install -c bioconda varscan
    vcftools            conda install -c bioconda vcftools
    bcftools            conda install -c bioconda bcftools
    seqkit              conda install -c bioconda seqkit
    bedtools            conda install -c bioconda bedtools
    fastqc              conda install -c bioconda fastqc
    ivar                conda install -c bioconda ivar
    seqtk               conda install -c bioconda seqtk
    csvkit              conda install -c anaconda csvkit
    
PIPELINE SETUP

Setup Multifasta References:

1. Multifasta references containing Viral genome sequences formatted with accession numbers only; Rhinovirus, Human Coronavirus, Influenza B, HPIV3, Respiratory Virus Panel.

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
        
        nextflow run HRV_Genome_Mapping_Pipeline --helpMsg helpMsg

Run Pipeline on Single-end sequence reads:
        
        nextflow run HRV_Genome_Mapping_Pipeline --reads '/Users/example/' --outdir '/Users/example/example_output/' --singleEnd 

Run Pipeline on Paired-end sequence reads:
        
        nextflow run HRV_Genome_Mapping_Pipeline --reads '/Users/example/' --outdir '/Users/example/example_output/'
