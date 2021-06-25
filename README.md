## HRV Pipeline v1.2

Github Repo:
https://github.com/greninger-lab/HRV_Genome_Mapping_Pipeline

Author:
Paul RK Cruz <kurtisc@uw.edu>
Michelle J Lin <Mjlin@uw.edu>
Alex L Greninger <agrening@uw.edu>

 
UW Medicine | Virology

Department of Laboratory Medicine and Pathology

University of Washington

Created: April, 2021

License: MIT

### Description:
Human Respiratory Virus (HRV) Pipeline was designed to run either single-end or paired end Illumina Next-Generation-Sequencing (NGS) sequence reads to identify Viral complete genomes for analysis and Genbank submission.

### PIPELINE OVERVIEW:
1. Trim Reads

    1.1. Trimmomatic - sequence read trimming of adaptors and low quality reads.
    
 2. Genome Mapping
 
 	2.1. BBMap - align to MultiFasta Reference Virus Genome.
 	
 	2.2. Samtools - SAM and BAM file processing.
 	
 3. Reference Fasta Generation
 
 	3.1. Generate a fasta reference from the genome mapping results.
 	
 4. Sort Bam
 
    4.1. Convert Sam to Bam
    
    4.2. Sort Bam file by coordinates
    
    4.3. Generate Statistics about Bam file
    
 5. Variant Calling
 
    5.1. Calculate the read coverage of positions in the genome
    
    5.2. Detect the single nucleotide polymorphisms (SNPs)
    
    5.3. Filter and report the SNP variants in variant calling format (VCF)
    
    5.4. CLI Command to view results:   less -S ${base}_final_variants.vcf
    
 6. Consensus
    6.1. Consensus generation using variants VCF, mapped reference fasta, and
    sorted bam.
   
 7. Final Consensus
    7.1. Creates the Final Consensus by editing the fasta header.
    
 8. FastQC
 	8.1. Sequence read quality control analysis.

Dependencies:
                        
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

Setup Multifasta Reference:

1. REFERENCE_FASTA (must be a multifasta/fasta containing  full length Viral genome sequences (6-10K bp) formatted with accession numbers only)

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
        
        nextflow run /Users/Kurtisc/Downloads/CURRENT/HRV_Genome_Mapping_Pipeline --helpMsg helpMsg

Run Pipeline on Single-end sequence reads:
        
        nextflow run /Users/Kurtisc/Downloads/CURRENT/HRV_Genome_Mapping_Pipeline --reads '/Users/example/' --outdir '/Users/example/example_output/'  --Reference_Fasta /Users/example/hrv_ref/hrv_ref_db01_accession_only.fa --singleEnd 

Run Pipeline on Paired-end sequence reads:
        
        nextflow run /Users/Kurtisc/Downloads/CURRENT/HRV_Genome_Mapping_Pipeline --reads '/Users/example/' --outdir '/Users/example/example_output/'  --Reference_Fasta /Users/example/hrv_ref/hrv_ref_db01_accession_only.fa
        
