#!/usr/bin/env nextflow
/*
 * STEP 1: Trimming
 * Trimming of low quality and short NGS sequences.
 */
process Trimming {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"
    errorStrategy 'retry'
    maxRetries 3

    input:
        file R1// from input_read_ch
        file ADAPTERS_SE
        val MINLEN
        val SETTING
        val LEADING
        val TRAILING
        val SWINDOW

    output:
        tuple env(base),file("*.trimmed.fastq.gz"), file("${R1}_num_trimmed.txt"),file("*summary.csv")// into Trim_out_SE, Trim_out_SE_FQC

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
/*
 * STEP 2: Aligning
 * Viral identification & mapping.
 */
process Aligning {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest" 
    // errorStrategy 'retry'
    // maxRetries 3
    // echo true

    input: 
        tuple val(base), file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_summary.csv")// from Trim_out_SE
        file Reference_rv
        file Reference_hcov
        file Reference_hpv
        file Reference_hpv_14        
        file Reference_inflb
        file Reference_hpiv3

    output:
        tuple val(base), file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"), file("${base}_summary2.csv"), file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),env(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),env(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt")// into Everything_ch
        tuple val(base), file("${base}_map1_histogram.txt"),file("${base}_map2_histogram.txt")// into BBmap_map1_hist_ch
        tuple val (base), file("*")// into Dump_map1_ch

    publishDir "${params.outdir}all_ref", mode: 'copy', pattern:'*all_ref.sam*'
    publishDir "${params.outdir}all_ref", mode: 'copy', pattern:'*all_ref*'    
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
    publishDir "${params.outdir}ref_ids", mode: 'copy', pattern:'*_all_ref_id.txt*'    
    publishDir "${params.outdir}ref_ids", mode: 'copy', pattern:'*_ids.txt*'

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


    # Rhinovirus`s
    if grep -q \$all_ref_id "${base}_rv_ids.txt";
    then
    echo "< Accession found in Rhinovirus multifasta file. hrv_ref_rhinovirus.fa will be used for mapping."

    # MAP 1
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${Reference_rv} threads=${task.cpus} covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map1_stats.txt 2>&1
    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${Reference_rv} \$id > ${base}_mapped_ref_genome.fa
    # MAP 2
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

    # MAP 1
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${Reference_hpv} threads=${task.cpus} covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map1_stats.txt 2>&1
    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${Reference_hpv} \$id > ${base}_mapped_ref_genome.fa
    # MAP 2    
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


    # Influenza B
    elif grep -q \$all_ref_id "${base}_inbflb_ids.txt";
    then
    echo "< Accession found in Influenza B multifasta file. hrv_ref_Influenza_b.fa will be used for mapping."

    # MAP 1
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${Reference_inflb} threads=${task.cpus} covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map1_stats.txt 2>&1
    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${Reference_inflb} \$id > ${base}_mapped_ref_genome.fa
    # MAP 2    
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

    # MAP 1
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${Reference_hcov} threads=${task.cpus} covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false maxindel=20 strictmaxindel -Xmx6g > ${base}_map1_stats.txt 2>&1
    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${Reference_hcov} \$id > ${base}_mapped_ref_genome.fa
    # MAP 2    
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

    # MAP 1
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${Reference_hpiv3} threads=${task.cpus} covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false maxindel=9 strictmaxindel -Xmx6g > ${base}_map1_stats.txt 2>&1
    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${Reference_hpiv3} \$id > ${base}_mapped_ref_genome.fa
    # MAP 2    
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
    
    # MAP 1    
    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map1.sam ref=${base}_all_ref.fa threads=${task.cpus} covstats=${base}_map1_bbmap_out.txt covhist=${base}_map1_histogram.txt local=true interleaved=false maxindel=20 strictmaxindel -Xmx6g > ${base}_map1_stats.txt 2>&1
    samtools view -S -b ${base}_map1.sam > ${base}_map1.bam
    samtools sort -@ 4 ${base}_map1.bam > ${base}.sorted.bam
    samtools index ${base}.sorted.bam
    samtools idxstats ${base}.sorted.bam > ${base}_idxstats.txt
    awk 'NR == 2 || \$5 > max {number = \$1; max = \$5} END {if (NR) print number, max}' < ${base}_map1_bbmap_out.txt > ${base}_most_mapped_ref.txt
    id=\$(awk 'FNR==1{print val,\$1}' ${base}_most_mapped_ref.txt)
    ref_coverage=\$(awk 'FNR==1{print val,\$2}' ${base}_most_mapped_ref.txt)
    samtools faidx ${base}_all_ref.fa \$id > ${base}_mapped_ref_genome.fa
    # MAP 2    
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
/*
 * STEP 2: Bam_Sorting
 * Sorting, indexing, and collecting of summary statistics from BAM files.
 */
process Bam_Sorting {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"
    // errorStrategy 'retry'        
	errorStrategy 'ignore'
    // maxRetries 3

    input: 
    tuple val(base), file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"), file("${base}_summary2.csv"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt")// from Everything_ch
    output:
    tuple val(base), file("${base}.bam"), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),env(bamsize),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"), file("${base}_summary.csv"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt")// into Consensus_ch

    publishDir "${params.outdir}bam_map2", mode: 'copy', pattern:'*.sorted.bam*'  
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
    cp ${base}_summary2.csv ${base}_summary3.csv
    printf ",\$mincoverage" >> ${base}_summary3.csv
    printf ",\$meancoverage" >> ${base}_summary3.csv
    printf ",\$maxcoverage" >> ${base}_summary3.csv
    printf ",\$bamsize" >> ${base}_summary3.csv
    cp ${base}_summary3.csv ${base}_summary.csv
    """
}
/*
 * STEP 3: Generate_Consensus
 * Consensus fasta generation.
 */
process Consensus_Generation {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"     
    // errorStrategy 'retry'
	// errorStrategy 'ignore'
    // maxRetries 3
    // echo true

    input:
    tuple val(base), file("${base}.bam"), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),val(bamsize),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"), file("${base}_summary.csv"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt")// from Consensus_ch
    
    output:
    tuple val(base), file("${base}.mpileup"), file("${base}.bam"), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),val(bamsize),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"), file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt"), file("${base}_final_summary.csv"),file("${base}.consensus_final.fa")// into Consensus_Fasta_ch

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
        --fasta-ref ${base}_mapped_ref_genome.fa \\
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
        --fasta-ref ${base}_mapped_ref_genome.fa \\
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
        --fasta-ref ${base}_mapped_ref_genome.fa \\
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
        --fasta-ref ${base}_mapped_ref_genome.fa \\
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
        --fasta-ref ${base}_mapped_ref_genome.fa \\
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
        --fasta-ref ${base}_mapped_ref_genome.fa \\
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
/*
 * STEP 4: Aligning_Final
 * Final round of mapping.
 */
process Aligning_Final {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"     
    // errorStrategy 'retry'	
    errorStrategy 'ignore'
    // maxRetries 3
    // echo true

    input:
    tuple val(base), file("${base}.mpileup"), file("${base}.bam"), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),val(bamsize),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"), file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt"), file("${base}_final_summary.csv"),file("${base}.consensus_final.fa")// from Consensus_Fasta_ch

    output:
    tuple val(base), file("${base}.mpileup"), file("${base}.bam"), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),val(bamsize),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt"), file("${base}_summary.csv"),file("${base}.consensus_final.fa"),file("${base}_map3.sam"),file("${base}_map3.sorted.bam")// into Mapping_Final_ch, Final_Processing_ch

    publishDir "${params.outdir}bam_map3", mode: 'copy', pattern:'*_map3.sorted.bam*'
    publishDir "${params.outdir}sam_map3", mode: 'copy', pattern:'*_map3.sam*'

    script:
    """
    #!/bin/bash
    # FINAL MAPPING - MAP3

    all_ref_id=\$(awk '{print \$1}' ${base}_all_ref_id.txt)

    # Rhinovirus
    if grep -q \$all_ref_id "${base}_rv_ids.txt"; 
    then
    echo "< Accession found in Rhinovirus multifasta file. hrv_ref_rhinovirus.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats_map3.txt 2>&1

    samtools view -S -b ${base}_map3.sam > ${base}_map3.bam
    samtools sort -@ 4 ${base}_map3.bam > ${base}_map3.sorted.bam
    samtools index ${base}_map3.sorted.bam
    # Rename summary file to retain file caching/-resume functionality
    cp ${base}_final_summary.csv ${base}_summary.csv


    # HPV
    elif grep -q \$all_ref_id "${base}_hpv_ids.txt";
    then
    echo "< Accession found in HPV multifasta file. hrv_ref_hpv.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats_map3.txt 2>&1

    samtools view -S -b ${base}_map3.sam > ${base}_map3.bam
    samtools sort -@ 4 ${base}_map3.bam > ${base}_map3.sorted.bam
    samtools index ${base}_map3.sorted.bam
    # Rename summary file to retain file caching/-resume functionality
    cp ${base}_final_summary.csv ${base}_summary.csv


    # Influenza B
    elif grep -q \$all_ref_id "${base}_inbflb_ids.txt";
    then
    echo "< Accession found in Influenza B multifasta file. hrv_ref_Influenza_b.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats_map3.txt 2>&1

    samtools view -S -b ${base}_map3.sam > ${base}_map3.bam
    samtools sort -@ 4 ${base}_map3.bam > ${base}_map3.sorted.bam
    samtools index ${base}_map3.sorted.bam
    # Rename summary file to retain file caching/-resume functionality
    cp ${base}_final_summary.csv ${base}_summary.csv


    # Human Coronavirus
    elif grep -q \$all_ref_id "${base}_hcov_ids.txt";
    then
    echo "Accession found in HCoVs multifasta file. hrv_ref_hcov.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats_map3.txt 2>&1

    samtools view -S -b ${base}_map3.sam > ${base}_map3.bam
    samtools sort -@ 4 ${base}_map3.bam > ${base}_map3.sorted.bam
    samtools index ${base}_map3.sorted.bam
    # Rename summary file to retain file caching/-resume functionality
    cp ${base}_final_summary.csv ${base}_summary.csv


    # HPIV3 - Human parainfluenza virus 3
    elif grep -q \$all_ref_id "${base}_hpiv3.txt";
    then
    echo "Accession found in HPIV3 multifasta file. hrv_ref_hpiv3.fa will be used for mapping."

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats_map3.txt 2>&1

    samtools view -S -b ${base}_map3.sam > ${base}_map3.bam
    samtools sort -@ 4 ${base}_map3.bam > ${base}_map3.sorted.bam
    samtools index ${base}_map3.sorted.bam
    # Rename summary file to retain file caching/-resume functionality
    cp ${base}_final_summary.csv ${base}_summary.csv

    else

    ${BBMAP_PATH}bbmap.sh in=${base}.trimmed.fastq.gz outm=${base}_map3.sam ref=${base}.consensus_final.fa threads=${task.cpus} local=true interleaved=false maxindel=9 -Xmx6g > ${base}_final_mapping_stats_map3.txt 2>&1

    samtools view -S -b ${base}_map3.sam > ${base}_map3.bam
    samtools sort -@ 4 ${base}_map3.bam > ${base}_map3.sorted.bam
    samtools index ${base}_map3.sorted.bam
    # Rename summary file to retain file caching/-resume functionality
    cp ${base}_final_summary.csv ${base}_summary.csv

    fi

    """  
}
/*
 * OPTIONAL: withMetadata - Summary_Generation
 * Integration of metadata to final run summary statistics.
 */
process Summary_Generation {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"
    // errorStrategy 'retry'        
    errorStrategy 'ignore'
    // maxRetries 3
    // echo true

    input:
    file SAMPLE_LIST from METADATA

    tuple val(base), file("${base}.mpileup"), file("${base}.bam"), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),val(bamsize),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt"), file("${base}_summary.csv"),file("${base}.consensus_final.fa"),file("${base}_map3.sam"),file("${base}_map3.sorted.bam")// from Mapping_Final_ch

    output:
    tuple val(base), file("${base}.mpileup"), file("${base}.bam"), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),val(bamsize),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_sample_id.txt"), file("${base}_pcr_ct.txt"), file("${base}_method.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt"), file("${base}_final_summary.csv"), file("${base}.consensus_final.fa"), file("${base}_map3.sam"), file("${base}_map3.sorted.bam"), file("${base}_sample_stats.csv")// into Final_Processing_final_ch

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
/*
 * OPTIONAL: withSerotype - Serotyping
 * Viral Serotyping/Genotyping
 */
process Serotyping {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"        
    // errorStrategy 'retry'  
    // maxRetries 3
    // echo true

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

    tuple val(base), file("${base}.mpileup"), file("${base}.bam"), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),val(bamsize),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_sample_id.txt"), file("${base}_pcr_ct.txt"), file("${base}_method.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt"), file("${base}_final_summary.csv"),file("${base}.consensus_final.fa"),file("${base}_map3.sam"),file("${base}_map3.sorted.bam"),file("${base}_sample_stats.csv")// from Final_Processing_final_ch  

    output:
    tuple val(base), file("${base}.mpileup"), file("${base}.bam"), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),val(bamsize),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_sample_id.txt"), file("${base}_pcr_ct.txt"), file("${base}_method.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt"), file("${base}.consensus_final.fa"),file("${base}_map3.sam"),file("${base}_map3.sorted.bam"),file("${base}_sample_stats.csv"), file("${base}_collection_year.txt"), file("${base}_country_collected.txt"), file("${base}_blast_db_vp1.txt"), file("${base}_blast_db_all_ref.txt"), file("${base}_sample_stats.csv"), file("${base}_all_ref_id.txt"), file("${base}_nomen.txt")// into Serotype_ch
    tuple val(base), file("${base}_summary_final.csv")// into Summary_cat_ch

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
    SAMPLEName=\${R1:2:14}
    char_edit=\$(echo "\${R1//_R1}")
    all_ref_id=\$(awk '{print \$1}' ${base}_all_ref_id.txt)

    # Rhinovirus
    if grep -q \$all_ref_id "${base}_rv_ids.txt";
    then
    echo "< Accession found in Rhinovirus multifasta file; RV typing."

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
    echo "< Accession found in HPV multifasta file; HPV typing."
    
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

    blastn -out ${base}_blast_db_vp1.txt -query ${base}.consensus_final.fa -db ${hpv_db_1} -outfmt 6 -task blastn -max_target_seqs 1 -evalue 1e-5

    serotype=\$(awk 'FNR==1{print val,\$2}' ${base}_blast_db_vp1.txt)
    cut -d "-" -f2- <<< "\$serotype" > ${base}_serotype-parse.txt
    serotype_parsed=\$(awk 'FNR==1{print val,\$1}' ${base}_serotype-parse.txt)
    rv='Human Papilloma Virus '
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


    # Human Coronavirus
    elif grep -q \$all_ref_id "${base}_hcov_ids.txt";
    then
    echo "Accession found in HCoVs multifasta file. hrv_ref_hcov.fa will be used for mapping."

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
    rv='HCoV_'
    serotype_parse="\${rv} \${serotype_parsed}" 
    echo \$serotype_parse > ${base}_sero.txt
    cat ${base}_sero.txt | tr -d " \t\n\r" > ${base}_serot.txt
    serotype_parsed2=\$(awk 'FNR==1{print val,\$1}' ${base}_serot.txt)
    space='/'
    city='Seattle'

    nomenclature="\${rv}\${serotype_parsed2} \${space} \${country_collected} \${space} \${city} \${space} \${NCBI_Name} \${space} \${collection_year}" 
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
    

    # Influenza B
    elif grep -q \$all_ref_id "${base}_inbflb_ids.txt";
    then
    echo "< Accession found in Influenza B multifasta file; IFB typing."


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
    echo "Accession found in HPIV3 multifasta file; HPIV3 typing."


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

    echo "Accession NOT found; using ALL_REF for typing"

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
    /*
 * STEP 5: Final_Processing
 * Final data summary statistic data processing.
 */
process Final_Processing {
    // container "docker.io/paulrkcruz/hrv-pipeline:latest"        
    // errorStrategy 'retry'
    errorStrategy 'ignore'
    // maxRetries 3
    // echo true

    input:
    file '*.csv' from Summary_cat_ch.collect()
    tuple val(base), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),val(bamsize),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt"), file("${base}_summary.csv"),file("${base}.consensus_final.fa"),file("${base}_map3.sam"),file("${base}_map3.sorted.bam")// from Final_Processing_ch

    output:
    file("Run_Summary_Final_cat.csv") into Final_summary_out_ch
    tuple val(base), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),val(bamsize),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt"), file("${base}_summary.csv"),file("${base}.consensus_final.fa"),file("${base}_map3.sam"),file("${base}_map3.sorted.bam")// into Final_Processing_out_ch

    publishDir "${params.outdir}summary_withserotype_cat", mode: 'copy', pattern:'*Run_Summary_Final_cat.csv*'

    script:
    """
    #!/bin/bash
    cp ${params.outdir}/ref_ids/${base}_all_ref_id.txt ${base}_all_ref_id.txt
    cp ${params.outdir}/ref_ids/${base}_rv_ids.txt ${base}_rv_ids.txt
    cp ${params.outdir}/ref_ids/${base}_hpv_ids.txt ${base}_hpv_ids.txt
    cp ${params.outdir}/ref_ids/${base}_inbflb_ids.txt ${base}_inbflb_ids.txt
    cp ${params.outdir}/ref_ids/${base}_hcov_ids.txt ${base}_hcov_ids.txt
    cp ${params.outdir}/ref_ids/${base}_hpiv3.txt ${base}_hpiv3.txt
    R1=${base}
    NCBI_Name=\${R1:4:6}
    SAMPLEName=\${R1:2:5}
    all_ref_id=\$(awk '{print \$1}' ${base}_all_ref_id.txt)
    # Rhinovirus
    if grep -q \$all_ref_id "${base}_rv_ids.txt"; 
    then
    echo "< Accession found in Rhinovirus multifasta file. hrv_ref_rhinovirus.fa will be used for final processing."
    awk '(NR == 1) || (FNR > 1)' *.csv >  Run_Summary_cat.csv
    sed '1d' Run_Summary_cat.csv > Run_Summary_catted.csv
    echo -e "Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Reference_Genome,Reference_Length,Mapped_Reads,Percent_Ref_Coverage,Min_Coverage,Mean_Coverage,Max_Coverage,Bam_Size,Consensus_Length,Percent_N,%_Reads_On_Target, PCR_CT,Method, NCBI_Name, Serotype, Nomenclature, Reference_Name, Reference_Genome, Biosample_name, Biosample_accession, SRA_Accession, Bioproject,  Release_date" | cat - Run_Summary_catted.csv > Run_Summary_Final_cat.csv
    # HPV
    elif grep -q \$all_ref_id "${base}_hpv_ids.txt";
    then
    echo "< Accession found in HPV multifasta file. HPV final processing."
    awk '(NR == 1) || (FNR > 1)' *.csv >  Run_Summary_cat.csv
    sed '1d' Run_Summary_cat.csv > Run_Summary_catted.csv
    echo -e "Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Reference_Genome,Reference_Length,Mapped_Reads,Percent_Ref_Coverage,Min_Coverage,Mean_Coverage,Max_Coverage,Bam_Size,Consensus_Length,Percent_N,%_Reads_On_Target, PCR_CT,Method, NCBI_Name, Serotype, Nomenclature, Reference_Name, Reference_Genome, Biosample_name, Biosample_accession, SRA_Accession, Bioproject,  Release_date" | cat - Run_Summary_catted.csv > Run_Summary_Final_cat.csv
    # Influenza B
    elif grep -q \$all_ref_id "${base}_inbflb_ids.txt";
    then
    echo "< Accession found in Influenza B multifasta file. hrv_ref_Influenza_b.fa will be used for final processing."
    awk '(NR == 1) || (FNR > 1)' *.csv >  Run_Summary_cat.csv
    sed '1d' Run_Summary_cat.csv > Run_Summary_catted.csv
    echo -e "Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Reference_Genome,Reference_Length,Mapped_Reads,Percent_Ref_Coverage,Min_Coverage,Mean_Coverage,Max_Coverage,Bam_Size,Consensus_Length,Percent_N,%_Reads_On_Target, PCR_CT,Method, NCBI_Name, Serotype, Nomenclature, Reference_Name, Reference_Genome, Biosample_name, Biosample_accession, SRA_Accession, Release_date, Bioproject" | cat - Run_Summary_catted.csv > Run_Summary_Final_cat.csv
    # Human Coronavirus
    elif grep -q \$all_ref_id "${base}_hcov_ids.txt";
    then
    echo "Accession found in HCoVs multifasta file. hrv_ref_hcov.fa will be used for final processing."
    awk '(NR == 1) || (FNR > 1)' *.csv >  Run_Summary_cat.csv
    sed '1d' Run_Summary_cat.csv > Run_Summary_catted.csv
    echo -e "Sample_Name,Raw_Reads,Trimmed_Reads,Percent_Trimmed,Reference_Genome,Reference_Length,Mapped_Reads,Percent_Ref_Coverage,Min_Coverage,Mean_Coverage,Max_Coverage,Bam_Size,Consensus_Length,Percent_N,%_Reads_On_Target, PCR_CT,Method, NCBI_Name, Serotype, Nomenclature, Reference_Name, Reference_Genome, Biosample_name, Biosample_accession, SRA_Accession, Release_date, Bioproject" | cat - Run_Summary_catted.csv > Run_Summary_Final_cat.csv
    # HPIV3 - Human parainfluenza virus 3
    elif grep -q \$all_ref_id "${base}_hpiv3.txt";
    then
    echo "Accession found in HPIV3 multifasta file. hrv_ref_hpiv3.fa will be used for final processing."
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
/*
 * OPTIONAL: withVapid - Viral_Annotation
 * Viral annotation and creation of *.sqn GenBank submission files.
 */
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

    tuple val(base), file("${base}.mpileup"), file("${base}.bam"), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),val(bamsize),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_sample_id.txt"), file("${base}_pcr_ct.txt"), file("${base}_method.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt"), file("${base}.consensus_final.fa"),file("${base}_map3.sam"),file("${base}_map3.sorted.bam"),file("${base}_sample_stats.csv"), file("${base}_collection_year.txt"), file("${base}_country_collected.txt"), file("${base}_blast_db_vp1.txt"), file("${base}_blast_db_all_ref.txt"), file("${base}_sample_stats.csv"), file("${base}_all_ref_id.txt"), file("${base}_nomen.txt")// from Serotype_ch
    output:
    tuple val(base), file("${base}.mpileup"), file("${base}.bam"), file("${base}.sorted.bam"),file("${base}_flagstats.txt"),val(bamsize),file("${base}_map2.sam"), file("${base}_most_mapped_ref.txt"),file("${base}_most_mapped_ref_size.txt"),file("${base}_most_mapped_ref_size_out.txt"),val(id_ref_size),file("${base}_idxstats.txt"),file("${base}_mapped_ref_genome.fa"),val(id),file("${base}_map1_bbmap_out.txt"),file("${base}_map2_bbmap_out.txt"),file("${base}_map1_stats.txt"),file("${base}_map2_stats.txt"),file("${base}_mapped_ref_genome.fa.fai"),file("${base}.trimmed.fastq.gz"), file("${base}_num_trimmed.txt"), file("${base}_num_mapped.txt"), file("${base}_sample_id.txt"), file("${base}_pcr_ct.txt"), file("${base}_method.txt"), file("${base}_rv_ids.txt"), file("${base}_hpv_ids.txt"), file("${base}_inbflb_ids.txt"), file("${base}_hcov_ids.txt"), file("${base}_hpiv3.txt"), file("${base}_all_ref_id.txt"), file("${base}.consensus_final.fa"),file("${base}_map3.sam"),file("${base}_map3.sorted.bam"),file("${base}_sample_stats.csv"), file("${base}_collection_year.txt"), file("${base}_country_collected.txt"), file("${base}_blast_db_vp1.txt"), file("${base}_blast_db_all_ref.txt"), file("${base}_nomen.txt"), file ("${base}_vapid_metadata.csv")// into All_files_ch

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
    sed "100s/.*/"\$sub"/" ${base}_95.txt > ${base}_100.txt
    
    # Edit Line:188 - str - replace nomenclature with ncbi_name
    subname_1_1="\$quote\$nomen\$quote"
    subname_2_2="                   subname \$subname_1_1 } ,"
    sub=\$subname_2_2
    sed "100s/.*/"\$sub"/" ${base}_95.txt > ${base}_100.txt




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
/*
 * OPTIONAL: withFastQC - FastQC
 * Evaluation of fastq reads. Ooutputs *.html files with fastq quality statistics.
 */    
process FastQC {
    container "docker.io/paulrkcruz/hrv-pipeline:latest"
	errorStrategy 'retry'
    maxRetries 3

    input:
        file R1// from Trim_out_fastqc_SE

    output:
	file '*_fastqc.{zip,html}'// into Fastqc_results_SE

    publishDir "${params.outdir}fastqc_results", mode: 'copy', pattern:'*_fastqc.{zip,html}*'  

    script:
    """
    #!/bin/bash

    /usr/local/bin/fastqc --quiet --threads ${task.cpus} ${base}.trimmed.fastq.gz
    """
}