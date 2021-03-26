#!/bin/bash -ue
trimmomatic PE -threads 1 -phred33 test_fastq test_fastq"_paired_R1.fastq" test_fastq"_unpaired_R1.fastq" test_fastq"_paired_R2.fastq" test_fastq"_unpaired_R2.fastq" ILLUMINACLIP:/Users/kurtiscruz/opt/anaconda3/pkgs/trimmomatic-0.39-0/share/trimmomatic-0.39-0/adapters/TruSeq3-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:101 2> test_fastq.log
