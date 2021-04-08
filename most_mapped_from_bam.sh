    #!/bin/bash
    # USE UNSORTED BAM FILE TO READ REF GENOME OF MOST MAPPED READS THEN SAVE INFO TO TEXT FILE
    bedtools bamtobed -i V338470881_S9_L001_R1_001.bam | head -1 > most_mapped.txt
    # Extract Accession id from text file
    filename='most_mapped.txt'
    while read -r line
    do
    id=$(cut -c-8 <<< "$line")
    echo $id
    #code for passing id to other script file as parameter
    done < "$filename"
    # USE SAM TOOLS TO EXTRACT GENOME REFERENCE FROM MULTI-FASTA - SAVE TO NEW FASTA FOR LATER USE
    samtools faidx rhv_ref_db01.fasta ${id} > MAPPED-REF-GENOME.fasta