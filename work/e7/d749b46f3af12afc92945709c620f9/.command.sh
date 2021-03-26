#!/bin/bash -ue
mkdir tmp
fastqc -t 1 -dir tmp test_fastq
rm -rf tmp
