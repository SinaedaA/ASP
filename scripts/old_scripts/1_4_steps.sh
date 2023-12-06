#!/bin/bash

mkdir run1; cd run1
mkdir logfiles

READ_DIR=$1
TIMESTAMP=$(date '+%Y%m%d-%H%M%S')

../../scripts/0_make_dir_branches.sh $READ_DIR 2>&1 | tee logfiles/0_make_dir_branches_${TIMESTAMP}.log
../../scripts/1_fastqc.sh ../sample-metadata.tsv 0_raw_reads/ --outpath 1_fastqc/ 2>&1 | tee logfiles/1_fastqc_${TIMESTAMP}.log
../../scripts/2_make_manifest.sh ../sample-metadata.tsv 0_raw_reads/ 2>&1 | tee logfiles/2_make_manifest_${TIMESTAMP}.log
../../scripts/3_cutadapt.sh ../primer_file.txt 0_raw_reads/ 2>&1 | tee logfiles/3_cutadapt_${TIMESTAMP}.log
Rscript ../../scripts/4_dada2.R 3_analysis/3.1_cutadapt/ 3_analysis/3.2_trimming/dada2/ --join dada2 --primerfile ../primer_file.txt --truncQ 2 2>&1 | tee logfiles/4_dada2_${TIMESTAMP}.log
Rscript ../../scripts/4_dada2.R 3_analysis/3.1_cutadapt/ 3_analysis/3.2_trimming/dada2_flash2/ --join flash2 --primerfile ../primer_file.txt --truncQ 2 2>&1 | tee logfiles/4_flash2_dada2_${TIMESTAMP}.log