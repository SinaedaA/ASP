#!/bin/bash

# Script that concatenates sequencing runs when there are several.
# Written by Stefania Concetta Quattro
# Modified by Sinaeda Anderssen
# 2 positional arguments (run files) as Infinite positional (could be more than 2, but min is 2)
# 1 positional argument = path to sample_names
# 1 optional argument = output_dir (otherwise, output in current dir, so default = PWD)

###fastqc on raw reads
run1_directory=$1
run2_directory=$2;
output_dir=$3;
sample_names=$4;

echo '#############';
echo 'Directory containing run 1 is ' $run1_directory;
echo 'Directory containing run 2 is ' $run2_directory;
echo 'Output directory of merged libraries is ' $output_dir;
echo 'Working on samples included in the file ' $sample_names;

cd $output_dir;

for i in $(cat $sample_names);
do
	echo $i;
	mkdir $i && cd "$_";
	echo "Working directory $PWD";
	echo $run1_directory/$i-*/*_R1_001.fastq.gz;
	echo $run1_directory/$i-*/*_R2_001.fastq.gz;
	echo $run2_directory/$i-*/*_R1_001.fastq.gz;
	echo $run2_directory/$i-*/*_R2_001.fastq.gz;
	cat $run1_directory/$i-*/*_R1_001.fastq.gz $run2_directory/$i-*/*_R1_001.fastq.gz > ${i}_R1_001.fastq.gz
	cat $run1_directory/$i-*/*_R2_001.fastq.gz $run2_directory/$i-*/*_R2_001.fastq.gz > ${i}_R2_001.fastq.gz
	mkdir fastqc;
	fastqc --outdir fastqc --extract -t 15 *R1_001.fastq.gz *R2_001.fastq.gz;
	cd -;
	echo 'Quality analysis (FASTQC) finished for sample' $i;
done
