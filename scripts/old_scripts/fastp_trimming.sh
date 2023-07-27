#!/bin/bash

###trimming reads
input_folder=$1;
input_samples=$2;

echo '#############';
echo 'main working directory is ' $input_folder;
echo 'working on samples included in the file ' $input_samples;

cd $input_folder;
cat $input_folder/$input_samples;

for i in $(cat $input_folder/$input_samples);
do
	echo $i;
	cd $input_folder/$i;
	echo 'working directory';
	pwd;
	fastp -i *R1_001.fastq.gz -I *R2_001.fastq.gz -o $(basename *R1_001.fastq.gz | cut -d "." -f 1).trimmed.R1.fastq.gz -O $(basename *R2_001.fastq.gz | cut -d "." -f 1).trimmed.R2.fastq.gz --thread 7;
	echo 'analysis finished for sample' $i;
done
