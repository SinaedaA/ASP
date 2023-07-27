#!/bin/bash

###fastqc on raw reads
input_folder=$1;
input_samples=$2;
min_overlap=$3;

echo '#############';
echo 'main working directory is ' $input_folder;
echo 'working on samples included in the file ' $input_samples;

cd $input_folder;
cat $input_folder/$input_samples;

for i in $(cat $input_folder/$input_samples);
do
	cd $input_folder;
	echo $i;
	mkdir $i;
	pwd;
	R1=$(realpath $(find $input_folder | grep $i"_" | grep "R1" | grep "\\.fastq\\.gz"));
	R2=$(realpath $(find $input_folder | grep $i"_" | grep "R2" | grep "\\.fastq\\.gz"));
	mv $R1 $R2 ./$i;
	pwd;
	cd $i;
	flash2 --min-overlap=$min_overlap --max-overlap=300 --compress $i*R1_001.fastq.gz $i*R2_001.fastq.gz --output-prefix=$i 2>&1 | tee flash.log;
	echo 'analysis finished for sample' $i;
done
