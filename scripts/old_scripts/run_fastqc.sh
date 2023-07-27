#!/bin/bash

###fastqc on raw reads
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
	mkdir fastqc;
        fastqc --outdir fastqc --extract -t 15 *R1_001.fastq.gz *R2_001.fastq.gz;
        echo 'analysis finished for sample' $i;
done
