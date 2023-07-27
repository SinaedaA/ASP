#!/bin/bash

SAMPLE_METADATA=$1
READS_DIR=$2

echo -e sample-id"\t"absolute-filepath > manifest.txt

for s in $(cat $SAMPLE_METADATA | cut -f1 ); do
  echo $s;
  R1=$(realpath $(find $READS_DIR | grep $s"." | grep "extendedFrags" | grep "\\.fastq\\.gz"));
  echo -e $s"\t"$R1 >> manifest.txt
done
