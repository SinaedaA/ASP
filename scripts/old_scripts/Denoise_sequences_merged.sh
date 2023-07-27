#!/bin/bash

TRUNC_LEN_F=$1
TRIM_LEFT_F=$2

qiime dada2 denoise-single \
--i-demultiplexed-seqs fj-joined-demux.qza \
--p-trunc-len $TRUNC_LEN_F --p-trim-left $TRIM_LEFT_F \
--p-trunc-q 0 \
--o-table table.qza --o-representative-sequences rep-seqs.qza --o-denoising-stats denoising-stats.qza \
--p-n-threads 0 \
--verbose

qiime feature-table summarize \
  --i-table table.qza \
  --o-visualization table.qzv \
  --m-sample-metadata-file sample-metadata.tsv

qiime feature-table tabulate-seqs \
  --i-data rep-seqs.qza \
  --o-visualization rep-seqs.qzv

qiime metadata tabulate \
  --m-input-file denoising-stats.qza \
  --o-visualization denoising-stats.qzv
