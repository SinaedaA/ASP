#!/bin/bash

MANIFEST=$1


qiime tools import \
--type 'SampleData[SequencesWithQuality]' \
--input-path $MANIFEST \
--output-path fj-joined-demux.qza \
--input-format SingleEndFastqManifestPhred33V2

qiime demux summarize \
--i-data fj-joined-demux.qza \
--o-visualization fj-joined-demux.qzv
