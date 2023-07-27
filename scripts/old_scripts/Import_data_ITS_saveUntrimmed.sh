#!/bin/bash

#
# Copyright 2021 Simone Maestri. All rights reserved.
# Simone Maestri <simone.maestri@univr.it>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.
#

MANIFEST=$1
FW_PRIMER=$2
RV_PRIMER=$3


qiime tools import \
--type 'SampleData[PairedEndSequencesWithQuality]' \
--input-path $MANIFEST \
--output-path sequences_untrimmed.qza \
--input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
--i-data sequences_untrimmed.qza \
--o-visualization demux_summary_untrimmed.qzv

#Note that we trim the forward primer and the reverse complement of the reverse primer from the forward reads
# (the forward primers have already been trimmed in the raw reads, but we will demonstrate forward + 
#reverse trimming here since attempting to trim the forward read will not hurt). 
#We trim the reverse primer and reverse complement of the forward primer from the reverse reads.

qiime cutadapt trim-paired \
  --i-demultiplexed-sequences sequences_untrimmed.qza \
  --p-cores 20 \
  --p-front-f $FW_PRIMER \
  --p-front-r $RV_PRIMER \
  --verbose \
  --o-trimmed-sequences sequences.qza


qiime demux summarize \
--i-data sequences.qza \
--o-visualization demux_summary.qzv

qiime tools export \
  --input-path sequences.qza \
  --output-path sequence_after_adapterTrimming
