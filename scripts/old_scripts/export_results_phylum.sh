#!/bin/bash

qiime taxa collapse \
--i-table table.qza --i-taxonomy taxonomy.qza \
--p-level 2 \
--o-collapsed-table table_collapsed_phylum.qza

mkdir phylum_classification

mv table_collapsed_phylum.qza ./phylum_classification

cd ./phylum_classification

qiime tools export \
--input-path table_collapsed_phylum.qza \
--output-path .

mv feature-table.biom feature-table_absfreq.biom

biom convert \
-i feature-table_absfreq.biom \
-o feature-table_absfreq.tsv \
--to-tsv \
--table-type 'Taxon table'

qiime feature-table relative-frequency \
--i-table table_collapsed_phylum.qza \
--o-relative-frequency-table table_collapsed_relfreq.qza

qiime metadata tabulate  \
--m-input-file table_collapsed_relfreq.qza  \
--o-visualization table_collapsed_relfreq.qzv

qiime tools export \
--input-path table_collapsed_relfreq.qza \
--output-path .

mv feature-table.biom feature-table_relfreq.biom

biom convert \
-i feature-table_relfreq.biom \
-o feature-table_relfreq.tsv \
--to-tsv \
--table-type 'Taxon table'
