#!/usr/bin/bash

mkdir -p target/

# Filter SNPs that are in ATAC peaks.
sh bin/snps_in_peaks.sh

# Merge the peaks and output population-specific statistics.
python bin/peak_merge.py -p 0.75 \
       -i data/peaks/all_peaks_sorted.bed.gz \
       > target/merged_peaks.bed
