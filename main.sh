#!/usr/bin/bash

mkdir -p target/

# Filter SNPs that are in ATAC peaks.
sh bin/snps_in_peaks.sh

# Merge the peaks and output population-specific statistics.
python bin/peak_merge.py \
       -i data/peaks/all_peaks_sorted.bed.gz \
       > target/merged_peaks.bed

# Compute continental variance between African and European populations.
# Do this for peak reads.
python bin/continental_variance.py target/pop_peak_reads.txt \
       > target/continental_variance/reads.txt
# Do this for max peak heights.
python bin/continental_variance.py target/pop_peak_heights.txt \
       > target/continental_variance/heights.txt

# Get list of SNP rsIDs corresponding to peaks for input into DEPICT.
# Do this for peak reads.
python bin/peak_to_rsid.py \
       data/snps_in_peaks.bed \
       target/continental_variance/reads.txt \
       > target/continental_variance/reads_rsids.txt
# Do this for max peak heights.
python bin/peak_to_rsid.py \
       data/snps_in_peaks.bed \
       target/continental_variance/heights.txt \
       > target/continental_variance/heights_rsids.txt

