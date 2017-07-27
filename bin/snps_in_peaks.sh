#!/usr/bin/bash

# Filter out SNPs that are not in ATAC peaks.
cat /godot/dbsnp/149_GRCh37/dbSNP149_GRCh37_chroms_only.bed | \
    tail -n+2 | \
    sort -k1,1 -k2,2n | \
    /modules/pkgs/bedtools/2.25.0/bin/bedtools \
        intersect \
        -a stdin \
        -b data/peaks/merged_peaks.bed \
        -wa -sorted \
        > data/snps_in_peaks.bed
