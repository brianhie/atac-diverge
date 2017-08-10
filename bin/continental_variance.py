from numpy import mean
from scipy.stats import ttest_ind
from statsmodels.stats.multitest import multipletests
from subprocess import check_output
import sys

from compare_pops import EURO_POPS, AFRO_POPS, iter_peaks

if __name__ == '__main__':
    infile_name = sys.argv[1]
    
    p_vals = []
    records = []
    with open(infile_name, 'r') as infile:
        for (chrom, start, end), pop_to_val in iter_peaks(infile):
            vals_afr = [
                pop_to_val[p] for p in AFRO_POPS
            ]
            vals_eur = [
                pop_to_val[p] for p in EURO_POPS
            ]
            
            t, p = ttest_ind(vals_afr, vals_eur)
            p_vals.append(p)

            record = [
                chrom, start, end,
                p, mean(vals_afr), mean(vals_eur)
            ]
            record += vals_afr + vals_eur
            records.append(record)

    reject, _, _, _ = multipletests(
        p_vals, alpha=0.05,
        #method='fdr_bh'
        method='bonferroni'
    )

    for p in range(len(p_vals)):
        if reject[p] == True:
            print('\t'.join([ str(f) for f in records[p] ]))
