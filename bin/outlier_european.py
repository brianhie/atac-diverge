from numpy import mean, std
from scipy.stats import norm
from statsmodels.stats.multitest import multipletests
import sys

from compare_pops import EURO_POPS, AFRO_POPS, iter_peaks

P_VAL_CUTOFF = 0.05
MULTI_TEST_METHOD = 'bonferroni'

if __name__ == '__main__':
    infile_name = sys.argv[1]
    pop_name = sys.argv[2].upper()

    p_vals = []
    records = []
    with open(infile_name, 'r') as infile:
        for (chrom, start, end), pop_to_val in iter_peaks(infile):
            val_pop = pop_to_val[pop_name]
            vals_eur = [
                pop_to_val[p] for p in EURO_POPS
            ]
            # Demand some amount of ATAC peak signal.
            if sum(vals_eur) == 0 or val_pop == 0:
                continue

            sigma = std(vals_eur)
            # z-score will always be positive.
            z = abs(val_pop - mean(vals_eur)) / sigma
            # One-sided p-value.
            p = 1 - norm.cdf(z)
            p_vals.append(p)

            record = [
                chrom, start, end,
                p, val_pop, mean(vals_eur), z,
            ]
            record += vals_eur
            records.append(record)

    reject, _, _, _ = multipletests(
        p_vals, alpha=P_VAL_CUTOFF,
        method=MULTI_TEST_METHOD
    )

    for p in range(len(p_vals)):
        if reject[p]:
            print('\t'.join([ str(f) for f in records[p] ]))
