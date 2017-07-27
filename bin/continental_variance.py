from numpy import mean
from scipy.stats import ttest_ind
from subprocess import check_output
import sys

from compare_pops import EURO_POPS, AFRO_POPS, iter_peaks

if __name__ == '__main__':
    infile_name = sys.argv[1]
    
    line_count = int(check_output(
        "wc -l {0} | sed 's/ .*//'".format(infile_name), shell=True)
                     .decode().rstrip())

    with open(infile_name, 'r') as infile:
        for (chrom, start, end), pop_to_val in iter_peaks(infile):
            vals_afr = [
                pop_to_val[p] for p in AFRO_POPS
            ]
            vals_eur = [
                pop_to_val[p] for p in EURO_POPS
            ]
            t, p = ttest_ind(vals_afr, vals_eur)

            if p * line_count < 0.05:
                fields = [
                    chrom, start, end,
                    p, mean(vals_afr), mean(vals_eur)
                ]
                fields += vals_afr + vals_eur
                print('\t'.join([ str(f) for f in fields ]))
