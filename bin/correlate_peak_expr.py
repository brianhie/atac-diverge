import json
import numpy as np
from scipy.stats import spearmanr, pearsonr
from statsmodels.stats.multitest import multipletests
import sys

from diff_expr import load_expr
from peak_merge import ALL_POPS
from peak_to_rsid import search_closest

GEUVADIS_POPS = [ 'CEU', 'FIN', 'TSI', 'YRI' ]
EXPR_AGG = np.mean # np.median
CORR_TYPE = pearsonr
P_VAL_CUTOFF = 0.05
MULTI_TEST_METHOD = 'fdr_bh'
DIST_CUTOFF = 50000

def load_tsss(tss_fname):
    tsss = {}
    with open(tss_fname, 'r') as tss_file:
        for line in tss_file:
            fields = line.rstrip().split()
            chrom = fields[0]
            tss = int(fields[4])
            ensid = fields[5]
            symbol = fields[6]

            if not chrom in tsss:
                tsss[chrom] = []
            tsss[chrom].append((tss, (ensid, symbol)))

    # Sort for binary search.
    for chrom in tsss:
        tsss[chrom] = sorted(tsss[chrom])
            
    return tsss

def peak_to_tss(tsss, peak_fname):
    with open(peak_fname, 'r') as peak_file:
        for line in peak_file:
            fields = line.rstrip().split('\t')
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            pops = [ float(f) for f in fields[3:] ]
            if chrom.startswith('chr'):
                chrom = chrom[len('chr'):]
            middle = (start + end) / 2

            # Search for TSS that is closest to the middle of the peak.
            closest, closest_idx = search_closest(middle, tsss[chrom])
            assert(closest != None)
            assert(closest == tsss[chrom][closest_idx])
            closest_pos = closest[0]
            ensid = closest[1][0]
            symbol = closest[1][1]

            # Forward search for TSSs within distance cutoff.
            tss_idx = closest_idx
            while tss_idx < len(tsss[chrom]) and \
                  tsss[chrom][tss_idx][0] - middle <= DIST_CUTOFF:
                ensid, symbol = tsss[chrom][tss_idx][1]
                yield (chrom, start, end, pops,
                       tsss[chrom][tss_idx][0], ensid, symbol)
                tss_idx += 1

            # Backward search.
            tss_idx = closest_idx - 1
            while tss_idx >= 0 and \
                  middle - tsss[chrom][tss_idx][0] <= DIST_CUTOFF:
                ensid, symbol = tsss[chrom][tss_idx][1]
                yield (chrom, start, end, pops,
                       tsss[chrom][tss_idx][0], ensid, symbol)
                tss_idx -= 1

def pop_expr(gene_expr, pops, pop_name):
    return [ float(gene_expr[indiv])
             for indiv in pops[pop_name]
             if indiv in gene_expr ]
                
if __name__ == '__main__':
    tss_fname = sys.argv[1]
    peak_fname = sys.argv[2]
    expr_fname = sys.argv[3]
    pops_fname = sys.argv[4]

    tsss = load_tsss(tss_fname)

    with open(pops_fname, 'r') as pops_file:
        pops = json.loads(pops_file.read())

    gene_to_expr = load_expr(expr_fname)

    pop_to_idx = {}
    for pop_name in GEUVADIS_POPS:
        pop_to_idx[pop_name] = ALL_POPS.index(pop_name)

    p_vals = []
    records = []
    for (chrom, start, end, pop_peaks,
         tss_pos, ensid, symbol) in peak_to_tss(tsss, peak_fname):
        
        if not ensid in gene_to_expr:
            continue
        
        peak_scores = []
        expr_values = []
        for pop_name in GEUVADIS_POPS:
            peak_scores.append(pop_peaks[pop_to_idx[pop_name]])
            expr_values.append(
                EXPR_AGG(pop_expr(gene_to_expr[ensid], pops, pop_name))
            )

        # Exclude list of all zeros since the correlation is undefined in
        # this case.
        if all(x == 0 for x in peak_scores) or \
           all(x == 0 for x in expr_values):
            continue

        rho, p = CORR_TYPE(peak_scores, expr_values)
        p_vals.append(p)
        records.append([
            chrom, start, end, tss_pos, ensid, symbol, rho, p
        ] + peak_scores + expr_values)

    reject, _, _, _ = multipletests(
        p_vals, alpha=P_VAL_CUTOFF,
        method=MULTI_TEST_METHOD
    )

    for p in range(len(p_vals)):
#        if reject[p]:
        print('\t'.join([ str(f) for f in records[p] ]))
