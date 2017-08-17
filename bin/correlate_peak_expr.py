import json
import numpy as np
from scipy.stats import spearmar
import sys

from diff_expr import load_expr
from peak_merge import ALL_POPS
from peak_to_rsid import search_closest

GEUVADIS_POPS = [ "CEU", "FIN", "GBR", "TSI", "YRI" ]
EXPR_AGG = np.median
CORR_TYPE = spearmanr

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
            tss_pos = closest_pos
            while tss_idx < len(tsss[chrom]) and \
                  tss_pos - middle <= 100000:
                yield (chrom, start, end, pops,
                       tss_pos, ensid, symbol)
                tss_idx += 1
                tss_pos = tsss[chrom][tss_idx][0]
                ensid, symbol = tsss[chrom][tss_idx][1]

            # Backward search.
            tss_idx = closest_idx - 1
            tss_pos = tsss[chrom][tss_idx][0]
            while tss_idx >= 0 and \
                  middle - tss_pos <= 100000:
                yield (chrom, start, end, pops,
                       tss_pos, ensid, symbol)
                tss_idx -= 1
                tss_pos = tsss[chrom][tss_idx][0]
                ensid, symbol = tsss[chrom][tss_idx][1]

def pop_expr(gene_to_expr, pops, pop_name):
    expr = []
    for indiv in pops[pop_name]:
        if indiv in gene_to_expr[gene]:
            expr.append(float(gene_to_expr[gene][indiv]))
    return expr
    
                
if __name__ == '__main__':
    tss_fname = sys.argv[1]
    peak_fname = sys.argv[2]
    expr_fname = sys.argv[3]
    pops_fname = sys.argv[4]

    tsss = load_tsss(tss_fname)

    with open(pops_fname, 'r') as pops_file:
        pops = json.loads(pops_file.read())

    gene_to_expr = load_expr(expr_fname)

    p_vals = []
    records = []
    for (chrom, start, end, pops,
         tss_pos, ensid, symbol) in peak_to_tss(tsss, peak_fname):
        peak_scores = []
        expr_values = []
        for pop_name in GEUVADIS_POPS:
            idx = ALL_POPS.index(pop_name)
            peak_scores.append(pops[idx])
            expr_values.append(
                EXPR_AGG(pop_expr(gene_to_expr, pops, pop_name))
            )

        rho, p = CORR_TYPE(peak_scores, expr_values)
        p_vals.append(p)
        records.append([
            chrom, start, end, tss_pos, ensid, symbol
        ] + peak_scores + expr_values)

    reject, _, _, _ = multipletests(
        p_vals, alpha=P_VAL_CUTOFF,
        method=MULTI_TEST_METHOD
    )

    for p in range(len(p_vals)):
        if reject[p]:
            print('\t'.join([ str(f) for f in records[p] ]))
