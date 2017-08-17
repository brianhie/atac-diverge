import json
from math import log
import numpy as np
from scipy.stats import ttest_ind
import sys

from diff_expr import load_col, load_expr

def load_fold_diffs(genes, pops, gene_to_expr):
    genes_fold_diffs = []
    
    for gene in genes:
        if not gene in gene_to_expr:
            sys.stderr.write('Warning: Could not find gene {}\n'
                             .format(gene))
            continue

        pops_expr = []
        for p in range(len(pops)):
            pop = list(set(pops[str(p)]))
            pop_expr = []
            for indiv in pop:
                if indiv in gene_to_expr[gene]:
                    pop_expr.append(float(gene_to_expr[gene][indiv]))
            pops_expr.append(pop_expr)

        try:
            genes_fold_diffs.append(
                log(np.mean(pops_expr[0]) / np.mean(pops_expr[1]))
            )
        except:
            continue
        
    return genes_fold_diffs

if __name__ == '__main__':
    expr_fname = sys.argv[1]
    pops_fname = sys.argv[2]
    gene_fname = sys.argv[3]
    if len(sys.argv) < 5:
        background_fname = 'depict/data/Genes.txt'
    else:
        background_fname = sys.argv[4]        

    gene_to_expr = load_expr(expr_fname)
    with open(pops_fname, 'r') as pops_file:
        pops = json.loads(pops_file.read())
    genes = load_col(gene_fname, 0)
    background = load_col(background_fname, 0)

    candidate_fold_diffs = load_fold_diffs(genes,
                                           pops, gene_to_expr)
    background_fold_diffs = load_fold_diffs(background,
                                            pops, gene_to_expr)

    for i, fd in enumerate([
            candidate_fold_diffs, background_fold_diffs
    ]):
        print('Pop {}:\t{:.2f}|{:.2f}|{:.2f}'.format(
            i,
            np.percentile(fd, 25),
            np.percentile(fd, 50),
            np.percentile(fd, 75)
        ))
    print('t-test p = {}'.format(
        ttest_ind(candidate_fold_diffs, background_fold_diffs,
                  equal_var=True)[1]
    ))
