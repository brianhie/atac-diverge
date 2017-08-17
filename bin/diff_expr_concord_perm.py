import json
import numpy as np
import random
import sys

from diff_expr import load_col, load_expr

N_PERMUTATIONS = 1000

if __name__ == '__main__':
    expr_fname = sys.argv[1]
    pops_fname = sys.argv[2]
    n_genes = int(sys.argv[3])
    n_one_direction = int(sys.argv[4])

    assert(n_one_direction <= n_genes)

    with open(pops_fname, 'r') as pops_file:
        pops = json.loads(pops_file.read())

    gene_to_expr = load_expr(expr_fname)

    n = 0
    for _ in range(N_PERMUTATIONS):
        n_one_direction_rand = 0
        for _ in range(n_genes):
            gene = random.choice(list(gene_to_expr.keys()))

            pops_expr = []
            for p in range(len(pops)):
                pop = list(set(pops[str(p)]))
                pop_expr = []
                for indiv in pop:
                    if indiv in gene_to_expr[gene]:
                        pop_expr.append(float(gene_to_expr[gene][indiv]))
                pops_expr.append(pop_expr)

            if np.median(pops_expr[0]) > np.median(pops_expr[1]):
                n_one_direction_rand += 1
        if n_one_direction_rand >= n_one_direction:
            n += 1

    print('p = {}'.format(n / N_PERMUTATIONS))
