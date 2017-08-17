import gzip
import json
import numpy as np
import pickle
from scipy.stats import ttest_ind
import sys

verbose = False

def load_col(fname, col_pos, col_type=str):
    col = []
    with open(fname, 'r') as f:
        for line in f:
            fields = line.rstrip().split()
            col.append(col_type(fields[col_pos]))
    return col

def load_expr(fname, gene_pos=0):
    gene_to_expr = {}

    is_gzip = fname.endswith('.gz')
    opener = gzip.open if is_gzip else open
    with opener(fname, 'r') as expr_file:
        header = None
        for line in expr_file:
            if is_gzip:
                line = line.decode('utf-8')
            
            if line.startswith('#'):
                continue
            if header == None:
                header = line.rstrip().split()
                continue

            fields = line.rstrip().split()
            gene = fields[gene_pos]
            if '.' in gene:
                gene = gene.split('.')[0]
                
            gene_to_expr[gene] = {}
            for pos, h in enumerate(header):
                if pos == gene_pos:
                    continue
                gene_to_expr[gene][h] = fields[pos]
                
    return gene_to_expr

if __name__ == '__main__':
    expr_fname = sys.argv[1]
    pops_fname = sys.argv[2]
    gene_fname = sys.argv[3]

    with open('data/ensg_to_gene_symbol.pickle', 'rb') as f:
        ensg_to_gene_symbol = pickle.load(f)

    with open(pops_fname, 'r') as pops_file:
        pops = json.loads(pops_file.read())

    genes = load_col(gene_fname, 0)

    gene_to_expr = load_expr(expr_fname)

    for gene in genes:
        if not gene in gene_to_expr:
            if verbose:
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

        if gene in ensg_to_gene_symbol:
            symbol = ensg_to_gene_symbol[gene]
        else:
            symbol = gene
            
        sys.stdout.write('{}\t{}'.format(gene, symbol))
        for pop_expr in pops_expr:
            sys.stdout.write('\t{:.2f}'.format(
                np.percentile(pop_expr, 50),
            ))
        sys.stdout.write('\t{}'.format(
            ttest_ind(pops_expr[0], pops_expr[1],
                      equal_var=True)[1]
        ))
        sys.stdout.write('\n')
