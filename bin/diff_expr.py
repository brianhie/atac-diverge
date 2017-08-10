import gzip
import json
import numpy as np
import sys

def load_col(fname, col_pos, col_type=str):
    col = []
    with open(fname, 'r') as f:
        for line in f:
            fields = line.rstrip().split()
            col.append(col_type(fields[col_pos]))
    return col

def load_expr(fname, gene_pos=0):
    gene_to_expr = {}

    is_gzip = expr_fname.endswith('.gz')
    opener = gzip.open if is_gzip else open
    with opener(expr_fname, 'r') as expr_file:
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

    with open(pops_fname, 'r') as pos_file:
        pops = json.loads(pops_file.read())

    genes = load_col(gene_fname, _)

    gene_to_expr = load_expr(expr_fname)

    for gene in genes:
        if not gene in gene_to_expr:
            sys.stderr.write('Warning: Could not find gene {}\n'
                             .format(gene))
            continue
        
        pops_expr = []
        for p in range(len(pops)):
            pop = pops[p]
            pop_expr = []
            for indiv in pops:
                if indiv in gene_to_expr[gene]:
                    pop_expr.append(float(gene_to_expr[gene][indiv]))
            pops_expr.append(pop_expr)

        sys.stdout.write('{}'.format('gene'))
        for pop_expr in pops_expr:
            sys.stdout.write('\t{}/{}/{}'.format(
                np.percentile(pop_expr, 25),
                np.percentile(pop_expr, 50),
                np.percentile(pop_expr, 75)
            ))
        sys.stdout.write('\n')

        