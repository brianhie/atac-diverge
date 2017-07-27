from peak_merge import ALL_POPS

AFRO_POPS = sorted([ 'ESN', 'GWD', 'LWK', 'YRI'])
EURO_POPS = sorted([ 'CEU', 'FIN', 'IBS', 'TSI'])

def iter_peaks(infile):
    for line in infile:
        base_idx = 3
        fields = line.rstrip().split('\t')
        chrom, start, end = fields[0], int(fields[1]), int(fields[2])
        
        pop_to_val = {}
        for i, pop in enumerate(ALL_POPS):
            pop_to_val[pop] = float(fields[base_idx + i])
            
        yield (chrom, start, end), pop_to_val

