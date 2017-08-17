import sys

SIGNAL_TYPE = 'reads'

def load_loci(loci_fname):
    ensid_to_rsids = {}
    with open(loci_fname, 'r') as loci_file:
        for line in loci_file:
            if line.startswith('#'):
                continue
            fields = line.rstrip().split()
            genes = fields[4].split(';')
            rsid = fields[0]
            for gene in genes:
                if not gene in ensid_to_rsids:
                    ensid_to_rsids[gene] = []
                ensid_to_rsids[gene].append(rsid)
    return ensid_to_rsids

def load_rsid_map(rsid_map_fname):
    rsid_to_peak = {}
    with open(rsid_map_fname, 'r') as rsid_map_file:
        for line in rsid_map_file:
            fields = line.rstrip().split()
            rsid = fields[0]
            chrom = fields[1].split(':')[0].replace('chr', '')
            start, end = fields[1].split(':')[1].split('-')
            start, end = int(start), int(end)
            assert(not rsid in rsid_to_peak)
            rsid_to_peak[rsid] = (chrom, start, end)
    return rsid_to_peak

def load_pops(pops_fname):
    peak_to_pops = {}
    with open(pops_fname, 'r') as pops_file:
        for line in pops_file:
            fields = line.rstrip().split()
            chrom = fields[0].replace('chr', '')
            start, end = int(fields[1]), int(fields[2])
            pops = [ float(f) for f in fields[3:] ]
            peak = (chrom, start, end)
            assert(not peak in peak_to_pops)
            peak_to_pops[peak] = pops
    return peak_to_pops

if __name__ == '__main__':
    ensid_fname = sys.argv[1]

    analysis_type = sys.argv[2]
    population = sys.argv[3]

    with open(ensid_fname, 'r') as ensid_file:
        ensids = ensid_file.read().rstrip().split()

    loci_fname = ('depict/results/{}_{}_biased_{}_loci.txt'
                  .format(analysis_type, SIGNAL_TYPE, population))
    ensid_to_rsids = load_loci(loci_fname)

    rsid_map_fname = ('target/{}/{}_biased_{}_rsids.txt'
                      .format(analysis_type, SIGNAL_TYPE, population))
    rsid_to_peak = load_rsid_map(rsid_map_fname)

    pops_fname = 'target/pop_peak_{}.txt'.format(SIGNAL_TYPE)
    peak_to_pops = load_pops(pops_fname)

    for ensid in ensids:
        rsids = ensid_to_rsids[ensid]
        for rsid in rsids:
            peak = rsid_to_peak[rsid]
            pops = peak_to_pops[peak]

            print('{}\t{}\t{}:{}-{}\t{}'.format(
                ensid, rsid, peak[0], peak[1], peak[2],
                '\t'.join([ str(p) for p in pops ])
            ))
    
