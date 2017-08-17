import sys

from peak_to_rsid import search_closest

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
                yield (chrom, start, end, tss_pos,
                       ensid, symbol)
                tss_idx += 1
                tss_pos = tsss[chrom][tss_idx][0]
                ensid, symbol = tsss[chrom][tss_idx][1]

            # Backward search.
            tss_idx = closest_idx - 1
            tss_pos = tsss[chrom][tss_idx][0]
            while tss_idx >= 0 and \
                  middle - tss_pos <= 100000:
                yield (chrom, start, end, tss_pos,
                       ensid, symbol)
                tss_idx -= 1
                tss_pos = tsss[chrom][tss_idx][0]
                ensid, symbol = tsss[chrom][tss_idx][1]
                
if __name__ == '__main__':
    tss_fname = sys.argv[1]
    peak_fname = sys.argv[2]

    tsss = load_tsss(tss_fname)

    for (chrom, start, end, tss_pos,
         ensid, symbol) in peak_to_tss(tssss, peak_fname):
        pass
