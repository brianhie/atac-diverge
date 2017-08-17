import gzip
import sys

def parse_meta(meta_str):
    meta = {}
    for elem in meta_str.rstrip(';').split(';'):
        name, val = elem.strip().split(' ')
        val = val.strip('"')
        meta[name] = val
    return meta

if __name__ == '__main__':
    gtf_fname = sys.argv[1]

    if gtf_fname.endswith('.gz'):
        opener = gzip.open
    else:
        opener = open

    with opener(gtf_fname) as gtf_file:
        for line in gtf_file:
            if gtf_fname.endswith('.gz'):
                line = line.decode('utf-8')
            if line.startswith('#'):
                continue
            fields = line.rstrip().split('\t')

            if not 'protein_coding' in fields[1]:
                continue

            if fields[2] == 'transcript':
                chrom = fields[0].replace('chr', '')
                start = int(fields[3])
                end = int(fields[4])
                strand = fields[6]
                meta = parse_meta(fields[8])
                ensid = meta['gene_id']
                symbol = meta['gene_name']
                
            elif fields[2] == 'start_codon':
                start_codon_start = int(fields[3])
                start_codon_end = int(fields[4])
                if strand == '+':
                    tss = start_codon_start
                elif strand == '-':
                    tss = start_codon_end
                else:
                    assert(False)
                assert(start <= tss <= end)

                ofields = [
                    chrom, start, end, strand, tss,
                    ensid, symbol
                ]
                print('\t'.join([ str(f) for f in ofields ]))
