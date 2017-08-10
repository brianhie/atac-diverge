import os.path
import pickle
import sys

def closest_dist(snp_a, snp_b, distance):
    if abs(distance - snp_a[0]) < abs(distance - snp_b[0]):
        return snp_a, 0
    else:
        return snp_b, 1

# Implements a binary search on a tuple of (distances, QTLs), where the
# user specifies a `distance' and we query `sorted_distances' for one of
# its elements with the closest distance.
# Implemented using base indices and lengths in order to not use Python's
# inefficient array slicing implementations.
def search_closest(distance, sorted_distances,
                   base=0, length=None):
    if length == None:
        length = len(sorted_distances)
    
    # Impossible to find closest distance in 0 length array,
    # just return None.
    if length == 0:
        return None, base

    # An array of length 1 must have the closest element regardless.
    elif length == 1:
        return sorted_distances[base], base

    # The distance can fall in-between or outside of the elements in an
    # array of length 2.
    elif length == 2:
        assert(sorted_distances[base][0] <= sorted_distances[base+1][0])        
        if sorted_distances[base][0] <= distance and \
           distance >= sorted_distances[base+1][0]:
            closest, idx = closest_dist(
                sorted_distances[base], sorted_distances[base+1],
                distance                
            )
        elif distance < sorted_distances[base][0]:
            closest, idx = sorted_distances[base], 0
        else:
            closest, idx = sorted_distances[base+1], 1
        return closest, base + idx

    # If the array is greater than length 2, do a binary search on half
    # the array.
    else:
        # Split down the middle.
        split = int(length / 2)
        assert(sorted_distances[base+split-1][0] <=
               sorted_distances[base+split][0])

        # Check if either of the elements that we are splitting on are the
        # closest to the distance.
        if sorted_distances[base+split-1][0] <= distance and \
           sorted_distances[base+split][0] >= distance:
            closest, idx = closest_dist(
                sorted_distances[base+split-1], sorted_distances[base+split],
                distance
            )
            if idx == 0:
                return closest, base + split - 1
            else:
                return closest, base + split 
        
        # Else, the distance lies in between elements somewhere in either
        # the left or right array.
        if distance < sorted_distances[base+split-1][0]:
            # Search left side of array.
            closest, idx = search_closest(
                distance, sorted_distances,
                base=base, length=split
            )
        else:
            # Search right side of array.
            closest, idx = search_closest(
                distance,
                sorted_distances,
                base=base+split, length=length-split
            )
        return closest, idx

def load_snps(dbsnp_fname):
    cached_fname = dbsnp_fname + '.pickle'
    # Cache this file to speed up performance.
    if os.path.isfile(cached_fname):
        with open(cached_fname, 'rb') as cached_file:
            snps = pickle.load(cached_file)
        return snps

    # Load a map from chromosomes to positions and rsIDs.
    snps = {}
    with open(dbsnp_fname, 'r') as dbsnp:
        for line in dbsnp:
            fields = line.rstrip().split('\t')
            chrom, pos, rsid = fields[0], int(fields[1]), fields[2]
            if chrom.startswith('chr'):
                chrom = chrom[len('chr'):]
            if not chrom in snps:
                snps[chrom] = []
            snps[chrom].append((pos, rsid))

    # Sort list of SNPs to enable binary search.
    for chrom in snps:
        snps[chrom] = sorted(snps[chrom])

    # Cache the file.
    with open(cached_fname, 'wb') as cached_file:
        pickle.dump(snps, cached_file,
                    protocol=pickle.HIGHEST_PROTOCOL)
    return snps


if __name__ == '__main__':
    dbsnp_fname = sys.argv[1]
    peak_fname = sys.argv[2]

    # Construct map from chromosome to a list of (position, rsID) tuples.
    snps = load_snps(dbsnp_fname)
    
    # Iterate through peak file, finding the closest SNP to the middle of
    # the peak and reporting the rsID of that SNP.
    with open(peak_fname, 'r') as peak_file:
        for line in peak_file:
            fields = line.rstrip().split('\t')
            chrom, start, end = fields[0], int(fields[1]), int(fields[2])
            if chrom.startswith('chr'):
                chrom = chrom[len('chr'):]
            middle = (start + end) / 2

            # Search for SNP that is closest to the middle of the peak.
            closest, closest_idx = search_closest(middle, snps[chrom])
            assert(closest != None)
            assert(closest == snps[chrom][closest_idx])
            closest_pos = closest[0]
            closest_rsid = closest[1]

            # Only allow SNPs within the peak or within 200 bp of the
            # middle of the peak.
            if start <= closest_pos <= end or \
               abs(closest_pos - middle) <= 200:
                print(closest_rsid)
                # Draw without replacement.
                snps[chrom].pop(closest_idx)

