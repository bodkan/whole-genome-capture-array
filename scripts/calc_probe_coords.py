#!/usr/bin/env python3

import gzip
import argparse

def tile_probes(start, end, probe_length, tiling_step, flank_length):
    """ Calculate coordinates of all probes overlapping a given region and
    return the whole thing as a list of pairs (probe_start, probe_end).
    """
    # the first probe overlapping a given region extends 'flank_length' bases
    # upstream from the start of that region...
    probe_start = start - flank_length
    # ...however, in case the region starts less than 'flank_length' bases
    # from the beginning of a chromosome, set the start of the first probe to 0
    if probe_start < 0:
        probe_start = 0

    probe_end = probe_start + probe_length

    probes = []

    while (probe_end - end) <= flank_length:
        probes.append((probe_start, probe_end))

        probe_start = probe_start + tiling_step
        probe_end = probe_end + tiling_step

    return probes

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--in_file", help="gzipped input BED file", required=True, type=str)
    parser.add_argument("--out_file", help="output BED file", required=True, type=str)
    parser.add_argument("--probe_length", help="probe size [bp]", required=True, type=int)
    parser.add_argument("--tiling_step", help="tiling density [bp]", required=True, type=int)
    parser.add_argument("--flank_length", help="max length of flanking regions [bp]", required=True, type=int)
    args = parser.parse_args()

    with gzip.open(args.in_file, "rt") as map_filter_file, open(args.out_file, "w") as output_file:

        for line in map_filter_file:
            fields = line.strip().split("\t")

            chrom = fields[0]
            start = int(fields[1])
            end = int(fields[2])

            probe_coords = tile_probes(start, end, args.probe_length, args.tiling_step, args.flank_length)
            for p in probe_coords:
                print(chrom, p[0], p[1], sep="\t", file=output_file)
