# Whole genome capture array design pipeline

This repository contains a pipeline for generating DNA capture probes for capturing
_whole_ chromosomes of human sequences from ancient DNA samples with high rates of
microbial contamination (such that it prohibits normal shotgun sequencing).

The main idea is to take the strictly mappable regions of the hg19 human reference
genome and lay 52 bp probe coordinates along the whole genome (see
`script/calc_probe_coords.py`), with 10 bp tiling gap between each probe. The
sequence of each probe is then extracted from the reference sequence and all probes
are saved to a tab-separated table that can be uploaded to the capture array
manufacturing service.
