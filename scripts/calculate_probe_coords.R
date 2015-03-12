#!/usr/bin/env Rscript

######################################################################
# argument parsing
suppressMessages(library(docopt))

"Calculate coordinates of whole-genome capture probes based on given parameters.

Usage:

    calculate_probe_coords.R [options]

Options:

    --chr=<CHROMOSOME>    chromosome id in 'chr[X,Y,0-9]' format

    --probe_length=<SIZE> probe size [# base pairs]

    --tiling_step=<TILING>     tiling density [# base pairs]

    --flank_length=<BY>   max length of flanking regions [# base pairs]

    -h --help             show this help

Output:

    Coordinates of capture probes in a BED file format." -> doc

my_opts <- docopt(doc)

# are all required arguments present?
if (with(my_opts, is.null(chr) | is.null(probe_length) | is.null(tiling_step) | is.null(flank_length))) {
  writeLines(doc)
  quit()
}

# extract command line arguments
chr <- my_opts$chr
probe_length  <- as.integer(my_opts$probe_length)
tiling_step <- as.integer(my_opts$tiling)
flank_length <- as.integer(my_opts$flank_length)


######################################################################
# initialization 

#chr <- "chr1"
#probe_length <- 52
#tiling_step <- 10
#flank_length <- 35

library(magrittr)
library(parallel)
suppressMessages(library(GenomicRanges))
suppressMessages(library(rtracklayer))

source("scripts/utils.R")

setwd("~/Projects/arrays")
input_dir <- "clean_data"
output_dir <- "output"

# load coordinates of mappable regions of the genome
hengs_filter_file <- file.path(input_dir, "mappable_regions.bed.gz")
mappable_regions <- import.bed(hengs_filter_file, genome = "hg19")

seqlevels(mappable_regions, force = TRUE) <- chr

# initialize the local cluster
cl <- makeCluster(detectCores())
clusterExport(cl, varlist = c("tile_probes", "probe_length", "tiling_step", "flank_length"))

# calculate coordinates of probes tiled over all regions of interest
probes_per_region <-
    parLapply(cl, mappable_regions,
              function(r) tile_probes(r, probe_length, tiling_step, flank_length))

stopCluster(cl)

# convert the list of GRanges (~list of probes per region) into GRangesList
# and flatten the GRangesList into a single GRanges object containing all
# probes for the whole genome
genome_wide_probes <- GRangesList(probes_per_region) %>% unlist

export.bed(genome_wide_probes, paste0("output/r_", chr, "_probes.bed"), ignore.strand = TRUE)

