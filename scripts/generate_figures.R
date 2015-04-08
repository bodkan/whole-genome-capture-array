args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 1) {
    stop("You must supply length of a tiling step (in base pairs)!")
}

tiling_step <- as.integer(args[1])
chromosomes <- c(as.character(1:22), "X", "Y")

library(magrittr)
library(dplyr)
library(rtracklayer)


########################################################################
# summary plots

# plot total number of probes
df <- read.delim(paste0("tmp/probe_count_", tiling_step, "bp_tiling.txt"))
png(paste0("figs/number_of_probes_", tiling_step, ".png"), width = 1280, height = 800, res = 100)
par(las = 2)
barplot(df[,2] / 1000000, names = df[,1], log = "y", ylab = "number of probes [millions]",
        border = NA, main = paste0("Number of probes per chromosome [", tiling_step, "bp tiling]"))
dev.off()


# plot distribution of lengths of covered regions
#quantile(fragment_lengths, seq(0.1, 1, 0.1))
fragments <- import.bed(paste0("tmp/merged_final_probes_", tiling_step, "bp_tiling.bed.gz"))
fragment_lengths <- width(fragments)

png(paste0("figs/fragment_lengths_", tiling_step, "bp.png"), width = 1280, height = 800, res = 100)
hist(fragment_lengths, col = "gray", border = NA, breaks = 300,
     main = "Size distribution of regions covered by probes", xlab = "length of a region [bp]")
dev.off()



########################################################################
# analysis of gaps between unique regions


# get the distance between neighboring regions in a GRanges object
get_gaps <- function(gr) {
    starts <- start(gr[2:length(gr)]) 
    ends <- end(gr[1:length(gr) - 1])

    starts - ends - 1
}


# load regions that passed Heng Li's alignability filter
# (bases where all overlapping 35mers do not match to any other position
# in the genome allowing for up to one mismatch)
map_filter <- import.bed("raw_data/hs37m_filt35_99.bed.gz")
seqlevels(map_filter, force = TRUE) <- chromosomes

# split the loaded GRanges object into GRangesList -- regions per chromosomes
map_filter <- split(map_filter, seqnames(map_filter))

# what is the minimal distance between neighboring unique regions?
min_gap <- data.frame(
    chr = chromosomes,
    min_gap = sapply(map_filter, function(chr) chr %>% get_gaps %>% min),
    row.names = NULL
)
# by definition, this has to be always more than 35, is it lower for some
# chromosomes?
png("figs/min_gaps_per_chr_before_TRF.png", width = 1280, height = 800, res = 100)
barplot(min_gap$min_gap, names.arg = min_gap$chr, border = NA, ylim = c(0, 35),
        main = "Minimum distance between any two unique regions per chromosome\n(before TRF removal)",
        xlab = "chromosome", ylab = "minimum distance between any two unique regions [bp]")
dev.off()

# chr21 and chr22 contain gaps lower than 35 -> what is the distribution
# gap lengths?
png("figs/gap_lengths_chr21_chr22_before_TRF.png", width = 1280, height = 800, res = 100)
par(mfrow = c(2, 1))

# distribution of gap lengths on chr21
chr21_gaps <- get_gaps(map_filter[["21"]]) %>% table
chr21_gaps <- as.data.frame(chr21_gaps) %>% head(25)
names(chr21_gaps) <- c("gap_length", "count")
barplot(chr21_gaps$count, names.arg = chr21_gaps$gap_length, border = NA,
        main = "Distribution of distances between any two unique regions for chromosome 21 (before TRF removal)",
        xlab = "distance between two unique regions [bp]", ylab = "count")

# distribution of gap lengths on chr22
chr22_gaps <- get_gaps(map_filter[["22"]]) %>% table
chr22_gaps <- as.data.frame(chr22_gaps) %>% head(25)
names(chr22_gaps) <- c("gap_length", "count")
barplot(chr22_gaps$count, names.arg = chr22_gaps$gap_length, border = NA,
        main = "Distribution of distances between any two unique regions for chromosome 22 (before TRF removal)",
        xlab = "distance between two unique regions [bp]", ylab = "count")

dev.off()


# load regions that passed Heng Li's alignability filter and which have
# not been found by a Tandem Repeat Finder
no_trf_map_filter <- import.bed("clean_data/unique_regions.bed.gz")

# split the loaded GRanges object into GRangesList -- regions per chromosomes
seqlevels(no_trf_map_filter, force = TRUE) <- chromosomes
no_trf_map_filter <- split(no_trf_map_filter, seqnames(no_trf_map_filter))

# what is the minimal distance between neighboring unique regions?
# how do these values change compared to just the original Heng's filter?
no_trf_min_gap <- data.frame(
    chr = chromosomes,
    min_gap = sapply(no_trf_map_filter, function(chr) chr %>% get_gaps %>% min),
    row.names = NULL
)
# min gaps between unique regions per chromosome (after TRF filtering)
png("figs/min_gaps_per_chr_after_TRF.png", width = 1280, height = 800, res = 100)
barplot(no_trf_min_gap$min_gap, names.arg = no_trf_min_gap$chr, border = NA, ylim = c(0, 35),
        main = "Minimum distance between any two unique regions per chromosome\n(after TRF removal)",
        xlab = "chromosome", ylab = "minimum distance between any two unique regions [bp]")
dev.off()

# what is the distribution of gap lengths after removing TRF regions?
png("figs/gap_lengths_chr21_chr22_after_TRF.png", width = 1280, height = 800, res = 100)
par(mfrow = c(2, 1))

# distribution of gap lengths on chr21 (after TRF filtering)
no_trf_chr21_gaps <- get_gaps(no_trf_map_filter[["21"]]) %>% table
no_trf_chr21_gaps <- as.data.frame(no_trf_chr21_gaps) %>% head(25)
names(no_trf_chr21_gaps) <- c("gap_length", "count")
barplot(no_trf_chr21_gaps$count, names.arg = no_trf_chr21_gaps$gap_length, border = NA,
        main = "Distribution of distances between any two unique regions for chromosome 21 (after TRF removal)",
        xlab = "distance between two unique regions [bp]", ylab = "count")

# distribution of gap lengths on chr22 (after TRF filtering)
no_trf_chr22_gaps <- get_gaps(no_trf_map_filter[["22"]]) %>% table
no_trf_chr22_gaps <- as.data.frame(no_trf_chr22_gaps) %>% head(25)
names(no_trf_chr22_gaps) <- c("gap_length", "count")
barplot(no_trf_chr22_gaps$count, names.arg = no_trf_chr22_gaps$gap_length, border = NA,
        main = "Distribution of distances between any two unique regions for chromosome 22 (after TRF removal)",
        xlab = "distance between two unique regions [bp]", ylab = "count")

dev.off()


########################################################################
# As shown above, the minimum distance between two neighboring Heng Li's
# unique regions gets reduced from 35bp to 25bp after subtraction of
# repetitive blocks found by TRF.
# The only possible reason for this (other than some weird bug in my
# code) is that some TR regions split unique regions in half, creating
# gaps of _exactly_ 25bp. Let's test this...

# load TRF filter data
trf <- import.bed("raw_data/simpleRepeat.bed.gz")

# find the first unique region on chr1 that is followed by 25bp gap
# (that is, probably a first half of a larger unique region before
# removal of TRF blocks)
i <- which(get_gaps(no_trf_map_filter[["1"]]) == 25)[1]
# start/end coordinate of this region
s <- start(no_trf_map_filter[["1"]][i])
e <- end(no_trf_map_filter[["1"]][i])

cat("Example of a region before subtraction of TRF blocks:\n")
map_filter[["1"]][start(map_filter[["1"]]) == s]
cat("The same region after subtraction of TRF blocks:\n")
no_trf_map_filter[["1"]][c(i, i+1)]
cat("Corresponding TRF block\n")
trf[start(trf) == e + 1]

# why 25bp gaps exactly? are there no shorter TRF blocks?
trf_count <- trf %>% width %>% table %>% as.data.frame %>% head(25)
names(trf_count) <- c("repeat_length", "count")
png("figs/trf_lengths.png", width = 1280, height = 800, res = 100)
barplot(trf_count$count, names.arg = trf_count$repeat_length, border = NA,
        main = "Distribution of lengths of TRF blocks",
        xlab = "repeat length [bp]", ylab = "count")
dev.off()

# whole probe-covered regions overlapped by TRFs
trf_intersect <- import.bed(paste0("tmp/intersect_with_trf_", tiling_step, "bp_tiling.bed.gz"))
trf_intersect <- width(trf_intersect) %>% table %>% as.data.frame
names(trf_intersect) <- c("intersect_length", "count")
# => the overlap with TRFs ranges from 1 to 67 at max
png("figs/trf_intersect.png", width = 1280, height = 800, res = 100)
barplot(trf_intersect$count, names.arg = trf_intersect$intersect_length, border = NA,
        main = "Degree of overlap of 'flanking' regions with TRF",
        xlab = "length of overlap with TRF region [bp]", ylab = "count")
dev.off()

