library(magrittr)
library(dplyr)
library(txtplot)
library(rtracklayer)

########################################################################
# get the distance between neighboring regions in a GRanges object
get_gaps <- function(gr) {
    starts <- start(gr[2:length(gr)]) 
    ends <- end(gr[1:length(gr) - 1])

    starts - ends - 1
}

chromosomes <- c(as.character(1:22), "X", "Y")

########################################################################
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
# chromosome?
cat("Min gaps between unique regions per chromosome\n")
min_gap

# chr21 and chr22 contain gaps lower than 35 -> what is the distribution
# gap lengths?
cat("Distribution of gap lengths on chr21\n")
chr21_gaps <- get_gaps(map_filter[["21"]]) %>% table %>% head(n=25) 
chr21_gaps
chr21_gaps %>% txtplot(xlab = "gap length", ylab = "count", width = 80)

cat("Distribution of gap lengths on chr22\n")
chr22_gaps <- get_gaps(map_filter[["22"]]) %>% table %>% head(n=25)  
chr22_gaps
chr22_gaps %>% txtplot(xlab = "gap length", ylab = "count", width = 80)


########################################################################
# load regions that passed Heng Li's alignability filter and which have
# not been found by a Tandem Repeat Finder
no_trf_map_filter <- import.bed("clean_data/mappable_regions.bed.gz")

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
cat("Min gaps between unique regions per chromosome (after TRF filtering)\n")
no_trf_min_gap

# what is the distribution of gap lengths after removing TRF regions?
cat("Distribution of gap lengths on chr21 (after TRF filtering)\n")
no_trf_chr21_gaps <- get_gaps(no_trf_map_filter[["21"]]) %>% table %>% head(n=25) 
no_trf_chr21_gaps 
no_trf_chr21_gaps %>% txtplot(xlab = "gap length", ylab = "count", width = 80)

cat("Distribution of gap lengths on chr22 (after TRF filtering)\n")
no_trf_chr22_gaps <- get_gaps(no_trf_map_filter[["22"]]) %>% table %>% head(n=25)  
no_trf_chr22_gaps
no_trf_chr22_gaps %>% txtplot(xlab = "gap length", ylab = "count", width = 80)

cat("Distribution of gap lengths on chr1 (after TRF filtering)\n")
no_trf_chr1_gaps <- get_gaps(no_trf_map_filter[["1"]]) %>% table %>% head(n=25) 
no_trf_chr1_gaps 
no_trf_chr1_gaps %>% txtplot(xlab = "gap length", ylab = "count", width = 80)


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

cat("Why 25bp gaps exactly? are there no shorter TRF blocks?\n")
trf_counts <- trf %>% width %>% table
trf_counts %>% head(20)
