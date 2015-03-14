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

chromosomes <- c(paste0("chr", 1:22), "chrX", "chrY")

########################################################################
# load regions that passed Heng Li's alignability filter
# (bases where all overlapping 35mers do not match to any other position
# in the genome allowing for up to one mismatch)
map_filter <- import.bed("raw_data/hs37m_filt35_99_with_chr.bed.gz", genome = "hg19")
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
chr21_gaps <- get_gaps(map_filter[["chr21"]]) %>% table %>% head(n=25) 
chr21_gaps
chr21_gaps %>% txtplot(xlab = "gap length", ylab = "count", width = 80)

cat("Distribution of gap lengths on chr22\n")
chr22_gaps <- get_gaps(map_filter[["chr22"]]) %>% table %>% head(n=25)  
chr22_gaps
chr22_gaps %>% txtplot(xlab = "gap length", ylab = "count", width = 80)


########################################################################
# load regions that passed Heng Li's alignability filter and which have
# not been found by a Tandem Repeat Finder
no_trf_map_filter <- import.bed("clean_data/mappable_regions.bed.gz", genome = "hg19")

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
no_trf_chr21_gaps <- get_gaps(no_trf_map_filter[["chr21"]]) %>% table %>% head(n=25) 
no_trf_chr21_gaps 
no_trf_chr21_gaps %>% txtplot(xlab = "gap length", ylab = "count", width = 80)

cat("Distribution of gap lengths on chr22 (after TRF filtering)\n")
no_trf_chr22_gaps <- get_gaps(no_trf_map_filter[["chr22"]]) %>% table %>% head(n=25)  
no_trf_chr22_gaps
no_trf_chr22_gaps %>% txtplot(xlab = "gap length", ylab = "count", width = 80)

cat("Distribution of gap lengths on chr1 (after TRF filtering)\n")
no_trf_chr1_gaps <- get_gaps(no_trf_map_filter[["chr1"]]) %>% table %>% head(n=25) 
no_trf_chr1_gaps 
no_trf_chr1_gaps %>% txtplot(xlab = "gap length", ylab = "count", width = 80)

