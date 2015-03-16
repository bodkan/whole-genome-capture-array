library(magrittr)
library(rtracklayer)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    stop("You must supply a file table of probe counts per chromosome
         and a tiling step size!")
}

input_table <- args[1]
tiling_step <- as.integer(args[2])

# total number of probes
df <- read.delim(input_table)
png(paste0("figs/number_of_probes_", tiling_step, ".png"), width = 1280, height = 800, res = 100)
par(las = 2)
barplot(df[,2] / 1000000, names = df[,1], log = "y", ylab = "number of probes [millions]",
        border = NA, main = paste0("Number of probes per chromosome [", tiling_step, "bp tiling]"))
dev.off()

# distribution of lengths of covered regions
#bedtools merge -i whole_genome_10bp_tiling.bed.gz | gzip > merged.bed.gz
fragments <- import.bed("output/merged.bed.gz")
fragment_lengths <- width(fragments)

png(paste0("figs/fragment_lengths_", tiling_step, "bp.png"), width = 1280, height = 800, res = 100)
hist(fragment_lengths, col = "black", border = NA, breaks = 300,
     main = "Size distribution of regions covered by probes", xlab = "length of a region [bp]")
dev.off()

quantile(fragment_lengths, seq(0.1, 1, 0.1))

