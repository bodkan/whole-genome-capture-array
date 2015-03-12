args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
    stop("You must supply a file table of probe counts per chromosome
         and a path to an output plot file!")
}

input_table <- args[1]
output_plot <- args[2]

df <- read.delim(input_table)

png(output_plot, width = 1280, height = 800)
par(las = 2)

barplot(df[,2] / 1000000, names = df[,1], log = "y", ylab = "number of probes [millions]",
        border = NA, main = "Number of probes per chromosome [$(tiling_step)bp tiling]")

dev.off()

