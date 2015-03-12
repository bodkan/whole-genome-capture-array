######################################################################
# Calculate coordinates of all probes overlapping a given region and
# return the whole thing as an GRanges object.
tile_probes <- function(region, probe_length = 52, tiling_step = 10, flank_length = 35) {
    # initialize vectors for storing start/end positions of all probes
    starts <- c()
    ends <- c()

    region_start <- start(region)
    region_end <- end(region)

    # calculate coordinates of the first probe, which extends by 'flank' bases
    # upstream of the start position of the whole region
    probe_start <- region_start - flank_length
    probe_end <- probe_start + probe_length - 1

    while (probe_end - region_end < flank_length) {
        starts <- c(starts, probe_start)
        ends <- c(ends, probe_end)

        probe_start <- probe_start + tiling_step
        probe_end <- probe_end + tiling_step
    }

    return(GRanges(seqnames = seqnames(region), IRanges(start = starts, end = ends)))
}

######################################################################
# Plot a given region (a single range in a GRanges object) along with
# all probes covering this region (GRanges object)
# -- this function is inspired by plotRanges function as shown in
# "An Introduction to IRanges", a tutorial on IRanges package from
# the Bioconductor project
show_probes <- function(region, probes, compact = TRUE, col = "black", sep = 0.5, ...)
{
    height <- 1

    if (compact)
        bins <- disjointBins(IRanges(start(probes), end(probes) + 1))
    else
        bins <- 1:length(probes)

    plot.new()

    # set xlim to [start_coord_of_first_probe, end_coord_of_last_probe]]
    xlim <- c(start(probes) %>% min, end(probes) %>% max)

    plot.window(xlim, ylim = c(0, (max(bins) + 1)*(height + sep)))

    ybottom <- bins * (sep + height)
    rect(start(probes) - 0.2, ybottom, end(probes) + 0.2, ybottom + height, border = NA, col = col, ...)

    ybottom <- min(ybottom) - sep - height
    rect(start(region) - 0.2, ybottom, end(region) + 0.2, ybottom + height, border = NA, col = "darkred", ...)

    title(paste0("Probes covering region ", seqnames(region), ":", start(region), "-", end(region)))
    axis(1)
}
