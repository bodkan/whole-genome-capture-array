# default probe-design parameters (in case they're not supplied by the user)
probe_length := 52
tiling_step := 10
flank_length := 34

chromosomes = $(shell seq 1 22) X Y

# input/output directories
raw_data_dir := ./raw_data
clean_data_dir := ./clean_data
output_dir := ./output
figs_dir := ./figs
scripts_dir := ./scripts
tmp_dir := ./tmp
DIRS := $(raw_data_dir) $(clean_data_dir) $(output_dir) $(figs_dir) $(tmp_dir)

# data files produced by the probe-design scrip
probe_coordinates := $(output_dir)/all_probes_$(tiling_step)bp_tiling.bed.gz
unique_regions := $(clean_data_dir)/unique_regions.bed.gz
probe_count := $(tmp_dir)/probe_count_$(tiling_step)bp_tiling.txt
merged_probes := $(tmp_dir)/merged_probes_$(tiling_step)bp_tiling.bed.gz
intersect_with_trf := $(tmp_dir)/intersect_with_trf_$(tiling_step)bp_tiling.bed.gz

# two essential input files
hengs_filter := $(raw_data_dir)/hs37m_filt35_99.bed.gz
trf := $(raw_data_dir)/simpleRepeat.bed.gz

script := $(scripts_dir)/calc_probe_coords.py

.PHONY: all probes figures clean

all: probes figures

probes: $(DIRS) $(probe_coordinates)

figures: $(DIRS) $(probe_count) $(intersect_with_trf)
	Rscript $(scripts_dir)/generate_figures.R $(tiling_step)

$(probe_coordinates): $(unique_regions)
	python3 $(script) \
	    --in_file=$(unique_regions) \
	    --out_file=$(probe_coordinates)_tmp \
	    --probe_length=$(probe_length) \
	    --tiling_step=$(tiling_step) \
	    --flank_length=$(flank_length)
	sort -k1,1V -k2,2n $(probe_coordinates)_tmp | gzip > $(probe_coordinates)
	rm $(probe_coordinates)_tmp

$(intersect_with_trf): $(merged_probes)
	bedtools intersect -a $(merged_probes) \
	                   -b $(trf) | gzip > $(intersect_with_trf)

$(merged_probes):
	bedtools merge -i $(probe_coordinates) | gzip > $@

$(probe_count): $(probe_coordinates)
	printf "chr\tprobes\nall\t" > $@
	zcat $< | wc -l >> $@
	for i in $(chromosomes); do \
	    printf "chr$${i}\t" >> $@; \
	    zgrep -w ^$${i} $< | wc -l >> $@; \
	done

# get coordinates of Heng's alignability regions that don't overlap TRF
$(unique_regions): $(hengs_filter) $(trf)
	bedtools subtract -a $(hengs_filter) -b $(trf) | gzip > $@

# create a local copy of Heng's alignability filter
$(hengs_filter):
	read -p "MPI login: " USERNAME; \
	scp $${USERNAME}@bio23.eva.mpg.de:/mnt/454/HCNDCAM/Hengs_Alignability_Filter/hs37m_filt35_99.bed.gz $@

# download tandem repeat filter from UCSC and extract only coordinate columns
$(trf):
	curl http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz | \
	gunzip | \
	cut -f2,3,4 | \
	sed 's/^chr//' > $@_tmp
	bedtools merge -i $@_tmp | gzip > $@
	rm $@_tmp

# download table of chromosome lengths from UCSC
$(clean_data_dir)/chrom_lengths.txt:
	curl http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/chromInfo.txt.gz | \
	gunzip | \
	cut -f1,2 | \
	grep -w "chr[X,Y,0-9]*" | \
	sort -k1,1V > $@

$(DIRS):
	mkdir $@

clean:
	rm -rf $(DIRS)
