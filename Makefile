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
DIRS := $(raw_data_dir) $(clean_data_dir) $(output_dir) $(figs_dir)

summary_file = $(output_dir)/summary_$(tiling_step)bp_tiling.txt
summary_fig = $(figs_dir)/summary_$(tiling_step)bp_tiling.png

.PHONY: probes figures $(DIRS) clean 

all: probes figures 

probes: $(DIRS) $(clean_data_dir)/mappable_regions.bed.gz
	python3 $(scripts_dir)/calc_probe_coords.py \
	    --in_file=$(clean_data_dir)/mappable_regions.bed.gz \
	    --out_file=$(output_dir)/whole_genome_$(tiling_step)bp_tiling.bed_tmp \
	    --probe_length=$(probe_length) \
	    --tiling_step=$(tiling_step) \
	    --flank_length=$(flank_length)
	sort -k1,1V -k2,2n $(output_dir)/whole_genome_$(tiling_step)bp_tiling.bed_tmp | \
	    gzip > $(output_dir)/whole_genome_$(tiling_step)bp_tiling.bed.gz
	rm $(output_dir)/whole_genome_$(tiling_step)bp_tiling.bed_tmp
	for i in $(chromosomes); do \
	    zgrep -w ^chr$${i} $(output_dir)/whole_genome_$(tiling_step)bp_tiling.bed.gz | \
	    gzip > $(output_dir)/chr$${i}_$(tiling_step)bp_tiling.bed.gz; \
	done

figures: $(DIRS) $(summary_file)
	Rscript $(scripts_dir)/plot_probes_summary.R $(summary_file) $(summary_fig)

$(summary_file):
	printf "chr\tprobes\n" > $(summary_file)
	printf "all\t" >> $(summary_file)
	zcat $(output_dir)/whole_genome_$(tiling_step)bp_tiling.bed.gz | \
	    wc -l >> $(summary_file)
	for i in $(chromosomes); do \
	    printf "chr$${i}\t" >> $(summary_file); \
	    zcat $(output_dir)/chr$${i}_$(tiling_step)bp_tiling.bed.gz | \
	        wc -l >> $(summary_file); \
	done

# get coordinates of Heng's alignability regions that don't overlap TRF
$(clean_data_dir)/mappable_regions.bed.gz: $(raw_data_dir)/hs37m_filt35_99_with_chr.bed.gz \
                                           $(raw_data_dir)/simpleRepeat.bed.gz
	bedtools subtract -a $(raw_data_dir)/hs37m_filt35_99_with_chr.bed.gz \
	                  -b $(raw_data_dir)/simpleRepeat.bed.gz | gzip > $@

# create a local copy of Heng's alignability filter
$(raw_data_dir)/hs37m_filt35_99_with_chr.bed.gz:
	read -p "MPI login: " USERNAME; \
	scp $${USERNAME}@bio23.eva.mpg.de:/mnt/454/HCNDCAM/Hengs_Alignability_Filter/hs37m_filt35_99.bed.gz $(raw_data_dir)
	zcat $(raw_data_dir)/hs37m_filt35_99.bed.gz | \
	sed 's/^/chr/' | gzip > $@

# download tandem repeat filter from UCSC and extract only coordinate columns
$(raw_data_dir)/simpleRepeat.bed.gz:
	curl http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/simpleRepeat.txt.gz | \
	gunzip | \
	cut -f2,3,4 | \
	gzip > $@

clean:
	rm -rf $(DIRS)

$(DIRS):
	mkdir -p $@
