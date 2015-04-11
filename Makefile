# default probe-design parameters (in case they're not supplied by the user)
probe_length := 52
tiling_step := 10
flank_length := 34

chromosomes = $(shell seq 1 22) X Y

scripts_dir := ./scripts
output_dir := ./output
raw_data_dir := ./raw_data
clean_data_dir := ./clean_data
figs_dir := ./figs
tmp_dir := ./tmp
DIRS := $(raw_data_dir) $(clean_data_dir) $(output_dir) $(figs_dir) $(tmp_dir)

probe_design_script := $(scripts_dir)/calc_probe_coords.py
plotting_script := $(scripts_dir)/generate_figures.R

# input files
ref_genome := /mnt/solexa/Genomes/hg19_evan/whole_genome.fa
alignability_filter := $(raw_data_dir)/hs37m_filt35_99.bed.gz
trf := $(raw_data_dir)/simpleRepeat.bed.gz
unique_regions := $(clean_data_dir)/unique_regions.bed.gz

# sequences of the final probe set
final_sequences := $(output_dir)/final_sequences_$(tiling_step)bp_tiling.txt.gz

# intermediate files produced during the calculation of probe coordinates
final_coordinates := $(tmp_dir)/final_coordinates_$(tiling_step)bp_tiling.bed.gz
unfiltered_coordinates := $(tmp_dir)/unfiltered_coordinates_$(tiling_step)bp_tiling.bed.gz

# temporary files used for summary and plotting
probe_count := $(tmp_dir)/probe_count_$(tiling_step)bp_tiling.txt
merged_final_probes := $(tmp_dir)/merged_final_probes_$(tiling_step)bp_tiling.bed.gz
merged_unfiltered_probes := $(tmp_dir)/merged_unfiltered_probes_$(tiling_step)bp_tiling.bed.gz
intersect_with_trf := $(tmp_dir)/intersect_with_trf_$(tiling_step)bp_tiling.bed.gz

.PHONY: all probes figures clean

all: probes figures

probes: $(DIRS) $(final_sequences)

figures: $(DIRS) $(probe_count) $(merged_final_probes) $(intersect_with_trf)
	Rscript $(plotting_script) $(tiling_step)

# get sequences of the final set of probes and exclude probes with 'N' bases
$(final_sequences): $(final_coordinates)
	bedtools getfasta -fi $(ref_genome) -bed $< -fo $@_withNs -tab
	grep -v "N" $@_withNs | sed 's/$$/CACTGCGG/' | gzip > $@
	rm $@_withNs

# exclude all probes that overlap repetitive regions found by TRF
$(final_coordinates): $(unfiltered_coordinates) $(trf)
	bedtools intersect -a $(unfiltered_coordinates) -b $(trf) -c | \
	    awk '($$4 == 0)' | cut -f1,2,3 | gzip > $@

# calculate coordinates of probes overlapping unique regions of human genome
$(unfiltered_coordinates): $(unique_regions)
	python3 $(probe_design_script) \
	    --in_file=$< \
	    --out_file=$@_tmp \
	    --probe_length=$(probe_length) \
	    --tiling_step=$(tiling_step) \
	    --flank_length=$(flank_length)
	sort -k1,1V -k2,2n $@_tmp | gzip > $@
	rm $@_tmp

$(intersect_with_trf): $(merged_unfiltered_probes)
	bedtools intersect -a $< -b $(trf) | gzip > $@

$(merged_final_probes): $(final_coordinates)
	bedtools merge -i $< | gzip > $@

$(merged_unfiltered_probes): $(unfiltered_coordinates)
	bedtools merge -i $< | gzip > $@

$(probe_count): $(final_coordinates)
	printf "chr\tprobes\nall\t" > $@
	zcat $< | wc -l >> $@
	for i in $(chromosomes); do \
	    printf "chr$${i}\t" >> $@; \
	    zgrep -w ^$${i} $< | wc -l >> $@; \
	done

# get coordinates of unique regions that don't overlap any TRF
$(unique_regions): $(alignability_filter) $(trf)
	bedtools subtract -a $(alignability_filter) -b $(trf) | gzip > $@

# create a local copy of the alignability filter
$(alignability_filter):
	read -p "MPI login: " USERNAME; \
	scp $${USERNAME}@bio23.eva.mpg.de:/mnt/454/HCNDCAM/Hengs_Alignability_Filter/hs37m_filt35_99.bed.gz $@

# download repetitive regions identified by Tandem Repeat Filter from UCSC
# and extract only columns with chromosome ID and start/end coordinates
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
