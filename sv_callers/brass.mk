# run BRASS
# https://github.com/cancerit/BRASS

include modules/Makefile.inc

LOGDIR = log/brass.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:

BAM_STATS = $(HOME)/share/usr/bin/bam_stats
BRASS = $(HOME)/share/usr/bin/brass.pl
BRASS_HIGH_DEPTH_BED = $(HOME)/share/reference/extremedepth.bed
BRASS_REPEATS = $(HOME)/share/reference/brassRepeats.bed.gz
BRASS_GENOME_CACHE ?= $(HOME)/share/reference/Homo_sapiens.GRCh37.74.vagrent.cache.gz
BRASS_PROTOCOL ?= RNA # WGS|WXS|RNA
BRASS_OPTS = -g $(REF_FASTA) -s HUMAN -as 37 -pr $(BRASS_PROTOCOL) -gc $(BRASS_GENOME_CACHE) -d $(BRASS_HIGH_DEPTH_BED) -r $(BRASS_REPEATS) 

PHONY += brass
brass : $(foreach sample,$(SAMPLES),brass/$(sample).brass_timestamp)

brass/bas/%.bas : bam/%.bam
	$(call LSCRIPT_MEM,6G,7G,"$(BAM_STATS) -i $< -o $@ -r $(REF_FASTA).fai")

brass/genome_window.bed :
	$(call LSCRIPT_MEM,2G,3G,"$(BEDTOOLS) makewindows -g $(REF_FASTA) -w 10000 > $@")

brass/%.brass_input_timestamp : bam/%.bam brass/bas/%.bas
	$(call LSCRIPT_PARALLEL_MEM,2,6G,10G,"$(BRASS) $(BRASS_OPTS) -t $< -t $* -c 2 -p input -o brass/$* && touch $@")

brass/%.brass_group_timestamp : bam/%.bam brass/%.brass_input_timestamp
	$(call LSCRIPT_MEM,7G,9G,"$(BRASS) $(BRASS_OPTS) -t $< -tn $* -o brass/$* -p group && touch $@")

brass/%.brass_filter_timestamp : bam/%.bam brass/%.brass_group_timestamp
	$(call LSCRIPT_MEM,7G,9G,"$(BRASS) $(BRASS_OPTS) -t $< -tn $* -o brass/$* -p filter && touch $@")

brass/%.brass_assemble_timestamp : bam/%.bam brass/%.brass_filter_timestamp
	$(call LSCRIPT_PARALLEL_MEM,8,2G,3G,"$(BRASS) $(BRASS_OPTS) -t $< -tn $* -o brass/$* -c 8 -p assemble && touch $@")

brass/%.brass_grass_timestamp : bam/%.bam brass/%.brass_assemble_timestamp
	$(call LSCRIPT_MEM,7G,9G,"$(BRASS) $(BRASS_OPTS) -t $< -tn $* -o brass/$* -p grass && touch $@")

brass/%.brass_tabix_timestamp : bam/%.bam bam/%.bam.bas brass/%.brass_grass_timestamp
	$(call LSCRIPT_MEM,2G,3G,"$(BRASS) $(BRASS_OPTS) -t $< -tn $* -o brass/$* -p tabix && touch $@")

brass/%.brass_timestamp : brass/%.brass_tabix_timestamp
	$(INIT) touch $@

.PHONY: $(PHONY)
