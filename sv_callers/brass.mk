# run lumpy

include modules/Makefile.inc

LOGDIR = log/brass.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:

BAM_STATS = $(HOME)/share/usr/bin/bam_stats
BRASS = $(HOME)/share/usr/bin/brass.pl
HIGH_DEPTH_BED = $(HOME)/share/reference/extremedepth.bed
BRASS_REPEATS = $(HOME)/share/reference/brassRepeats.bed.gz
GENOME_CACHE = $(HOME)/share/reference/Homo_sapiens.GRCh37.74.vagrent.cache.gz
BRASS_PROTOCOL = WGS # WGS|WXS|RNA
BRASS_OPTS = -g $(REF_FASTA) -s HUMAN -as 37 -pr $(BRASS_PROTOCOL) -gc $(GENOME_CACHE) -d $(HIGH_DEPTH_BED) -r $(BRASS_REPEATS) 

PHONY += brass
brass : $(foreach pair,$(SAMPLE_PAIRS),brass/$(pair).brass_timestamp)

bam/%.bam.bas : bam/%.bam
	$(call LSCRIPT_MEM,6G,7G,"$(BAM_STATS) -i $< -o $@ -r $(REF_FASTA).fai")

brass/genome_window.bed :
	$(call LSCRIPT_MEM,2G,3G,"$(BEDTOOLS) makewindows -g $(REF_FASTA) -w 10000 > $@")

define brass-tumor-normal
brass/$1_$2.brass_input_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bas bam/$2.bam.bas
	$$(call LSCRIPT_PARALLEL_MEM,2,6G,10G,"$$(BRASS) $$(BRASS_OPTS) -t $$< -n $$(<<) -tn $1 -nn $2 -c 2 -p input -o brass/$1_$2 && touch $$@")

brass/$1_$2.brass_group_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bas bam/$2.bam.bas brass/$1_$2.brass_input_timestamp
	$$(call LSCRIPT_MEM,7G,9G,"$$(BRASS) $$(BRASS_OPTS) -t $$< -n $$(<<) -tn $1 -nn $2 -o brass/$1_$2 -p group && touch $$@")

brass/$1_$2.brass_filter_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bas bam/$2.bam.bas brass/$1_$2.brass_group_timestamp
	$$(call LSCRIPT_MEM,7G,9G,"$$(BRASS) $$(BRASS_OPTS) -t $$< -n $$(<<) -tn $1 -nn $2 -o brass/$1_$2 -p filter && touch $$@")

brass/$1_$2.brass_assemble_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bas bam/$2.bam.bas brass/$1_$2.brass_filter_timestamp
	$$(call LSCRIPT_PARALLEL_MEM,8,2G,3G,"$$(BRASS) $$(BRASS_OPTS) -t $$< -n $$(<<) -tn $1 -nn $2 -o brass/$1_$2 -c 8 -p assemble && touch $$@")

brass/$1_$2.brass_grass_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bas bam/$2.bam.bas brass/$1_$2.brass_assemble_timestamp
	$$(call LSCRIPT_MEM,7G,9G,"$$(BRASS) $$(BRASS_OPTS) -t $$< -n $$(<<) -tn $1 -nn $2 -o brass/$1_$2 -p grass && touch $$@")

brass/$1_$2.brass_tabix_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bas bam/$2.bam.bas brass/$1_$2.brass_grass_timestamp
	$$(call LSCRIPT_MEM,2G,3G,"$$(BRASS) $$(BRASS_OPTS) -t $$< -n $$(<<) -tn $1 -nn $2 -o brass/$1_$2 -p tabix && touch $$@")

brass/$1_$2.brass_timestamp : brass/$1_$2.brass_tabix_timestamp
	$$(INIT) touch $$@
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call brass-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

.PHONY: $(PHONY)
