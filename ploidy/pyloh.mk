# vim: set ft=make :
# pyloh loss of heterozygosity
# uses exomecnv segment calls

include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR := log/pyloh.$(NOW)

PYLOH = PYTHONPATH=$(HOME)/usr/lib/python:$(HOME)/usr/lib/python2.7 $(HOME)/usr/bin/python $(HOME)/usr/bin/PyLOH.py

BICSEQ = $(HOME)/usr/BICseq_1.1.2/BIC-seq/BICseq/BIC-seq
SAMTOOLS_UNIQ = $(HOME)/usr/BICseq_1.1.2/SAMgetUnique/samtools-0.1.7a_getUnique-0.1.1/samtools

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all run_models heatmaps

all : run_models
run_models : $(foreach pair,$(SAMPLE_PAIRS),pyloh/$(pair).run_timestamp)
heatmaps : $(foreach pair,$(SAMPLE_PAIRS),pyloh/$(pair).PyLOH.heatmap.pkl)

pyloh/%.readpos_timestamp : bam/%.bam
	$(call LSCRIPT_MEM,6G,8G,"mkdir -p pyloh/$*_readpos && $(SAMTOOLS_UNIQ) view -U BWA$(,)pyloh/$*_readpos/$(,)N$(,)N $< && touch $@")

define bicseq-tumor-normal-chr
pyloh/$1_$2_chr_segments/$3.bicseq.txt : pyloh/$1.readpos_timestamp pyloh/$2.readpos_timestamp
	$$(call LSCRIPT_MEM,6G,8G,"$$(BICSEQ) -l 10 pyloh/$1_readpos/$3.seq pyloh/$2_readpos/$3.seq > $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach chr,$(CHROMOSOMES),\
	$(eval $(call bicseq-tumor-normal-chr,$(tumor.$(pair)),$(normal.$(pair)),$(chr)))))

define pyloh-tumor-normal
pyloh/$1_$2.segments.bed : $$(foreach chr,$$(CHROMOSOMES),pyloh/$1_$2_chr_segments/$$(chr).bicseq.txt)
	$$(INIT) for chr in $$(CHROMOSOMES); do \
		awk -v chr=$$$$chr 'BEGIN { OFS = "\t" } $$$$5 - $$$$4 > 1000000 { print chr, $$$$4, $$$$5 }' pyloh/$1_$2_chr_segments/$$$${chr}.bicseq.txt >> $$@; \
	done

pyloh/$1_$2.preprocess_timestamp : bam/$1.bam bam/$2.bam pyloh/$1_$2.segments.bed
	$$(call LSCRIPT_PARALLEL_MEM,6,1G,2G,"$$(PYLOH) preprocess $$(REF_FASTA) $$< $$(<<) pyloh/$1_$2 --segments_bed $$(<<<) --min_depth 20 --min_base_qual 10 --min_map_qual 10 --process_num 6 && touch $$@")

pyloh/$1_$2.run_timestamp : pyloh/$1_$2.preprocess_timestamp
	$$(call LSCRIPT_MEM,7G,10G,"$$(PYLOH) run_model pyloh/$1_$2 --allelenumber_max 2 --max_iters 100 --stop_value 1e-7 && touch $$@")

pyloh/$1_$2.PyLOH.heatmap.pkl : pyloh/$1_$2.preprocess_timestamp
	$$(call LSCRIPT_MEM,6G,8G,"$$(PYLOH) BAF_heatmap pyloh/$1_$2")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call pyloh-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
