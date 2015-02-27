include ~/share/modules/Makefile.inc

LOGDIR = log/hmmcopy.$(NOW)

HMMCOPY_WINDOW_SIZE ?= 1000
HMMCOPY = $(RSCRIPT) $(HOME)/share/scripts/hmmCopy.R

READ_COUNTER = $(HOME)/share/usr/bin/readCounter
MAP_COUNTER = $(HOME)/share/usr/bin/mapCounter
GC_COUNTER = $(HOME)/share/usr/bin/gcCounter

MAP_BW = $(HOME)/share/references/genomes/wgEncodeCrgMapabilityAlign100mer.bigWig

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : hmmcopy

hmmcopy : $(foreach pair,$(SAMPLE_PAIRS),hmmcopy/results.w$(WINDOW_SIZE)_wig/$(tumor.$(pair)).$(normal.$(pair)).hmmcopy_seg.txt)

hmmcopy/wig/%.w$(HMMCOPY_WINDOW_SIZE).wig : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,6G,8G,"$(READ_COUNTER) -w $(HMMCOPY_WINDOW_SIZE) -c $(subst $( ),$(,),$(strip $(CHROMOSOMES))) $< > $@")

hmmcopy/wig/gc.w$(HMMCOPY_WINDOW_SIZE).wig :
	$(call LSCRIPT_MEM,6G,8G,"$(GC_COUNTER) -w $(HMMCOPY_WINDOW_SIZE) -c $(subst $( ),$(,),$(strip $(CHROMOSOMES))) $(REF_FASTA) > $@")

hmmcopy/wig/map.w$(HMMCOPY_WINDOW_SIZE).wig :
	$(call LSCRIPT_MEM,6G,8G,"$(MAP_COUNTER) -w $(HMMCOPY_WINDOW_SIZE) -c $(subst $( ),$(,),$(strip $(CHROMOSOMES))) $(MAP_BIGWIG) > $@")

define hmmcopy-tumor-normal
hmmcopy/results.w$$(HMMCOPY_WINDOW_SIZE)/$1_$2.hmmcopy_seg.txt : hmmcopy/wig/$1.w$$(HMMCOPY_WINDOW_SIZE).wig hmmcopy/wig/$2.w$$(HMMCOPY_WINDOW_SIZE).wig hmmcopy/wig/gc.w$$(HMMCOPY_WINDOW_SIZE).wig hmmcopy/wig/map.w$$(HMMCOPY_WINDOW_SIZE).wig 
	$$(call LSCRIPT_MEM,8G,12G,"$$(HMMCOPY) --normalWig $$(<<) --gcWig $$(<<<) --mapWig $$(<<<<) --outPrefix $$(@D)/$1_$2 $$<")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call hmmcopy-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
