include ~/share/modules/Makefile.inc

QUALIMAP_BAMQC_OPTS = -gd HUMAN 
QUALIMAP = unset DISPLAY; $(JAVA) -Xmx18G -classpath $(HOME)/share/usr/qualimap/qualimap.jar:$(HOME)/share/usr/qualimap/lib/* org.bioinfo.ngs.qc.qualimap.main.NgsSmartMain 


SAMPLE_FILE = samples.txt
SAMPLES = $(shell cat $(SAMPLE_FILE))

LOGDIR = log/qualimap.$(NOW)


ifeq ($(EXOME),true)
TARGETS_FILE = $(HOME)/share/reference/SureSelect_50MB_S02972011_Regions_nochr.bed
endif

ifdef TARGETS_FILE
QUALIMAP_BAMQC_OPTS += -gff $(TARGETS_FILE) -os
endif

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : all

all : $(foreach sample,$(SAMPLES),qualimap/$(sample)_bamqc.timestamp)

qualimap/%_bamqc.timestamp : bam/%.clean.bam
	$(call INIT_PARALLEL_MEM,10,2G,2.5G) $(QUALIMAP) bamqc $(QUALIMAP_BAMQC_OPTS) -bam $< -nr 6 -nt 10 -outdir qualimap/$*_bamqc &> $(LOG) && touch $@


include ~/share/modules/processBam.mk
