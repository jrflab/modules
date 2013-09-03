include ~/share/modules/Makefile.inc

QUALIMAP_BAMQC_OPTS = -gd HUMAN 
QUALIMAP = unset DISPLAY; $(JAVA) -Xmx16G -classpath $(HOME)/share/usr/qualimap/qualimap.jar:$(HOME)/share/usr/qualimap/lib/* org.bioinfo.ngs.qc.qualimap.main.NgsSmartMain 


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

qualimap/%_bamqc.timestamp : bam/%.bam
	$(call LSCRIPT_PARALLEL_MEM,4,4.5G,5G,"$(QUALIMAP) bamqc $(QUALIMAP_BAMQC_OPTS) -bam $< -nr 6 -nt 8 -outdir qualimap/$*_bamqc &> $(LOG) && touch $@")


include ~/share/modules/processBamMD5.mk
