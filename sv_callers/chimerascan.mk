# Chimerascan

##### MAKE INCLUDES #####
include modules/Makefile.inc

LOGDIR = log/chimscan.$(NOW)

CHIMSCAN_PYTHONPATH := /home/limr/share/usr/lib/python:/home/limr/share/usr/lib/python2.7
CHIMSCAN_PYTHON := $(HOME)/share/usr/bin/python
CHIMSCAN_INDEX := $(HOME)/share/reference/chimerascan_index
CHIMSCAN_NORMAL_FILTER = modules/sv_callers/normalFilterChimerascan.pl

RECURRENT_FUSIONS = $(RSCRIPT) modules/sv_callers/recurrentFusions.R

ONCOFUSE_MEM = $(JAVA) -Xmx$1 -jar $(ONCOFUSE_JAR)
ONCOFUSE_TISSUE_TYPE ?= EPI


USE_BIG_MEM ?= false
ifeq ($(USE_BIG_MEM),true)
CHIMSCAN_MEM := 240G
CHIMSCAN_CORES := 4
LAUNCH_CHIMSCAN = $(LSCRIPT_PARALLEL_MEM_BIG)
else
CHIMSCAN_MEM := 20G
CHIMSCAN_CORES := 4
LAUNCH_CHIMSCAN = $(LSCRIPT_PARALLEL_MEM)
endif

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all 

ALL = $(foreach sample,$(SAMPLES),chimscan/bedpe/$(sample).chimscan.bedpe)
ifdef NORMAL_CHIMSCAN_RESULTS
ALL += $(foreach sample,$(SAMPLES),chimscan/bedpe/$(sample).chimscan.nft.bedpe)
ALLTABLE = chimscan/alltables/all.chimscan.nft.oncofuse.merged.txt
ALL += chimscan/recur_tables/recurFusions.chimscan.nft.gene.txt
else 
ALLTABLE = chimscan/alltables/all.chimscan.oncofuse.merged.txt
ALL += chimscan/recur_tables/recurFusions.chimscan.gene.txt
endif
ALL += $(ALLTABLE)

all : $(ALL)

CHIMERASCAN = PYTHONPATH=$(CHIMSCAN_PYTHONPATH) $(CHIMSCAN_PYTHON) /home/limr/share/usr/lib/python/chimerascan/chimerascan_run.py
CHIMERASCAN_OPTS ?= 

#chimerascan/tables/all.chimscan_results.txt : $(foreach sample,$(SAMPLES),chimerascan/$(sample).chimscan_timestamp)
#$(INIT) head -1 $(basename $<)/chimeras.bedpe > $@ && for x in $(addsuffix /chimeras.bedpe,$(basename $^)); do sed '1d' $$x >> $@; done

chimscan/bedpe/%.chimscan.bedpe : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(call LAUNCH_CHIMSCAN,$(CHIMSCAN_CORES),$(CHIMSCAN_MEM),$(CHIMSCAN_MEM),"$(CHIMERASCAN) $(CHIMERASCAN_OPTS) -v --quals illumina -p $(CHIMSCAN_CORES) $(CHIMSCAN_INDEX) $^ chimscan/$* && cp -f chimscan/$*/chimeras.bedpe $@ && rm -r chimscan/$*")
# was 8 and 20
%.chimscan.nft.bedpe : %.chimscan.bedpe
	$(call LSCRIPT_MEM,2G,4G,"$(PERL) $(CHIMSCAN_NORMAL_FILTER) -w 1000 $(NORMAL_CHIMSCAN_RESULTS) $< > $@")

chimscan/alltables/all.chimscan%txt : $(foreach sample,$(SAMPLES),chimscan/bedpe/$(sample).chimscan%bedpe)
	$(INIT) head -1 $< | sed 's/^/Sample\t/; s/#//' > $@ && for i in $^; do sed "1d; s/^/$$(basename $${i%%.chimscan$*bedpe})\t/" $$i >> $@; done

chimscan/recur_tables/recurFusions.%.gene.txt : chimscan/alltables/all.%.txt
	$(INIT) $(RECURRENT_FUSIONS) --geneCol1 genes5p --geneCol2 genes3p --sampleCol Sample --outPrefix chimscan/recur_tables/recurFusions.$*  $< 

chimscan/alltables/all.chimscan%coord.txt : chimscan/alltables/all.chimscan%txt
	$(INIT) perl -lane 'if ($$. > 1) { $$coord5 = (($$F[9] eq "+")? $$F[3] + 1 : $$F[2] - 1) + 1; $$coord3 = (($$F[10] eq "+")? $$F[5] - 1 : $$F[6] + 1) + 1; print "$$F[1]\t$$coord5\t$$F[4]\t$$coord3\tEPI"; }' $< > $@

include modules/sv_callers/oncofuse.mk
