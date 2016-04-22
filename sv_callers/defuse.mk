# This module is used for running defuse
# input: $(SAMPLES) 
# Options: BAM_PHRED64 = true/false
# Authors: Fong Chun Chan <fongchunchan@gmail.com>
#
include modules/Makefile.inc

DEFUSE_CONFIG_FILE = $(HOME)/share/usr/defuse-0.6.1/scripts/config.txt
#DEFUSE_CONFIG_FILE = /opt/common/defuse/defuse-0.6.1/scripts/config.txt

DEFUSE_FILTER = $(PERL) modules/sv_callers/filterDefuse.pl
DEFUSE_NORMAL_FILTER = $(PERL) modules/sv_callers/normalFilterDefuse.pl

RECURRENT_FUSIONS = $(RSCRIPT) modules/sv_callers/recurrentFusions.R
#EXTRACT_COORDS = $(PERL) modules/sv_callers/extractCoordsFromDefuse.pl
DEFUSE_ONCOFUSE = $(RSCRIPT) modules/sv_callers/defuseOncofuse.R
DEFUSE_ONCOFUSE_OPTS = --oncofuseJar $(ONCOFUSE_JAR) --oncofuseTissueType $(ONCOFUSE_TISSUE_TYPE) --java $(JAVA_BIN) 
ONCOFUSE_TISSUE_TYPE ?= EPI

LOGDIR = log/defuse.$(NOW)

# Runs defuse locally on the same node
LOCAL ?= FALSE

# Only applies if LOCAL is set to TRUE
NUM_CORES ?= 2

ifeq ($(LOCAL),true)
	DEFUSE_OPTS = -p $(NUM_CORES)
else
	DEFUSE_OPTS = -s sge -p 10
endif

.PHONY : all tables


#all : $(foreach sample,$(SAMPLES),defuse/$(sample).defuse_timestamp)
ALL = $(foreach sample,$(SAMPLES),defuse/tables/$(sample).defuse.txt)
ifdef NORMAL_DEFUSE_RESULTS
ALLTABLE = defuse/alltables/all.defuse.nft.oncofuse.txt
ALLTABLE += defuse/alltables/all.defuse_ft.nft.oncofusetxt
ALL += defuse/recur_tables/recurFusions.defuse.nft.gene.txt
ALL += defuse/recur_tables/recurFusions.defuse_ft.nft.gene.txt
else
ALLTABLE = defuse/alltables/all.defuse.oncofuse.txt
ALLTABLE += defuse/alltables/all.defuse_ft.oncofuse.txt
ALL += defuse/recur_tables/recurFusions.defuse.gene.txt
ALL += defuse/recur_tables/recurFusions.defuse_ft.gene.txt
endif
all : $(ALLTABLE) $(ALL)


defuse/tables/%.defuse.txt defuse/tables/%.defuse_ft.txt : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(INIT) $(DEFUSE) -c $(DEFUSE_CONFIG_FILE) -1 $(word 1,$(^M)) -2 $(word 2,$(^M)) -o defuse/$* $(DEFUSE_OPTS) &> $(LOG) && \
	$(DEFUSE_FILTER) defuse/$*/results.filtered.tsv > defuse/tables/$*.defuse_ft.txt 2>> $(LOG) && \
	$(DEFUSE_FILTER) defuse/$*/results.tsv > defuse/tables/$*.defuse.txt 2>> $(LOG) \
	&& $(RMR) defuse/$*

defuse/alltables/all.%.txt : $(foreach sample,$(SAMPLES),defuse/tables/$(sample).%.txt)
	$(INIT) head -1 $< > $@ && for x in $^; do sed '1d' $$x >> $@; done

defuse/alltables/%.nft.txt : defuse/alltables/%.txt
	$(INIT) $(DEFUSE_NORMAL_FILTER) -w 1000 $(NORMAL_DEFUSE_RESULTS) $< > $@

defuse/recur_tables/recurFusions.%.gene.txt : defuse/alltables/all.%.txt
	$(INIT) $(RECURRENT_FUSIONS) --geneCol1 upstream_gene --geneCol2 downstream_gene --sampleCol library_name --outPrefix $(@D)/recurFusions.$* $< 

defuse/alltables/%.coord.txt : defuse/alltables/%.txt
	$(INIT) $(EXTRACT_COORDS) -t $(ONCOFUSE_TISSUE_TYPE) $< > $@ 2> $(LOG)

defuse/alltables/%.oncofuse.txt : defuse/alltables/%.txt
	$(call LSCRIPT_CHECK_MEM,7G,8G,"$(DEFUSE_ONCOFUSE) --outPrefix $(@D)/$* $(DEFUSE_ONCOFUSE_OPTS) $<")

include modules/fastq_tools/fastq.mk
