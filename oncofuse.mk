# run oncofuse
# hg19 only
include ~/share/modules/Makefile.inc

LOGDIR = log/oncofuse.$(NOW)

EXTRACT_COORDS = $(PERL) $(HOME)/share/scripts/extractCoordsFromDefuse.pl

ONCOFUSE_MEM = $(JAVA7) -Xmx$1 -jar $(HOME)/share/usr/oncofuse-v1.0.3/Oncofuse.jar

ONCOFUSE_TISSUE_TYPE ?= EPI

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all

all : $(foreach sample,$(SAMPLES),oncofuse/tables/$(sample).oncofuse_results.txt)

oncofuse/coord.txt : defuse/tables/all.defuse_results.txt
	$(INIT) $(EXTRACT_COORDS) -t $(ONCOFUSE_TISSUE_TYPE) $< > $@ 2> $(LOG)

oncofuse/oncofuse_results.txt : oncofuse/coord.txt
	$(call LSCRIPT_MEM,8G,12G,"$(call ONCOFUSE_MEM,7G) $< coord $(ONCOFUSE_TISSUE_TYPE) $@")

oncofuse/defuse_oncofuse_results.txt : defuse/tables/all.defuse_results.txt oncofuse/oncofuse_results.txt
	head -1 $< | sed 's/^/RowID/' > $<.tmp && awk 'NR > 1 { print NR-1, $$0 }' $< >> $<.tmp \
		$(RSCRIPT) $(MERGE) $<.tmp $(word 2,$^) 
