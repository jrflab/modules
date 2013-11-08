# run oncofuse
# hg19 only
include ~/share/modules/Makefile.inc

LOGDIR = log/oncofuse.$(NOW)

EXTRACT_COORDS = $(PERL) $(HOME)/share/scripts/extractCoordsFromDefuse.pl

ONCOFUSE_MEM = $(JAVA7) -Xmx$1 -jar $(HOME)/share/usr/oncofuse-v1.0.3/Oncofuse.jar

ONCOFUSE_TISSUE_TYPE ?= EPI

DEFUSE_RESULTS := defuse/tables/all.defuse_results.txt
CHIMSCAN_RESULTS := chimscan/tables/all.chimscan_results.txt
ifdef NORMAL_CHIMSCAN_RESULTS
CHIMSCAN_RESULTS := chimscan/tables/all.chimscan_results.nft.txt
endif
ifdef NORMAL_DEFUSE_RESULTS
DEFUSE_RESULTS := defuse/tables/all.defuse_results.nft.txt
endif


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all

all : oncofuse/defuse.merged_oncofuse_results.txt oncofuse/chimscan.merged_oncofuse_results.txt

oncofuse/defuse.coord.txt : $(DEFUSE_RESULTS)
	$(INIT) $(EXTRACT_COORDS) -t $(ONCOFUSE_TISSUE_TYPE) $< > $@ 2> $(LOG)

oncofuse/%.oncofuse_results.txt : oncofuse/%.coord.txt
	$(call LSCRIPT_MEM,8G,12G,"$(call ONCOFUSE_MEM,7G) $< coord $(ONCOFUSE_TISSUE_TYPE) $@")

.SECONDEXPANSION: 
oncofuse/%.merged_oncofuse_results.txt : $$*/tables/all.$$*_results.nft.txt oncofuse/%.oncofuse_results.txt
	head -1 $< | sed 's/^/RowID\t/' > $<.tmp && awk 'BEGIN {OFS = "\t" } NR > 1 { print NR-1, $$0 }' $< >> $<.tmp ;\
		$(RSCRIPT) $(MERGE) -X --byColX 1 --byColY 1 -H $<.tmp $(word 2,$^) > $@

oncofuse/chimscan.coord.txt : $(CHIMSCAN_RESULTS)
	$(INIT) cut -f 2,4,5,6 $< | awk 'BEGIN { OFS = "\t" } { print $$0, "$(ONCOFUSE_TISSUE_TYPE)" }' | sed '1d' > $@

