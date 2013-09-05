# run oncofuse
# hg19 only
include ~/share/modules/Makefile.inc

LOGDIR = log/oncofuse.$(NOW)

EXTRACT_COORDS = $(PERL) $(HOME)/share/scripts/extractCoordsFromDefuse.pl

ONCOFUSE_MEM = $(JAVA) -Xmx$1 -jar $(JARDIR)/Oncofuse.jar

ONCOFUSE_TISSUE_TYPE ?= AVG


oncofuse/coord/%.coord : defuse/tables/%.defuse_results.txt
	$(INIT) $(EXTRACT_COORDS) $< > $@ 2> $(LOG)

oncofuse/tables/%.oncufuse_results.txt : oncofuse/coord/%.coord
	$(LSCRIPT_MEM,4G,6G,"$(call ONCOFUSE_MEM,4G) $< coord $(ONCOFUSE_TISSUE_TYPE) $@")

