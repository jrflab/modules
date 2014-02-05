# run absCNseq on varscan segmentation data
include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

ABS_CN_SEQ = $(RSCRIPT) $(HOME)/share/scripts/absCNseq.R

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all

all : $(foreach pair,$(SAMPLE_PAIRS),absCN/$(pair).absCN.txt)

absCN/%.absCN.txt absCN/%.absSNV.txt : varscan/segment/%.varscan2copynumber.txt tables/%.mutect.som_ad_ft.pass.dbsnp.nsfp.chasm.fathmm.eff.opl_tab.txt
	$(call LSCRIPT_MEM,4G,6G,"$(ABS_CN_SEQ) --outPrefix absCN/$* $^")
