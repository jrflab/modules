# run absCNseq on varscan segmentation data
include ~/share/modules/Makefile.inc
include ~/share/modules/variant_callers/gatk.inc

LOGDIR = log/absCN.$(NOW)
ABS_CN_SEQ = $(RSCRIPT) $(HOME)/share/scripts/absCNseq.R

# file containing positions: chr:start-stop
SNV_POS_FILE = snv_posn.intervals

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all

all : $(foreach pair,$(SAMPLE_PAIRS),absCN/$(pair).absCN.txt)

define abs-gatk-tumor-normal
absCN/$1_$2.gatk.vcf : bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_MEM,9G,12G,"$$(call GATK_MEM,8G) -T UnifiedGenotyper -o $$@ -I $$< -I $$(<<) -R $$(REF_FASTA)  --output_mode EMIT_ALL_SITES -L $$(SNV_POS_FILE)")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call abs-gatk-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

absCN/%.absCN.txt absCN/%.absSNV.txt : varscan/segment/%.varscan2copynumber.txt absCN/%.gatk.vcf
	$(call LSCRIPT_MEM,4G,6G,"$(ABS_CN_SEQ) --genome $(REF) --outPrefix absCN/$* $^")
