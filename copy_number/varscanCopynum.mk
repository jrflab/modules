# Run VarScan to detect copynumber
# Detect copy number
##### DEFAULTS ######

REF ?= hg19
LOGDIR = log/varscan.$(NOW)

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc

SEGMENTCNV = $(HOME)/share/scripts/segmentCNV2.R
CGHCALL = $(HOME)/share/scripts/cghCall.R

SEG_SDS = 1 2 3
SEG_SMOOTHS = 10 100
SEG_ALPHAS = 0.05 0.01 0.0001 0.000001

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: all copycalls segments cnv

all : copycalls segments cghcalls


segments : $(foreach pair,$(SAMPLE_PAIRS),varscan/segment/$(pair).segment.Rdata)
copycalls : $(foreach pair,$(SAMPLE_PAIRS),varscan/copycall/$(pair).copycall)
cghcalls : $(foreach pair,$(SAMPLE_PAIRS),\
				$(foreach sd,$(SEG_SDS),\
					$(foreach alpha,$(SEG_ALPHAS),\
						$(foreach smooth,$(SEG_SMOOTHS),\
							varscan/segment_sd$(sd)_alpha$(alpha)_smooth$(smooth)/$(pair).cgh_call.txt))))

define varscan-copynum-tumor-normal
varscan/copynum/$1_$2.copynumber :  bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_MEM,9G,12G,"$$(SAMTOOLS) mpileup -q 1 -f $$(REF_FASTA) $$(word 2,$$^) $$< | awk 'NF == 9 { print }' |  $$(VARSCAN) copynumber - $$(basename $$@) --mpileup 1")
endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call varscan-copynum-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))

varscan/copycall/%.copycall : varscan/copynum/%.copynumber
	$(call LSCRIPT_MEM,9G,12G,"n=`awk '{ total += \$$7 } END { print total / NR }' $<`; \
	if [ \$$(bc <<< \"\$$n > 0\") -eq 1 ]; then \
		recenter_opt=\"--recenter-up \$$n\"; \
	else \
		n=\$$(bc <<< \"\$$n*-1\"); \
		recenter_opt=\"--recenter-down \$$n\"; \
	fi; \
	$(VARSCAN) copyCaller $< --output-file $@ \$$recenter_opt")


define varscan-segment-sd-alpha-smooth
varscan/segment_sd$1_alpha$2_smooth$3/%.segment.Rdata : varscan/copycall/%.copycall
	$$(call LSCRIPT_MEM,4G,6G,"$$(RSCRIPT) $$(SEGMENTCNV) --undoSD $1 --alpha $2 --smoothRegion $3 --centromereFile=$$(CENTROMERE_TABLE2) --prefix=$$(@D)/$$* $$<")

varscan/segment_sd$1_alpha$2_smooth$3/%.cgh_call.txt : varscan/segment_sd$1_alpha$2_smooth$3/%.segment.Rdata
	$$(call LSCRIPT_MEM,4G,6G,"$$(RSCRIPT) $$(CGHCALL) --centromereFile=$$(CENTROMERE_TABLE2) --prefix=$$(@D)/$$* $$<")
endef
$(foreach sd,$(SEG_SDS),\
	$(foreach alpha,$(SEG_ALPHAS),\
		$(foreach smooth,$(SEG_SMOOTHS),\
			$(eval $(call varscan-segment-sd-alpha-smooth,$(sd),$(alpha),$(smooth))))))

