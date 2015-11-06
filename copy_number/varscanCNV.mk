# Run VarScan to detect copynumber
# Detect copy number
##### DEFAULTS ######

LOGDIR = log/varscan.$(NOW)

##### MAKE INCLUDES #####
include modules/Makefile.inc

SEGMENTCNV = modules/copy_number/segmentVarscanCNV.R
CGHCALL = modules/copy_number/cghCall.R

SEG_SD ?= 2
SEG_SMOOTH ?= 10
SEG_ALPHA ?= 0.000001

MULTIPARAM_SEGMENT ?= false
SEG_SDS ?= 1 2 3
SEG_SMOOTHS ?= 10 5
SEG_ALPHAS ?= 0.01 0.0001 0.000001 0.0000000001

VARSCAN_GENE_CN = $(RSCRIPT) modules/copy_number/varscanCNVGeneCN.R
VARSCAN_GENE_CN_OPTS = $(if $(GENES_FILE),--genesFile $(GENES_FILE))

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: all copycalls segments cnv

all : copycalls segments geneCN #cghcalls

CGHCALLS := $(foreach pair,$(SAMPLE_PAIRS),varscan/segment/$(pair).cgh_call.txt)
ifeq ($(MULTIPARAM_SEGMENT),true) 
CGHCALLS += $(foreach pair,$(SAMPLE_PAIRS),\
				$(foreach sd,$(SEG_SDS),\
					$(foreach alpha,$(SEG_ALPHAS),\
						$(foreach smooth,$(SEG_SMOOTHS),\
							varscan/segment_sd$(sd)_alpha$(alpha)_smooth$(smooth)/$(pair).cgh_call.txt))))
endif

segments : $(foreach pair,$(SAMPLE_PAIRS),varscan/segment/$(pair).segment.Rdata)
copycalls : $(foreach pair,$(SAMPLE_PAIRS),varscan/copycall/$(pair).copycall)
cghcalls : $(CGHCALLS)
geneCN : varscan/segment/geneCN.txt

define varscan-copynum-tumor-normal
varscan/copynum/$1_$2.copynumber :  bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_CHECK_MEM,9G,12G,"$$(SAMTOOLS) mpileup -q 1 -l $$(TARGETS_FILE) -f $$(REF_FASTA) $$(word 2,$$^) $$< | awk 'NF == 9 { print }' |  $$(VARSCAN) copynumber - $$(basename $$@) --mpileup 1")
endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call varscan-copynum-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))

varscan/copycall/%.copycall : varscan/copynum/%.copynumber
	$(call LSCRIPT_CHECK_MEM,9G,12G,"n=`awk '{ total += \$$7 } END { print total / NR }' $<`; \
	if [ \$$(bc <<< \"\$$n > 0\") -eq 1 ]; then \
		recenter_opt=\"--recenter-up \$$n\"; \
	else \
		n=\$$(bc <<< \"\$$n*-1\"); \
		recenter_opt=\"--recenter-down \$$n\"; \
	fi; \
	$(VARSCAN) copyCaller $< --output-file $@ \$$recenter_opt")

varscan/segment/%.segment.Rdata : varscan/copycall/%.copycall
	$(call LSCRIPT_CHECK_MEM,4G,6G,"$(RSCRIPT) $(SEGMENTCNV) --alpha $(SEG_ALPHA) --smoothRegion $(SEG_SMOOTH) --undoSD $(SEG_SD) --centromereFile=$(CENTROMERE_TABLE2) --prefix=$(@D)/$* $<")

varscan/segment/geneCN.txt : $(foreach pair,$(SAMPLE_PAIRS),varscan/segment/$(pair).collapsed_seg.txt)
	$(call LSCRIPT_CHECK_MEM,4G,6G,"$(VARSCAN_GENE_CN) $(VARSCAN_GENE_CN_OPTS) --outFile $@ $^")	

varscan/segment/%.cgh_call.txt : varscan/segment/%.segment.Rdata
	$(call LSCRIPT_CHECK_MEM,4G,6G,"$(RSCRIPT) $(CGHCALL) --centromereFile=$(CENTROMERE_TABLE2) --prefix=$(@D)/$* $<")

define varscan-segment-sd-alpha-smooth
varscan/segment_sd$1_alpha$2_smooth$3/%.segment.Rdata : varscan/copycall/%.copycall
	$$(call LSCRIPT_CHECK_NAMED_MEM,$$*_seg_$1_$2_$3,4G,6G,"$$(RSCRIPT) $$(SEGMENTCNV) --undoSD $1 --alpha $2 --smoothRegion $3 --centromereFile=$$(CENTROMERE_TABLE2) --prefix=$$(@D)/$$* $$<")

varscan/segment_sd$1_alpha$2_smooth$3/%.cgh_call.txt : varscan/segment_sd$1_alpha$2_smooth$3/%.segment.Rdata
	$$(call LSCRIPT_CHECK_NAMED_MEM,$$*_cgh_$1_$2_$3,4G,6G,"$$(RSCRIPT) $$(CGHCALL) --centromereFile=$$(CENTROMERE_TABLE2) --prefix=$$(@D)/$$* $$<")
endef
$(foreach sd,$(SEG_SDS),\
	$(foreach alpha,$(SEG_ALPHAS),\
		$(foreach smooth,$(SEG_SMOOTHS),\
			$(eval $(call varscan-segment-sd-alpha-smooth,$(sd),$(alpha),$(smooth))))))

