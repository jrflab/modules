# vim: set ft=make :
# amplicon qc using bams and gatk vcf results

include modules/Makefile.inc

INTERVAL_FILE ?= intervals.bed

VPATH ?= bam

TEQC = modules/qc/TEQC.R
INTERVAL_BAM_QC = modules/qc/intervalBamQC.R
VARIANT_EVAL_REPORT = modules/qc/variantEvalGatkReport.R

LOGDIR ?= log/interval_qc.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: rdata all

all : interval_qc/coverage/index.html # amplicon_qc/variant_eval/index.html
rdata : $(foreach sample,$(SAMPLES),amplicon_qc/rdata/$(sample).Rdata)

# load each bam file into R and create R data files
interval_qc/rdata/%.Rdata : %.bam
	$(call RUN,-s 8G -m 18G,"$(RSCRIPT) $(TEQC) --outFile $@ --ref $(REF) $< $(INTERVAL_FILE)")

# GATK variant eval for each sample variants.vcf file
# stratified by intervals and filter
interval_qc/variantEval.grp : $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).variants.vcf)
	$(call RUN,-n 4 -s 1G -m 1.5G,"$(call GATK_MEM,4G) -T VariantEval -nt 4 -o $@ -R $(REF_FASTA) --stratIntervals $(INTERVAL_FILE) --dbsnp $(DBSNP) $(foreach vcf,$^, --eval:$(call strip-suffix,$(notdir $(vcf))) $(vcf)) -ST IntervalStratification -ST Filter")

# Create amplicon coverage plots using R script AMPLICON_BAM_QC
interval_qc/coverage/index.html : $(foreach sample,$(SAMPLES),interval_qc/rdata/$(sample).Rdata)
	$(call RUN,-s 2G -m 4G,"$(RSCRIPT) $(INTERVAL_BAM_QC) --outDir $(@D) $^")

# Create variant evaluation plots using R script VARIANT_EVAL_REPORT
interval_qc/variant_eval/index.html : interval_qc/variantEval.grp
	$(call RUN,-s 2G -m 4G,"$(RSCRIPT) $(VARIANT_EVAL_REPORT) --outDir $(@D) $< &> $(LOG)")
