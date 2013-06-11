# vim: set ft=make :
# amplicon qc using bams and gatk vcf results

include ~/share/modules/Makefile.inc

SAMPLE_FILE ?= samples.txt
SAMPLES ?= $(shell cat $(SAMPLE_FILE))

INTERVAL_FILE ?= intervals.bed

VPATH ?= bam

TEQC = $(HOME)/share/scripts/TEQC.R
AMPLICON_BAM_QC = $(HOME)/share/scripts/ampliconBamQC.R
VARIANT_EVAL_REPORT = $(HOME)/share/scripts/variantEvalGatkReport.R

LOGDIR ?= log/amplicon_qc.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: rdata all

all : amplicon_qc/coverage/index.html amplicon_qc/variant_eval/index.html
rdata : $(foreach sample,$(SAMPLES),amplicon_qc/rdata/$(sample).Rdata)

# load each bam file into R and create R data files
amplicon_qc/rdata/%.Rdata : %.bam
	$(call INIT_MEM,4G,8G) $(RSCRIPT) $(TEQC) --outFile $@ --ref $(REF) $< $(INTERVAL_FILE) &> $(LOG)

# GATK variant eval for each sample variants.vcf file
# stratified by intervals and filter
amplicon_qc/variantEval.grp : $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).variants.vcf)
	$(call INIT_PARALLEL_MEM,4,1G,1.5G) $(call GATK_MEM,4G) -T VariantEval -nt 4 -o $@ -R $(REF_FASTA) --stratIntervals $(INTERVAL_FILE) --dbsnp $(DBSNP) $(foreach vcf,$^, --eval:$(call strip-suffix,$(notdir $(vcf))) $(vcf)) -ST IntervalStratification -ST Filter &> $(LOG)

# Create amplicon coverage plots using R script AMPLICON_BAM_QC
amplicon_qc/coverage/index.html : $(foreach sample,$(SAMPLES),amplicon_qc/rdata/$(sample).Rdata)
	$(call INIT_MEM,2G,4G) $(RSCRIPT) $(AMPLICON_BAM_QC) --outDir $(@D) $^ &> $(LOG)

# Create variant evaluation plots using R script VARIANT_EVAL_REPORT
amplicon_qc/variant_eval/index.html : amplicon_qc/variantEval.grp
	$(call INIT_MEM,2G,4G) $(RSCRIPT) $(VARIANT_EVAL_REPORT) --outDir $(@D) $< &> $(LOG)
