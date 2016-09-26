include modules/Makefile.inc
include modules/config.inc

LOGDIR = log/tseq_workflow.$(NOW)

TSEQ_WORKFLOW += bwamem
TSEQ_WORKFLOW += mutect
TSEQ_WORKFLOW += strelka
TSEQ_WORKFLOW += varscan
TSEQ_WORKFLOW += facets
TSEQ_WORKFLOW += fastqc
TSEQ_WORKFLOW += bam_interval_metrics
TSEQ_WORKFLOW += mutation_summary
TSEQ_WORKFLOW += strelka_varscan_merge
TSEQ_WORKFLOW += gatk

PHONY += tseq_workflow
tseq_workflow : $(TSEQ_WORKFLOW)

include modules/aligners/bwamemAligner.mk
include modules/variant_callers/gatkVariantCaller.mk
include modules/variant_callers/somatic/mutect.mk
include modules/variant_callers/somatic/strelka.mk
include modules/variant_callers/somatic/varscanTN.mk
include modules/variant_callers/somatic/mergeStrelkaVarscanIndels.mk
include modules/copy_number/facets.mk
include modules/summary/mutationSummary.mk
include modules/qc/fastqc.mk
include modules/qc/bamIntervalMetrics.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)

