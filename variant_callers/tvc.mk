# use torrent variant caller on ion torrent bams

LOGDIR ?= log/tvc.$(NOW)

include modules/Makefile.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc

VPATH ?= bam

TVC_OPTS ?= -r $(REF_FASTA) $(if $(TARGETS_FILE),-t $(TARGETS_FILE)) 

.DELETE_ON_ERROR:
.SECONDARY: 

PHONY += tvc tvc_vcfs

tvc : tvc_vcfs

tvc_vcfs : $(foreach sample,$(SAMPLES),vcf/$(sample).tvc_snps_indels.vcf)

vcf/%.tvc_snps_indels.vcf : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_PARALLEL_MEM,4,1G,2G,"$(TVC) $(TVC_OPTS) -n 4 -o $(@) -b $<")

.PHONY: $(PHONY)

