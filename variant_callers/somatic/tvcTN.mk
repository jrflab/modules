# use torrent variant caller on ion torrent bams

LOGDIR ?= log/tvcTN.$(NOW)

include modules/Makefile.inc
include modules/variant_callers/gatk.inc

VPATH ?= bam

TVC_OPTS ?= -r $(REF_FASTA) $(if $(TARGETS_FILE),-t $(TARGETS_FILE)) 

.DELETE_ON_ERROR:
.SECONDARY: 

PHONY += tvc tvc_vcfs tvc_tables

tvc : tvc_vcfs tvc_tables

tvc_vcfs : $(call SOMATIC_VCFS,tvc_snps_indels)
tvc_tables : $(call SOMATIC_TABLES,tvc_snps_indels)

%.contig.vcf : %.vcf 
	$(INIT) awk '{print "##contig=<ID=" $$1 ",length=" $$2 ",assembly=$(REF)>"'} $(REF_FASTA).fai | $(BCFTOOLS2) annotate -h - $< > $@ 

%.vcf.gz : %.vcf
	$(call LSCRIPT,"bgzip -c $< > $@")

define tvc-tumor
vcf/$1.tvc_snps_indels.vcf : bam/$1.bam bam/$1.bam.bai
	$$(call LSCRIPT_PARALLEL_MEM,4,1G,2G,"$$(TVC) $$(TVC_OPTS) -b $$< -o $$@ -n 4")
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call tvc-tumor,$(tumor))))

define tvc-normal
vcf/$1.tumor_tvc_snps_indels.vcf : $$(foreach tumor,$2,vcf/$$(tumor).tvc_snps_indels.vcf.gz vcf/$$(tumor).tvc_snps_indels.vcf.gz.tbi)
	$$(INIT) $$(BCFTOOLS2) merge $$(filter %.vcf.gz,$$^) > $$@

vcf/$1.tvc_snps_indels.vcf : bam/$1.bam vcf/$1.tumor_tvc_snps_indels.vcf bam/$1.bam.bai
	$$(call LSCRIPT_PARALLEL_MEM,4,1G,2G,"$$(TVC) $$(TVC_OPTS) -c $$(<<) -b $$< -o $$@ -n 4")
endef
$(foreach normal,$(NORMAL_SAMPLES),$(eval $(call tvc-normal,$(normal),$(tumor.$(normal)))))

define tvc-tumor-normal
vcf/$1_$2.tvc_snps_indels.vcf : vcf/$1.tvc_snps_indels.vcf.gz vcf/$2.tvc_snps_indels.vcf.gz vcf/$1.tvc_snps_indels.vcf.gz.tbi vcf/$2.tvc_snps_indels.vcf.gz.tbi
	$$(INIT) $$(BCFTOOLS2) merge $$(filter %.vcf.gz,$$^) > $$@
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call tvc-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

.PHONY: $(PHONY)

include modules/vcf_tools/vcftools.mk
