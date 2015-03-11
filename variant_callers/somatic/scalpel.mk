# scalpel variant detection
# vim: set ft=make :

include modules/Makefile.inc
include modules/variant_callers/gatk.inc
include modules/variant_callers/somatic/scalpel.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR = log/scalpel.$(NOW)

..DUMMY := $(shell mkdir -p version; echo "$(SCALPEL) $(SCALPEL_OPTS) > version/scalpel.txt")

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all scalpel_vcfs scalpel_tables

scalpel : scalpel_vcfs scalpel_tables

scalpel_vcfs : $(call VCFS,scalpel_indels)
scalpel_tables : $(call TABLES,scalpel_indels)

ifdef BED_FILES
define scalpel-bed-tumor-normal
scalpel/$2_$3/$1/somatic.$(SCALPEL_MIN_COV)x.indel.annovar : bam/$2.bam.md5 bam/$3.bam.md5
	$$(call LSCRIPT_NAMED_PARALLEL_MEM,$2_$3_$1_scalpel,2,4G,7G,"$$(SCALPEL) --somatic --numprocs 2 --tumor $$(<M) --normal $$(<<M) $$(SCALPEL_OPTS) --bed $$(BED_DIR)/$1 --dir $$(@D)")
endef
$(foreach bed,$(BED_FILES),\
	$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call scalpel-bed-tumor-normal,$(bed),$(tumor.$(pair)),$(normal.$(pair))))))

define merge-scalpel-tumor-normal
scalpel/$1_$2/somatic.$(SCALPEL_MIN_COV)x.indel.annovar : $$(foreach bed,$$(BED_FILES),scalpel/$1_$2/$$(bed)/somatic.$(SCALPEL_MIN_COV)x.indel.annovar)
	$$(INIT) grep '^#chr' $$< > $$@ && sed '/^#/d' $$^ >> $$@
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call merge-scalpel-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

else # dont split across exome bed files

define scalpel-tumor-normal
scalpel/$1_$2/somatic.$(SCALPEL_MIN_COV)x.indel.annovar : bam/$1.dcov.bam.md5 bam/$2.dcov.bam.md5
	$$(call LSCRIPT_NAMED_PARALLEL_MEM,$1_$2_scalpel,8,1G,4G,"$$(SCALPEL) --somatic --numprocs 8 --tumor $$(<M) --normal $$(<<M) $$(SCALPEL_OPTS) --dir $$(@D)")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call scalpel-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif

define scalpel2vcf-tumor-normal
vcf/$1_$2.scalpel.vcf : scalpel/$1_$2/somatic.$(SCALPEL_MIN_COV)x.indel.annovar
	$$(INIT) $$(SCALPEL2VCF) -f $$(REF_FASTA) -t $1 -n $2 < $$< > $$@ 2> $$(LOG)
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call scalpel2vcf-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include modules/vcf_tools/vcftools.mk
include modules/bam_tools/processBam.mk
