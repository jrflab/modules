# scalpel variant detection
# vim: set ft=make :

include modules/Makefile.inc
include modules/variant_callers/gatk.inc
include modules/variant_callers/somatic/scalpel.inc

LOGDIR = log/scalpel.$(NOW)

..DUMMY := $(shell mkdir -p version; echo "$(SCALPEL) $(SCALPEL_OPTS) > version/scalpel.txt")

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all scalpel_vcfs scalpel_mafs

scalpel : scalpel_vcfs scalpel_mafs


scalpel_vcfs : $(foreach pair,$(SAMPLE_PAIR),vcf/$(pair).scalpel_indels.vcf)
scalpel_mafs : $(foreach pair,$(SAMPLE_PAIR),maf/$(pair).scalpel_indels.maf)

ifdef BED_FILES
define scalpel-bed-tumor-normal
scalpel/$2_$3/$1/somatic.$(SCALPEL_MIN_COV)x.indel.annovar : bam/$2.bam bam/$3.bam
	$$(call LSCRIPT_NAMED_PARALLEL_MEM,$2_$3_$1_scalpel,2,4G,7G,"$$(SCALPEL) --somatic --numprocs 2 --tumor $$(<) --normal $$(<<) $$(SCALPEL_OPTS) --bed $$(BED_DIR)/$1 --dir $$(@D) && \
		if [! -e $$@ ]; then cp $$(@D)/main/$$(@F) $$@; fi")
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
scalpel/$1_$2/somatic.$(SCALPEL_MIN_COV)x.indel.annovar : bam/$1.dcov.bam bam/$2.dcov.bam
	$$(call LSCRIPT_NAMED_PARALLEL_MEM,$1_$2_scalpel,8,1G,4G,"$$(SCALPEL) --somatic --numprocs 8 --tumor $$(<M) --normal $$(<<M) $$(SCALPEL_OPTS) --dir $$(@D) && \
		if [! -e $$@ ]; then cp $$(@D)/main/$$(@F) $$@; fi")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call scalpel-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif

# convert annovar table to vcf and remove dcov bams
define scalpel2vcf-tumor-normal
vcf/$1_$2.scalpel_indels.vcf : scalpel/$1_$2/somatic.$(SCALPEL_MIN_COV)x.indel.annovar bam/$1.dcov.bam bam/$2.dcov.bam
	$$(INIT) tempf=`mktemp` && \
		$$(SCALPEL2VCF) -f $$(REF_FASTA) -t $1 -n $2 < $$< > $$$$tempf 2> $$(LOG) \
		&& (grep '^#' $$$$tempf; grep -v '^#' $$$$tempf | sort -V)  > $$@ 2>> $$(LOG) \
		&& rm $$$$tempf 2>> $$(LOG) \
		&& rm $$(<<) $$(<<<) 2>> $$(LOG)
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call scalpel2vcf-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

include modules/vcf_tools/vcftools.mk
include modules/bam_tools/processBam.mk
