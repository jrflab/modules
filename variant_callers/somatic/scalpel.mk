# scalpel variant detection
# vim: set ft=make :

include ~/share/modules/Makefile.inc
include ~/share/modules/variant_callers/gatk.inc

LOGDIR = log/scalpel.$(NOW)

SCALPEL_DIR = $(HOME)/share/usr/scalpel-0.3.2
SCALPEL = export LD_LIBRARY_PATH=$(LD_LIBRARY_PATH):$(SCALPEL_DIR)/bamtools-2.3.0/lib/; $(PERL) $(SCALPEL_DIR)/scalpel
SCALPEL_OPTS = --ref $(REF_FASTA)
ifeq ($(EXOME),true)
BED_DIR = $(HOME)/share/reference/splitExonBed/
BED_FILES = $(shell ls $(BED_DIR))
endif
ifdef TARGETS_FILE
SCALPEL_OPTS += --bed $(TARGETS_FILE)
endif

SCALPEL2VCF = $(PERL) $(HOME)/share/scripts/scalpelToVcf.pl

FILTER_SUFFIX := dbsnp.eff
EFF_TYPES = silent missense nonsilent_cds nonsilent
TABLE_SUFFIXES = $(foreach eff,$(EFF_TYPES),$(FILTER_SUFFIX).tab.$(eff).pass.novel)

..DUMMY := $(shell mkdir -p version; echo "$(SCALPEL) $(SCALPEL_OPTS) > version/scalpel.txt")

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all vcfs tables

all : vcfs tables

vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).scalpel.$(FILTER_SUFFIX).vcf)
tables : $(foreach pair,$(SAMPLE_PAIRS),$(foreach suff,$(TABLE_SUFFIXES),tables/$(pair).scalpel.$(suff).txt)) $(foreach suff,$(TABLE_SUFFIXES),alltables/allTN.scalpel.$(suff).txt)

ifdef BED_FILES
define scalpel-bed-tumor-normal
scalpel/$2_$3/$1/somatic.5x.indel.vcf : bam/$2.dcov.bam.md5 bam/$3.dcov.bam.md5
	$$(call LSCRIPT_NAMED_PARALLEL_MEM,$2_$3_$1_scalpel,2,4G,7G,"$$(SCALPEL) --somatic --numprocs 2 --tumor $$(word 1,$$^) --normal $$(word 2,$$^) $$(SCALPEL_OPTS) --bed $$(BED_DIR)/$1 --dir $$(@D)")
endef
$(foreach bed,$(BED_FILES),\
	$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call scalpel-bed-tumor-normal,$(bed),$(tumor.$(pair)),$(normal.$(pair))))))

define merge-scalpel-tumor-normal
scalpel/vcf/$1_$2.scalpel.vcf : $$(foreach bed,$$(BED_FILES),scalpel/$1_$2/$$(bed)/somatic.5x.indel.vcf)
	$$(INIT) grep '^#' $$< > $$@ && for x in $$^; grep -v '^#' $$$$x >> $$@; done
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call merge-scalpel-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

#define scalpel2vcf-tumor-normal
#vcf/$1_$2.scalpel.vcf : scalpel/tables/$1_$2.scalpel.txt
#$$(INIT) $$(SCALPEL2VCF) -f $$(REF_FASTA) -t $1 -n $2 < $$< > $$@ 2> $$(LOG)
#endef
#$(foreach i,$(SETS_SEQ),\
#$(foreach tumor,$(call get_tumors,$(set.$i)), \
#$(eval $(call scalpel2vcf-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))
else

define scalpel-tumor-normal
scalpel/$1_$2/somatic.5x.indel.vcf : bam/$1.dcov.bam.md5 bam/$2.dcov.bam.md5
	$$(call LSCRIPT_NAMED_PARALLEL_MEM,$1_$2_scalpel,8,1G,2.5G,"$$(SCALPEL) --somatic --numprocs 8 --tumor $$(word 1,$$^) --normal $$(word 2,$$^) $$(SCALPEL_OPTS) --dir $$(@D)")

#vcf/$1_$2.scalpel.vcf : scalpel/$1_$2/somatic.5x.indel.txt
#$$(INIT) $$(SCALPEL2VCF) -f $$(REF_FASTA) -t $1 -n $2 < $$< > $$@ 2> $$(LOG)
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call scalpel-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif

vcf/%.scalpel.vcf : scalpel/%/somatic.5x.indel.vcf
	$(INIT) cp $< $@

include ~/share/modules/vcf_tools/vcftools.mk
include ~/share/modules/bam_tools/processBam.mk
