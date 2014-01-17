# scalpel variant detection
# vim: set ft=make :

include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

LOGDIR = log/scalpel.$(NOW)

SCALPEL = $(HOME)/share/usr/scalpel-0.1.1/scalpel
SCALPEL_OPTS = --ref $(REF_FASTA)
ifeq ($(EXOME),true)
BED_DIR = $(HOME)/share/reference/splitExonBed/
BED_FILES = $(shell ls $(BED_DIR))
endif
ifdef TARGETS_FILE
SCALPEL_OPTS += --bed $(TARGETS_FILE)
endif

SCALPEL2VCF = $(PERL) $(HOME)/share/scripts/scalpelToVcf.pl

FILTER_SUFFIX := dbsnp.nsfp.chasm.fathmm
EFF_TYPES = silent missense nonsilent_cds nonsilent
TABLE_SUFFIXES = $(foreach eff,$(EFF_TYPES),$(FILTER_SUFFIX).tab.$(eff).pass.novel)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all vcfs tables

all : vcfs tables

vcfs : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).scalpel.$(FILTER_SUFFIX).vcf)
tables : $(foreach pair,$(SAMPLE_PAIRS),$(foreach suff,$(TABLE_SUFFIXES),tables/$(pair).scalpel.$(suff).txt))

ifdef BED_FILES
define scalpel-bed-tumor-normal
scalpel/$2_$3/$1/somatic.5x.indel.txt : bam/$2.bam bam/$3.bam
	$$(call LSCRIPT_NAMED_PARALLEL_MEM,$2_$3_$1_scalpel,2,4G,7G,"$$(SCALPEL) --somatic --numprocs 2 --tumor $$(word 1,$$^) --normal $$(word 2,$$^) $$(SCALPEL_OPTS) --bed $$(BED_DIR)/$1 --dir $$(@D)")
endef
$(foreach bed,$(BED_FILES),\
	$(foreach i,$(SETS_SEQ),\
		$(foreach tumor,$(call get_tumors,$(set.$i)), \
			$(eval $(call scalpel-bed-tumor-normal,$(bed),$(tumor),$(call get_normal,$(set.$i)))))))

define merge-scalpel-tumor-normal
scalpel/tables/$1_$2.scalpel.txt : $$(foreach bed,$$(BED_FILES),scalpel/$1_$2/$$(bed)/somatic.5x.indel.txt)
	$$(INIT) head -1 $$< > $$@ && foreach x in $$^; do sed '1d' $$$$x >> $$@; done
endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call merge-scalpel-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))

define scalpel2vcf-tumor-normal
vcf/$1_$2.scalpel.vcf : scalpel/tables/$1_$2.scalpel.txt
	$$(INIT) $$(SCALPEL2VCF) -f $$(REF_FASTA) -t $1 -n $2 < $$< > $$@ 2> $$(LOG)
endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call scalpel2vcf-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))
else
define scalpel-tumor-normal
scalpel/$1_$2/somatic.5x.indel.txt : bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_NAMED_PARALLEL_MEM,$1_$2_scalpel,8,1G,2G,"$$(SCALPEL) --somatic --numprocs 8 --tumor $$(word 1,$$^) --normal $$(word 2,$$^) $$(SCALPEL_OPTS) --dir $$(@D)")
endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call scalpel-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))
define scalpel2vcf-tumor-normal
vcf/$1_$2.scalpel.vcf : scalpel/$1_$2/somatic.5x.indel.txt
	$$(INIT) $$(SCALPEL2VCF) -f $$(REF_FASTA) -t $1 -n $2 < $$< > $$@ 2> $$(LOG)
endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call scalpel2vcf-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))
endif


include ~/share/modules/vcftools.mk


