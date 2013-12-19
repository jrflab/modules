# scalpel variant detection
# vim: set ft=make :

include ~/share/modules/Makefile.inc

LOGDIR = log/scalpel.$(NOW)
SCALPEL = $(HOME)/share/usr/scalpel-0.1.1/scalpel
SCALPEL_OPTS = --ref $(REF_FASTA)
ifeq ($(EXOME),true)
SCALPEL_OPTS += --bed $(EXOME_BED)
endif
ifdef TARGETS_FILE
SCALPEL_OPTS += --bed $(TARGETS_FILE)
endif


.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all

all : $(foreach pair,$(SAMPLE_PAIRS),scalpel/$(pair)/somatic.indel.txt)

define scalpel-tumor-normal
scalpel/$1_$2/somatic.indel.txt : bam/$1.bam bam/$2.bam
	$$(SCALPEL) --somatic --tumor $$(word 1,$$^) --normal $$(word 2,$$^) $$(SCALPEL_OPTS)
endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call scalpel-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))

