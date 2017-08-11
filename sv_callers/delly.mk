# run delly
LOGDIR = log/delly.$(NOW)

include modules/Makefile.inc

DELLY_ENV = $(HOME)/share/usr/anaconda-envs/delly-0.7.6

SUAVE_BAM_TO_H5 = python $(HOME)/share/usr/delly/vis/suave/suave_bam_to_h5.py

DELLY_TYPES = DEL DUP INV TRA INS

.SECONDARY:
.DELETE_ON_ERROR:

.PHONY: delly

delly: $(foreach pair,$(SAMPLE_PAIRS),$(foreach type,$(DELLY_TYPES),delly/bcf/$(pair).delly_$(type).bcf)) \
	$(foreach sample,$(SAMPLES),delly/h5/$(sample).h5)

define delly-pair-type
delly/bcf/$1_$2.delly_$3.bcf : bam/$1.bam bam/$2.bam bam/$1.bam.bai bam/$2.bam.bai
	$$(call RUN,-v $$(DELLY_ENV) -w 256:00:00,8G,12G,"delly call -t $3 -o $$@ -g $$(REF_FASTA) $$< $$(<<)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach type,$(DELLY_TYPES),\
	$(eval $(call delly-pair-type,$(tumor.$(pair)),$(normal.$(pair)),$(type)))))

delly/h5/%.h5 : bam/%.bam
	$(call RUN,-v $(DELLY_ENV) -s 8G -m 10G,"$(SUAVE_BAM_TO_H5) -s $* -c gzip -o $@ $<")


