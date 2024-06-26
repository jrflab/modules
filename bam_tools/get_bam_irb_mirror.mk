include modules/Makefile.inc

LOGDIR = log/getbam_irb_mirror.$(NOW)

get_bam : $(foreach sample,$(SAMPLES),bam/$(sample).bam) \
	  $(foreach sample,$(SAMPLES),bam/$(sample).bam.bai) \
	  $(foreach sample,$(SAMPLES),bam/$(sample).bai)

define get-bam
bam/$1.bam :
	$$(call RUN,-c -n 1 -s 2G -m 4G, "set -o pipefail && \
					  scp $(USER)@juno-xfer01.mskcc.org:/juno/dmp/share/irb12_245/`echo $1 | cut -c 1-1`/`echo $1 | cut -c 2-2`/$1.bam \
					  bam/")
					  
bam/$1.bam.bai : bam/$1.bam
	$$(call RUN,-c -n 1 -s 2G -m 4G, "set -o pipefail && \
					  $(SAMTOOLS) index $$(<)")
					  
bam/$1.bai : bam/$1.bam bam/$1.bam.bai
	$$(call RUN,-c -n 1 -s 2G -m 4G, "set -o pipefail && \
					  cp $$(<<) $$(@)")


endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call get-bam,$(sample))))

..DUMMY := $(shell mkdir -p version; \
             which scp > version/getbam_irb_mirror.txt)
.SECONDARY: 
.DELETE_ON_ERROR:
.PHONY: get_bam