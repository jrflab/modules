include modules/Makefile.inc

LOGDIR = log/get_bam.$(NOW)

get_bam : $(foreach sample,$(SAMPLES),bam/$(sample).bam)

define get-bam
bam/$1.bam :
	$$(call RUN,-c -n 1 -s 2G -m 4G, "set -o pipefail && \
					  scp $(USER)@selene.mskcc.org:/res/dmpcollab/dmpshare/share/irb12_245/$(echo $1 | cut -c 1-1)/`echo $1 | cut -c 2-2`/$1.bam \
					  bam/")

endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call get-bam,$(sample))))

..DUMMY := $(shell mkdir -p version; \
             which scp > version/get_bam.txt)
.SECONDARY: 
.DELETE_ON_ERROR:
.PHONY: get_bam