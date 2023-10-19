include modules/Makefile.inc

LOGDIR = log/getbam_data_mirror.$(NOW)

get_bam : $(foreach sample,$(SAMPLES),bam/$(sample).bam) \
	  $(foreach sample,$(SAMPLES),bam/$(sample).bam.bai) \
	  $(foreach sample,$(SAMPLES),bam/$(sample).bai)
	  
PROJECT_NAME = $(shell basename $(PWD))

define get-bam
bam/$1.bam :
	$$(call RUN,-c -n 1 -s 2G -m 4G, "set -o pipefail && \
					  rsync -aP -e ssh $(USER)@swan.mskcc.org:/oscar/warm/reis-filho/by_user/$(USER)/$(PROJECT_NAME)/$1.bam \
					  bam/")
					  
bam/$1.bam.bai :
	$$(call RUN,-c -n 1 -s 2G -m 4G, "set -o pipefail && \
					  rsync -aP -e ssh $(USER)@swan.mskcc.org:/oscar/warm/reis-filho/by_user/$(USER)/$(PROJECT_NAME)/$1.bam.bai \
					  bam/")
					  
bam/$1.bai :
	$$(call RUN,-c -n 1 -s 2G -m 4G, "set -o pipefail && \
					  rsync -aP -e ssh $(USER)@swan.mskcc.org:/oscar/warm/reis-filho/by_user/$(USER)/$(PROJECT_NAME)/$1.bai \
					  bam/")


endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call get-bam,$(sample))))

..DUMMY := $(shell mkdir -p version; \
             which scp > version/getbam_data_mirror.txt)
.SECONDARY: 
.DELETE_ON_ERROR:
.PHONY: get_bam