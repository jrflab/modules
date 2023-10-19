include modules/Makefile.inc

LOGDIR = log/putbam_data_mirror.$(NOW)

put_bam : $(foreach sample,$(SAMPLES),bam/$(sample).taskcomplete)
	  
PROJECT_NAME = $(basename $(PWD))

define put-bam
bam/$1.taskcomplete : bam/$1.bam
	$$(call RUN,-c -n 1 -s 2G -m 4G, "set -o pipefail && \
					  rsync -aP -e ssh bam/$1.bam $(USER)@swan.mskcc.org:/oscar/warm/reis-filho/by_user/$(USER)/$(PROJECT_NAME)/$1.bam && \
					  rsync -aP -e ssh bam/$1.bam.bai $(USER)@swan.mskcc.org:/oscar/warm/reis-filho/by_user/$(USER)/$(PROJECT_NAME)/$1.bam.bai && \
					  rsync -aP -e ssh bam/$1.bai $(USER)@swan.mskcc.org:/oscar/warm/reis-filho/by_user/$(USER)/$(PROJECT_NAME)/$1.bam.bai && \
					  echo 'finished!' > $$(@)")
					  
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call put-bam,$(sample))))

..DUMMY := $(shell mkdir -p version; \
             which scp > version/putbam_data_mirror.txt)
.SECONDARY: 
.DELETE_ON_ERROR:
.PHONY: put_bam