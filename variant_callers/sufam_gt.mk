include modules/Makefile.inc

LOGDIR ?= log/sufam_gt.$(NOW)

SUFAM_ENV = $(HOME)/share/usr/anaconda-envs/sufam-dev
SUFAM_OPTS = --mpileup-parameters='-A -q 15 -Q 15 -d 15000'

sufam_gt : $(foreach sample,$(TUMOR_SAMPLES),sufam/$(sample).vcf)

define tsv-2-vcf
sufam/$1.vcf : summary/tsv/all.tsv
	$$(call RUN,-c -n 1 -s 4G -m 8G,"set -o pipefail && \
					 $(RSCRIPT) $(SCRIPTS_DIR)/sufam_gt.R \
					 --option 1 \
					 --sample $1 \
					 --sample_set '$(set.$1)' \
					 --normal_sample '$(normal.$1)' \
					 --input_file $$(<) \
					 --output_file $$(@)")

endef
$(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call tsv-2-vcf,$(sample))))
		
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY:
