include modules/Makefile.inc

LOGDIR ?= log/mimsi.$(NOW)

mimsi: $(foreach pair,$(SAMPLE_PAIRS),mimsi/$(pair)/$(pair).txt)

MICROSATELLITES_LIST = $(HOME)/share/lib/resource_files/mimsi/microsatellites_impact_only.list
MODEL = $(HOME)/share/lib/resource_files/mimsi/mi_msi_v0_4_0_200x.model

define mimsi-tumor-normal
mimsi/$1_$2/$1_$2.txt : bam/$1.bam bam/$2.bam
	$$(call RUN,-c -n 8 -s 1G -m 2G -v $(MIMSI_ENV),"set -o pipefail && \
							 mkdir -p mimsi/$1_$2/ && \
							 analyze \
							 --tumor-bam $$(<) \
							 --normal-bam $$(<<) \
							 --case-id $1 \
							 --norm-case-id $2 \
							 --microsatellites-list $$(MICROSATELLITES_LIST) \
							 --save-location mimsi/$1_$2/ \
							 --model $$(MODEL) \
							 --save && \
							 mv mimsi/$1_$2/BATCH_results.txt $$(@)")

endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call mimsi-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: mimsi
