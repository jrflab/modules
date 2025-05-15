include modules/Makefile.inc

LOGDIR ?= log/mimsi.$(NOW)

mimsi: $(foreach pair,$(SAMPLE_PAIRS),mimsi/$(pair)/$(pair).txt) \
       mimsi/summary.txt

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
	
mimsi/summary.txt : $(foreach pair,$(SAMPLE_PAIRS),mimsi/$(pair)/$(pair).txt)
	$(call RUN, -c -n 1 -s 8G -m 12G -v $(INNOVATION_ENV),"set -o pipefail && \
							       $(RSCRIPT) $(SCRIPTS_DIR)/mimsi.R --option 1 --sample_names '$(SAMPLE_PAIRS)'")


.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: mimsi
