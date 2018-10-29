include modules/Makefile.inc

LOGDIR ?= log/hotspot_summary.$(NOW)
PHONY += hotspot

HOTSPOT ?= $(wildcard $(foreach set,$(SAMPLE_PAIRS),hotspot/$(set).txt))

hotspot_summary : hotspot/hotspot_summary.txt

hotspot/hotspot_summary.txt : $(wildcard hotspot/$(SAMPLE_PAIRS).txt)
	$(call RUN,-n 1 -s 4G -m 4G,"$(RSCRIPT) modules/summary/hotspotsummary.R --in_file '$(HOTSPOT)' --out_file hotspot/hotspot_summary.txt")
		

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
