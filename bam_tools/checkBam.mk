include modules/Makefile.inc

LOGDIR ?= log/check_bam.$(NOW)
PHONY += check_bam

check_bam : $(foreach sample,$(SAMPLES),check_bam/$(sample).txt)

CHECK ?= $(wildcard $(foreach set,$(SAMPLES),check_bam/$(set).txt))

define check-bam
check_bam/%.txt : bam/%.bam
	$$(call RUN,-c -n 1 -s 2G -m 4G," if [ -f $$(<) ]; then touch $$(<<); fi")
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call check-bam,$(sample))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)