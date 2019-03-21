include modules/Makefile.inc

LOGDIR ?= log/titrations.$(NOW)
PHONY += titrations

pyclone : $(foreach pair,$(SAMPLE_PAIRS),titrations/$(pair).timestamp)

define make-titrations
titrations/$1_$2.timestamp : bam/$1.bam bam/$2.bam
	$$(call RUN, -s 12G -m 16G,"mkdir -p titrations && \
							    $(RSCRIPT) modules/test/copy_number/titrations.R --tumor_normal $1_$2")
							    
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call make-titrations,$(tumor.$(pair)),$(normal.$(pair)))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
