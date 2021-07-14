include modules/Makefile.inc

LOGDIR ?= log/ntai_score.$(NOW)

ntai_score : $(foreach pair,$(SAMPLE_PAIRS),genome_stats/$(pair).ntai)

define ntai-score
genome_stats/$1_$2.ntai : facets/cncf/$1_$2.txt
	$$(call RUN,-n 1 -s 3G -m 6G,"$(RSCRIPT) modules/copy_number/ntaiscore.R --file_in $$< --file_out $$(@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ntai-score,$(tumor.$(pair)),$(normal.$(pair)))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: ntai_score
