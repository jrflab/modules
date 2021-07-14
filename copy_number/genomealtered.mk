include modules/Makefile.inc

LOGDIR ?= log/genome_altered.$(NOW)

genome_altered : $(foreach pair,$(SAMPLE_PAIRS),genome_stats/$(pair).fga)

define fraction-genome-altered
genome_stats/$1_$2.fga : facets/cncf/$1_$2.Rdata
	$$(call RUN,-n 1 -s 3G -m 6G,"$(RSCRIPT) modules/copy_number/genomealtered.R --file_in $$(<) --file_out $$(@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call fraction-genome-altered,$(tumor.$(pair)),$(normal.$(pair)))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: genome_altered
