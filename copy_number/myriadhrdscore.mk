include modules/Makefile.inc

LOGDIR ?= log/myriad_score.$(NOW)
PHONY += genome_stats

myriad_score : $(foreach pair,$(SAMPLE_PAIRS),genome_stats/$(pair).mrs)

define myriad-score
genome_stats/$1_$2.mrs : facets/cncf/$1_$2.txt
	$$(call RUN,-n 1 -s 3G -m 6G,"$(RSCRIPT) modules/copy_number/myriadhrdscore.R --file_in $$< --file_out genome_stats/$1_$2.mrs")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call myriad-score,$(tumor.$(pair)),$(normal.$(pair)))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
