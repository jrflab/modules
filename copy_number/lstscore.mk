include modules/Makefile.inc

LOGDIR ?= log/lst_score.$(NOW)
PHONY += genome_stats

lst_score : $(foreach pair,$(SAMPLE_PAIRS),genome_stats/$(pair).lst)

define lst-score
genome_stats/$1_$2.lst : facets/cncf/$1_$2.txt
	$$(call RUN,-n 1 -s 3G -m 6G,"$(RSCRIPT) modules/copy_number/lstscore.R --file_in $$< --file_out genome_stats/$1_$2.lst")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call lst-score,$(tumor.$(pair)),$(normal.$(pair)))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
