include modules/Makefile.inc

LOGDIR ?= log/genome_summary.$(NOW)


genome_summary : $(foreach pair,$(SAMPLE_PAIRS),genome_stats/$(pair).fga) \
		 genome_stats/genome_altered.tsv \
		 $(foreach pair,$(SAMPLE_PAIRS),genome_stats/$(pair).lst) \
		 genome_stats/lst_score.tsv \
		 $(foreach pair,$(SAMPLE_PAIRS),genome_stats/$(pair).ntai) \
		 genome_stats/ntai_score.tsv \
		 $(foreach pair,$(SAMPLE_PAIRS),genome_stats/$(pair).mrs) \
		 genome_stats/myriad_score.tsv
#		 summary/tsv/genome_summary.tsv \
#		 summary/genome_summary.xlsx
		 
define fraction-genome-altered
genome_stats/$1_$2.fga : facets/cncf/$1_$2.Rdata
	$$(call RUN,-n 1 -s 3G -m 6G,"$(RSCRIPT) modules/copy_number/genomealtered.R --file_in $$(<) --file_out $$(@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call fraction-genome-altered,$(tumor.$(pair)),$(normal.$(pair)))))
		
define lst-score
genome_stats/$1_$2.lst : facets/cncf/$1_$2.txt
	$$(call RUN,-n 1 -s 3G -m 6G,"$(RSCRIPT) modules/copy_number/lstscore.R --file_in $$< --file_out $$(@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call lst-score,$(tumor.$(pair)),$(normal.$(pair)))))
		
define ntai-score
genome_stats/$1_$2.ntai : facets/cncf/$1_$2.txt
	$$(call RUN,-n 1 -s 3G -m 6G,"$(RSCRIPT) modules/copy_number/ntaiscore.R --file_in $$< --file_out $$(@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ntai-score,$(tumor.$(pair)),$(normal.$(pair)))))
		
define myriad-score
genome_stats/$1_$2.mrs : facets/cncf/$1_$2.txt
	$$(call RUN,-n 1 -s 3G -m 6G,"$(RSCRIPT) modules/copy_number/myriadhrdscore.R --file_in $$< --file_out $$(@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call myriad-score,$(tumor.$(pair)),$(normal.$(pair)))))

#genome_stats/genome_altered.tsv : $(GENOME_ALTERED)
#	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
#				     mkdir -p genome_stats && \
#				     cat $$(GENOME_ALTERED) > $$(@)")
#							 
#genome_stats/lst_score.tsv : $(LST_SCORE)
#	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
#				     mkdir -p genome_stats && \
#				     cat $(LST_SCORE) > $$(@)")
#				     
#genome_stats/ntai_score.tsv : $(NTAI_SCORE)
#	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
#				     mkdir -p genome_stats && \
#				     cat $(NTAI_SCORE) > $$(@)")
#
#genome_stats/myriad_score.tsv : $(MYRIAD_SCORE)
#	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
#				     mkdir -p genome_stats && \
#				     cat $(MYRIAD_SCORE) > $$(@)")
#
#summary/tsv/genome_summary.tsv : genome_stats/genome_altered.tsv genome_stats/lst_score.tsv genome_stats/ntai_score.tsv genome_stats/myriad_score.tsv
#	$(call RUN,-n 1 -s 6G -m 8G,"set -o pipefail && \
#				     mkdir -p genome_stats && \
#				     $(RSCRIPT) modules/summary/genomesummary.R")
#
#summary/genome_summary.xlsx : summary/tsv/genome_summary.tsv
#	$(call RUN,-n 1 -s 4G -m 4G,"python modules/summary/genome_summary_excel.py")

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: genome_sumary
