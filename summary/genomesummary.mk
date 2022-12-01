include modules/Makefile.inc

LOGDIR ?= log/genome_summary.$(NOW)


genome_summary : $(foreach pair,$(SAMPLE_PAIRS),genome_summary/genome_altered/$(pair).txt) \
		 $(foreach pair,$(SAMPLE_PAIRS),genome_summary/lst/$(pair).txt) \
		 $(foreach pair,$(SAMPLE_PAIRS),genome_summary/ntai/$(pair).txt) \
		 $(foreach pair,$(SAMPLE_PAIRS),genome_summary/myriad_score/$(pair).txt)

#		 genome_summary/genome_altered/summary.txt
#		 genome_stats/myriad_score.tsv \
#		 summary/tsv/genome_summary.tsv \
#		 summary/genome_summary.xlsx
		 
define fraction-genome-altered
genome_summary/genome_altered/$1_$2.txt : facets/cncf/$1_$2.Rdata
	$$(call RUN,-n 1 -s 3G -m 6G,"set -o pipefail && \
				      $(RSCRIPT) modules/summary/genomesummary.R \
				      --option 1 \
				      --sample_name $1_$2 \
				      --file_in $$(<) \
				      --file_out $$(@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call fraction-genome-altered,$(tumor.$(pair)),$(normal.$(pair)))))
		
define lst-score
genome_summary/lst/$1_$2.txt : facets/cncf/$1_$2.txt
	$$(call RUN,-n 1 -s 3G -m 6G,"set -o pipefail && \
				      $(RSCRIPT) modules/summary/genomesummary.R \
				      --option 2 \
				      --sample_name $1_$2 \
				      --file_in $$(<) \
				      --file_out $$(@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call lst-score,$(tumor.$(pair)),$(normal.$(pair)))))
		
define ntai-score
genome_summary/ntai/$1_$2.txt : facets/cncf/$1_$2.txt
	$$(call RUN,-n 1 -s 3G -m 6G,"set -o pipefail && \
				      $(RSCRIPT) modules/summary/genomesummary.R \
				      --option 3 \
				      --sample_name $1_$2 \
				      --file_in $$(<) \
				      --file_out $$(@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ntai-score,$(tumor.$(pair)),$(normal.$(pair)))))
		
define myriad-score
genome_summary/myriad_score/$1_$2.txt : facets/cncf/$1_$2.txt
	$$(call RUN,-n 1 -s 3G -m 6G,"set -o pipefail && \
				      $(RSCRIPT) modules/summary/genomesummary.R \
				      --option 4 \
				      --sample_name $1_$2 \
				      --file_in $$(<) \
				      --file_out $$(@)")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call myriad-score,$(tumor.$(pair)),$(normal.$(pair)))))

genome_summary/genome_altered/summary.txt : $(foreach pair,$(SAMPLE_PAIRS),genome_summary/genome_altered/$(pair).txt)
	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
				     ")
							 
genome_stats/lst_score.tsv : $(LST_SCORE)
	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
				     mkdir -p genome_stats && \
				     cat $(LST_SCORE) > $(@)")
				     
genome_stats/ntai_score.tsv : $(NTAI_SCORE)
	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
				     mkdir -p genome_stats && \
				     cat $(NTAI_SCORE) > $(@)")

genome_stats/myriad_score.tsv : $(MYRIAD_SCORE)
	$(call RUN,-n 1 -s 4G -m 4G,"set -o pipefail && \
				     mkdir -p genome_stats && \
				     cat $(MYRIAD_SCORE) > $(@)")

summary/tsv/genome_summary.tsv : genome_stats/genome_altered.tsv genome_stats/lst_score.tsv genome_stats/ntai_score.tsv genome_stats/myriad_score.tsv
	$(call RUN,-n 1 -s 6G -m 8G,"set -o pipefail && \
				     mkdir -p genome_stats && \
				     $(RSCRIPT) modules/summary/genomesummary.R --option 5")

summary/genome_summary.xlsx : summary/tsv/genome_summary.tsv
	$(call RUN,-n 1 -s 4G -m 4G,"python modules/summary/genome_summary_excel.py")

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: genome_summary
