include modules/Makefile.inc

LOGDIR ?= log/hla_polysolver.$(NOW)


hla_polysolver : $(foreach pair,$(SAMPLE_PAIRS),hla_polysolver/$(pair)/winners.hla.txt) \
		 $(foreach pair,$(SAMPLE_PAIRS),hla_polysolver/$(pair)/hla.intervals) \
		 $(foreach pair,$(SAMPLE_PAIRS),hla_polysolver/$(pair)/$(pair).mutect.unfiltered.annotated) \
		 $(foreach pair,$(SAMPLE_PAIRS),hla_polysolver/$(pair)/$(pair).strelka_indels.unfiltered.annotated) \
		 hla_polysolver/summary/hla_summary.txt \
		 hla_polysolver/summary/mutect_summary.txt \
		 hla_polysolver/summary/strelka_summary.txt
		 

define hla-polysolver
hla_polysolver/$1_$2/winners.hla.txt : bam/$1.bam bam/$2.bam
	$$(call RUN,-c -n 8 -s 2G -m 4G -v $(POLYSOLVER_ENV) -w 72:00:00, "set -o pipefail && \
									   export CONDA_PREFIX=$$(POLYSOLVER_ENV) && \
									   export PERL5LIB=$$(POLYSOLVER_ENV)/lib/perl5/5.22.0 && \
									   shell_call_hla_type \
									   $$(<<) \
									   Unknown \
									   1 \
									   hg19 \
									   STDFQ \
									   0 \
									   hla_polysolver/$1_$2")
								 	  
hla_polysolver/$1_$2/hla.intervals : bam/$1.bam bam/$2.bam hla_polysolver/$1_$2/winners.hla.txt
	$$(call RUN,-c -n 8 -s 2G -m 4G -v $(POLYSOLVER_ENV) -w 72:00:00, "set -o pipefail && \
									   export CONDA_PREFIX=$$(POLYSOLVER_ENV) && \
									   export PERL5LIB=$$(POLYSOLVER_ENV)/lib/perl5/5.22.0 && \
									   shell_call_hla_mutations_from_type \
									   $$(<<) \
									   $$(<) \
									   $$(<<<) \
									   hg19 \
									   STDFQ \
									   hla_polysolver/$1_$2")
								 	 
hla_polysolver/$1_$2/$1_$2.mutect.unfiltered.annotated : hla_polysolver/$1_$2/hla.intervals
	$$(call RUN,-c -n 8 -s 2G -m 4G -v $(POLYSOLVER_ENV) -w 72:00:00, "set -o pipefail && \
									   export CONDA_PREFIX=$$(POLYSOLVER_ENV) && \
									   export PERL5LIB=$$(POLYSOLVER_ENV)/lib/perl5/5.22.0 && \
									   shell_annotate_hla_mutations \
									   $1_$2 \
									   hla_polysolver/$1_$2")

hla_polysolver/$1_$2/$1_$2.strelka_indels.unfiltered.annotated : hla_polysolver/$1_$2/$1_$2.mutect.unfiltered.annotated
	

endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call hla-polysolver,$(tumor.$(pair)),$(normal.$(pair)))))
		
hla_polysolver/summary/hla_summary.txt : $(foreach pair,$(SAMPLE_PAIRS),hla_polysolver/$(pair)/$(pair).mutect.unfiltered.annotated) $(foreach pair,$(SAMPLE_PAIRS),hla_polysolver/$(pair)/$(pair).strelka_indels.unfiltered.annotated)
	$(call RUN,-c -s 12G -m 24G,"set -o pipefail && \
				     $(RSCRIPT) modules/variant_callers/somatic/hla_summary.R --option 1 --sample_names '$(SAMPLE_PAIRS)'")

hla_polysolver/summary/mutect_summary.txt : $(foreach pair,$(SAMPLE_PAIRS),hla_polysolver/$(pair)/$(pair).mutect.unfiltered.annotated) $(foreach pair,$(SAMPLE_PAIRS),hla_polysolver/$(pair)/$(pair).strelka_indels.unfiltered.annotated)
	$(call RUN,-c -s 12G -m 24G,"set -o pipefail && \
				     $(RSCRIPT) modules/variant_callers/somatic/hla_summary.R --option 2 --sample_names '$(SAMPLE_PAIRS)'")

hla_polysolver/summary/strelka_summary.txt : $(foreach pair,$(SAMPLE_PAIRS),hla_polysolver/$(pair)/$(pair).mutect.unfiltered.annotated) $(foreach pair,$(SAMPLE_PAIRS),hla_polysolver/$(pair)/$(pair).strelka_indels.unfiltered.annotated)
	$(call RUN,-c -s 12G -m 24G,"set -o pipefail && \
				     $(RSCRIPT) modules/variant_callers/somatic/hla_summary.R --option 3 --sample_names '$(SAMPLE_PAIRS)'")


..DUMMY := $(shell mkdir -p version; \
	     $(POLYSOLVER_ENV)/bin/shell_call_hla_type --help &> version/hla_polysolver.txt; \
	     $(POLYSOLVER_ENV)/bin/shell_call_hla_mutations_from_type --help &>> version/hla_polysolver.txt; \
	     $(POLYSOLVER_ENV)/bin/shell_annotate_hla_mutations --help &>> version/hla_polysolver.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: hla_polysolver
