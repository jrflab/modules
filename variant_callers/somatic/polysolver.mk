include modules/Makefile.inc

LOGDIR ?= log/hla_polysolver.$(NOW)
PHONY += hla_polysolver

hla_polysolver : $(foreach pair,$(SAMPLE_PAIRS),hla_polysolver/$(pair)/$(pair).taskcomplete)

define hla-polysolver
hla_polysolver/$1_$2/winners.hla.txt : bam/$1.bam bam/$2.bam
	$$(call RUN,-c -n 8 -s 12G -m 24G -v /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/, \
									  "source activate /home/${USER}/share/usr/anaconda-envs/hla-polysolver && \
									   export CONDA_PREFIX=/home/${USER}/share/usr/anaconda-envs/hla-polysolver && \
									   export PERL5LIB=/home/${USER}/share/usr/anaconda-envs/hla-polysolver/lib/perl5/5.22.0/ && \
								 	   if [ -d hla_polysolver/$1_$2 ]; then mkdir hla_polysolver/$1_$2; fi && \
								 	   shell_call_hla_type bam/$2.bam Unknown 1 hg19 STDFQ 0 hla-polysolver/$1_$2")

hla_polysolver/$1_$2/$1_$2.taskcomplete : hla_polysolver/$1_$2/winners.hla.txt
	$$(call RUN,-c -s 1G -m 1G,"echo $$< $$(<<) > hla_polysolver/$1_$2/$1_$2.taskcomplete")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call hla-polysolver,$(tumor.$(pair)),$(normal.$(pair)))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
