include modules/Makefile.inc

LOGDIR ?= log/hla_polysolver.$(NOW)
PHONY += hla_polysolver

hla_polysolver : $(foreach pair,$(SAMPLE_PAIRS),hla_polysolver/$(pair)/$(pair).taskcomplete)

define hla-polysolver
hla_polysolver/$1_$2/winners.hla.txt : bam/$1.bam bam/$2.bam
	$$(call RUN,-n 8 -s 12G -m 24G, "source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate \
									 /home/${USER}/share/usr/anaconda-envs/hla-polysolver && \
									 export CONDA_PREFIX=/home/${USER}/share/usr/anaconda-envs/hla-polysolver && \
									 export PERL5LIB=/home/${USER}/share/usr/anaconda-envs/hla-polysolver/lib/perl5/5.22.0 && \
								 	 if [ ! -d hla_polysolver/$1_$2 ]; then mkdir hla_polysolver/$1_$2; fi && \
								 	 shell_call_hla_type bam/$2.bam Unknown 1 hg19 STDFQ 0 hla_polysolver/$1_$2")
								 	  
hla_polysolver/$1_$2/hla.intervals : hla_polysolver/$1_$2/winners.hla.txt
	$$(call RUN,-n 8 -s 12G -m 24G, "source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate \
									 /home/${USER}/share/usr/anaconda-envs/hla-polysolver && \
									 export CONDA_PREFIX=/home/${USER}/share/usr/anaconda-envs/hla-polysolver && \
									 export PERL5LIB=/home/${USER}/share/usr/anaconda-envs/hla-polysolver/lib/perl5/5.22.0 && \
								 	 shell_call_hla_mutations_from_type bam/$2.bam bam/$1.bam hla_polysolver/$1_$2/winners.hla.txt hg19 STDFQ hla_polysolver/$1_$2")
								 	 
hla_polysolver/$1_$2/$1_$2.mutect.unfiltered.annotated hla_polysolver/$1_$2/$1_$2.strelka_indels.unfiltered.annotated : hla_polysolver/$1_$2/hla.intervals
	$$(call RUN,-n 8 -s 12G -m 24G, "source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate \
									 /home/${USER}/share/usr/anaconda-envs/hla-polysolver && \
									 export CONDA_PREFIX=/home/${USER}/share/usr/anaconda-envs/hla-polysolver && \
									 export PERL5LIB=/home/${USER}/share/usr/anaconda-envs/hla-polysolver/lib/perl5/5.22.0 && \
								 	 shell_annotate_hla_mutations $1_$2 hla_polysolver/$1_$2")

hla_polysolver/$1_$2/$1_$2.taskcomplete : hla_polysolver/$1_$2/$1_$2.mutect.unfiltered.annotated hla_polysolver/$1_$2/$1_$2.strelka_indels.unfiltered.annotated
	$$(call RUN,-s 1G -m 1G,"touch hla_polysolver/$1_$2/$1_$2.taskcomplete")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call hla-polysolver,$(tumor.$(pair)),$(normal.$(pair)))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
