include modules/Makefile.inc

LOGDIR ?= log/setup_pyclone.$(NOW)
PHONY += pyclone

pyclone : $(foreach pair,$(SAMPLE_PAIRS),pyclone/$(pair)/trace/alpha.tsv.bz2)

define make-pyclone
pyclone/$1_$2/config.yaml : summary/tsv/mutation_summary.tsv
	$$(call RUN, -s 16G -m 24G,"mkdir -p pyclone/$1_$2 && \
							    $(RSCRIPT) modules/test/clonality/tsvtopyclone.R --sample_name $1_$2 --sample_fraction '$${pyclone_fraction.$1}'")
							    
pyclone/$1_$2/trace/alpha.tsv.bz2 : pyclone/$1_$2/config.yaml
	$$(call RUN,-s 8G -m 16G -w 7200,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
							 		  PyClone run_analysis --config_file pyclone/$1_$2/config.yaml --seed 0")

endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call make-pyclone,$(tumor.$(pair)),$(normal.$(pair)))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
