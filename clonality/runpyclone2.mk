include modules/Makefile.inc

LOGDIR ?= log/run_pyclone2.$(NOW)
PHONY += pyclone

run_pyclone2 : $(foreach pair,$(SAMPLE_PAIRS),pyclone/$(pair)/pyclone.tsv)

define run-pyclone
pyclone/$1_$2/trace/alpha.tsv.bz2 : pyclone/$1_$2/config.yaml
	$$(call RUN,-s 4G -m 6G,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
							 PyClone run_analysis --config_file pyclone/$1_$2/config.yaml --seed 0")

pyclone/$1_$2/pyclone.tsv : pyclone/$1_$2/trace/alpha.tsv.bz2
	$$(call RUN,-s 4G -m 6G,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
							 PyClone build_table --config_file pyclone/$1_$2/config.yaml --out_file pyclone/$1_$2/pyclone.tsv --max_cluster 5 --table_type old_style --burnin 1000")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call run-pyclone,$(tumor.$(pair)),$(normal.$(pair)))))
