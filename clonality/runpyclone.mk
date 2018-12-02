include modules/Makefile.inc

LOGDIR ?= log/run_pyclone.$(NOW)
PHONY += pyclone

run_pyclone : $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/pyclone.tsv)

define run-pyclone
pyclone/%/trace/alpha.tsv.bz2 : pyclone/%/config.yaml
	$$(call RUN,-s 4G -m 6G -w 7200,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
							 PyClone run_analysis --config_file pyclone/$$*/config.yaml --seed 0")

pyclone/%/pyclone.tsv : pyclone/%/trace/alpha.tsv.bz2
	$$(call RUN,-s 4G -m 6G -w 7200,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
							 PyClone build_table --config_file pyclone/$$*/config.yaml --out_file pyclone/$$*/pyclone.tsv --max_cluster 5 --table_type old_style --burnin 1000")
							 
endef
$(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call run-pyclone,$(sample))))
