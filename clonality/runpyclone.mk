include modules/Makefile.inc

LOGDIR ?= log/run_pyclone.$(NOW)
PHONY += pyclone

run_pyclone : $(foreach set,$(SAMPLE_SETS),pyclone/$(set)/report/pyclone.tsv)

define run-pyclone
pyclone/%/trace/alpha.tsv.bz2 : pyclone/%/config.yaml
	$$(call RUN,-s 4G -m 6G -w 7200,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
							 		 PyClone run_analysis --config_file pyclone/$$*/config.yaml --seed 0")

pyclone/%/report/pyclone.tsv : pyclone/%/trace/alpha.tsv.bz2
	$$(call RUN,-s 4G -m 6G -w 7200,"make -p pyclone/$$*/report && \
									 source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
							 		 PyClone build_table --config_file pyclone/$$*/config.yaml --out_file pyclone/$$*/report/pyclone.tsv --max_cluster 10 --table_type old_style --burnin 5000")
							 
endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call run-pyclone,$(set))))
