include modules/Makefile.inc

LOGDIR ?= log/run_pyclone.$(NOW)
PHONY += pyclone

run_pyclone : $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/pyclone.tsv) $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/pyclone_density.pdf) $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/pyclone_matrix.pdf) $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/pyclone_vaf_coordinates.pdf) $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/pyclone_scatter.pdf)

define run-pyclone
pyclone/%/trace/alpha.tsv.bz2 : pyclone/%/config.yaml
	$$(call RUN,-s 4G -m 6G,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
							 PyClone run_analysis --config_file pyclone/$$*/config.yaml --seed 0")

pyclone/%/pyclone.tsv : pyclone/%/trace/alpha.tsv.bz2
	$$(call RUN,-s 4G -m 6G,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
							 PyClone build_table --config_file pyclone/$$*/config.yaml --out_file pyclone/$$*/pyclone.tsv --max_cluster 5 --table_type old_style --burnin 1000")
							 
pyclone/%/pyclone_density.pdf : pyclone/%/trace/alpha.tsv.bz2
	$$(call RUN,-s 4G -m 6G,"export DISPLAY='localhost:10.0' && \
							 source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
							 PyClone plot_loci --config_file pyclone/$$*/config.yaml --plot_file pyclone/$$*/pyclone_density.pdf --plot_type density --burnin 1000 --min_cluster_size 2 --max_clusters 5")

pyclone/%/pyclone_cellular_coordinates.pdf : pyclone/%/trace/alpha.tsv.bz2
	$$(call RUN,-s 4G -m 6G,"export DISPLAY='localhost:10.0' && \
							 source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
							 PyClone plot_loci --config_file pyclone/$$*/config.yaml --plot_file pyclone/$$*/pyclone_cellular_coordinates.pdf --plot_type parallel_coordinates --burnin 1000 --min_cluster_size 2 --max_clusters 5")
							 
pyclone/%/pyclone_matrix.pdf : pyclone/%/trace/alpha.tsv.bz2
	$$(call RUN,-s 4G -m 6G,"export DISPLAY='localhost:10.0' && \
							 source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
							 PyClone plot_loci --config_file pyclone/$$*/config.yaml --plot_file pyclone/$$*/pyclone_matrix.pdf --plot_type similarity_matrix --burnin 1000 --min_cluster_size 2 --max_clusters 5")

pyclone/%/pyclone_vaf_coordinates.pdf : pyclone/%/trace/alpha.tsv.bz2
	$$(call RUN,-s 4G -m 6G,"export DISPLAY='localhost:10.0' && \
							 source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
							 PyClone plot_loci --config_file pyclone/$$*/config.yaml --plot_file pyclone/$$*/pyclone_vaf_coordinates.pdf --plot_type vaf_parallel_coordinates --burnin 1000 --min_cluster_size 2 --max_clusters 5")

pyclone/%/pyclone_scatter.pdf : pyclone/%/trace/alpha.tsv.bz2
	$$(call RUN,-s 4G -m 6G,"export DISPLAY='localhost:10.0' && \
							 source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
							 PyClone plot_loci --config_file pyclone/$$*/config.yaml --plot_file pyclone/$$*/pyclone_scatter.pdf --plot_type vaf_scatter --burnin 1000 --min_cluster_size 2 --max_clusters 5")

endef
$(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call run-pyclone,$(sample))))
