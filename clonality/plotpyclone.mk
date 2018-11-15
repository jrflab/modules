include modules/Makefile.inc

LOGDIR ?= log/plot_pyclone.$(NOW)
PHONY += pyclone

plot_pyclone : $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/plots/by_loci_density.pdf) $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/plots/all_loci_density.pdf) $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/plots/all_loci_scatter.pdf)

#$(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/pyclone_loci_coordinates.pdf) $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/pyclone_loci_matrix.pdf) $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/pyclone_vaf_coordinates.pdf) $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/pyclone_cluster_density.pdf) $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/parallel_cluster_coordinates.pdf) $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/pyclone_cluster_scatter.pdf)

define plot-pyclone
pyclone/%/plots/by_loci_density.pdf : pyclone/%/pyclone.tsv
	$$(call RUN,-s 4G -m 6G,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.6/bin/activate /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.6/ && \
							 $(RSCRIPT) modules/clonality/pyclonelocidensity.R --sample_name $$* --burnin 5000")

pyclone/%/plots/all_loci_density.pdf : pyclone/%/pyclone.tsv
	$$(call RUN,-s 4G -m 6G,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.6/bin/activate /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.6/ && \
							 $(RSCRIPT) modules/clonality/pyclonealldensity.R --sample_name $$* --burnin 5000")
							 
pyclone/%/plots/all_loci_scatter.pdf : pyclone/%/pyclone.tsv
	$$(call RUN,-s 4G -m 6G,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.6/bin/activate /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.6/ && \
							 $(RSCRIPT) modules/clonality/pyclonelociscatter.R --sample_name $$* --burnin 5000")


#pyclone/%/pyclone_loci_coordinates.pdf : pyclone/%/pyclone.tsv
#	$$(call RUN,-s 4G -m 6G,"export DISPLAY=localhost:10.0 && \
#							 source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
#							 PyClone plot_loci --config_file pyclone/$$*/config.yaml --plot_file pyclone/$$*/pyclone_loci_coordinates.pdf --plot_type parallel_coordinates --burnin 1000 --min_cluster_size 2 --max_clusters 5")
#							 
#pyclone/%/pyclone_loci_matrix.pdf : pyclone/%/pyclone.tsv
#	$$(call RUN,-s 4G -m 6G,"export DISPLAY=localhost:10.0 && \
#							 source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
#							 PyClone plot_loci --config_file pyclone/$$*/config.yaml --plot_file pyclone/$$*/pyclone_loci_matrix.pdf --plot_type similarity_matrix --burnin 1000 --min_cluster_size 2 --max_clusters 5")

#pyclone/%/pyclone_vaf_coordinates.pdf : pyclone/%/pyclone.tsv
#	$$(call RUN,-s 4G -m 6G,"export DISPLAY=localhost:10.0 && \
#							 source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
#							 PyClone plot_loci --config_file pyclone/$$*/config.yaml --plot_file pyclone/$$*/pyclone_vaf_coordinates.pdf --plot_type vaf_parallel_coordinates --burnin 1000 --min_cluster_size 2 --max_clusters 5")
#
#pyclone/%/pyclone_loci_scatter.pdf : pyclone/%/pyclone.tsv
#	$$(call RUN,-s 4G -m 6G,"export DISPLAY=localhost:10.0 && \
#							 source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
#							 PyClone plot_loci --config_file pyclone/$$*/config.yaml --plot_file pyclone/$$*/pyclone_loci_scatter.pdf --plot_type vaf_scatter --burnin 1000 --min_cluster_size 2 --max_clusters 5")
#
#pyclone/%/pyclone_cluster_density.pdf : pyclone/%/pyclone.tsv
#	$$(call RUN,-s 4G -m 6G,"export DISPLAY=localhost:10.0 && \
#							 source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
#							 PyClone plot_clusters --config_file pyclone/$$*/config.yaml --plot_file pyclone/$$*/pyclone_cluster_density.pdf --plot_type density --burnin 1000 --min_cluster_size 2 --max_clusters 5")
#
#pyclone/%/parallel_cluster_coordinates.pdf : pyclone/%/pyclone.tsv
#	$$(call RUN,-s 4G -m 6G,"export DISPLAY=localhost:10.0 && \
#							 source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
#							 PyClone plot_clusters --config_file pyclone/$$*/config.yaml --plot_file pyclone/$$*/parallel_cluster_coordinates.pdf --plot_type parallel_coordinates --burnin 1000 --min_cluster_size 2 --max_clusters 5")
#
#pyclone/%/pyclone_cluster_scatter.pdf : pyclone/%/pyclone.tsv
#	$$(call RUN,-s 4G -m 6G,"export DISPLAY=localhost:10.0 && \
#							 source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
#							 PyClone plot_clusters --config_file pyclone/$$*/config.yaml --plot_file pyclone/$$*/pyclone_cluster_scatter.pdf --plot_type scatter --burnin 1000 --min_cluster_size 2 --max_clusters 5")

endef
$(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call plot-pyclone,$(sample))))
