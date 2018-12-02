include modules/Makefile.inc

LOGDIR ?= log/plot_pyclone.$(NOW)
PHONY += pyclone

plot_pyclone : $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/plots/by_loci_density.pdf) $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/plots/all_loci_density.pdf) $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/plots/all_loci_scatter.pdf) $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/summary.tsv) $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/plots/all_loci_scatter_filtered.pdf)

define plot-pyclone
pyclone/%/plots/by_loci_density.pdf : pyclone/%/pyclone.tsv
	$$(call RUN,-s 4G -m 6G -w 7200,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.6/bin/activate /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.6/ && \
							 $(RSCRIPT) modules/clonality/pyclonelocidensity.R --sample_name $$* --burnin 5000")

pyclone/%/plots/all_loci_density.pdf : pyclone/%/pyclone.tsv
	$$(call RUN,-s 4G -m 6G -w 7200,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.6/bin/activate /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.6/ && \
							 $(RSCRIPT) modules/clonality/pyclonealldensity.R --sample_name $$* --burnin 5000")
							 
pyclone/%/plots/all_loci_scatter.pdf pyclone/%/summary.tsv : pyclone/%/pyclone.tsv sufam/%.tsv
	$$(call RUN,-s 4G -m 6G -w 7200,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.6/bin/activate /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.6/ && \
							 $(RSCRIPT) modules/clonality/pyclonelociscatter.R --sample_name $$* --burnin 5000")
							 
pyclone/%/plots/all_loci_scatter_filtered.pdf : pyclone/%/summary.tsv
	$$(call RUN,-s 4G -m 6G -w 7200,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.6/bin/activate /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.6/ && \
							 $(RSCRIPT) modules/clonality/pyclonelociscatterupdated.R --sample_name $$*")

endef
$(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call plot-pyclone,$(sample))))
