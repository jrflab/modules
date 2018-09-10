include modules/Makefile.inc

LOGDIR ?= log/sufam_multisample.$(NOW)
PHONY += sufam

sufam_multisample : $(foreach sample,$(NORMAL_SAMPLES),sufam/$(sample).tsv)

define combine-samples
sufam/%.txt : summary/tsv/mutation_summary.tsv
	$$(call RUN,-c -s 4G -m 6G,"source /home/${USER}/share/usr/opt/miniconda/bin/activate /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/ && \
								$(RSCRIPT) modules/variant_callers/combineSamples.R --patient $$*")

sufam/%.tsv : sufam/%.txt
	$$(call RUN,-c -s 4G -m 6G -v /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5,"source /home/${USER}/share/usr/opt/miniconda/bin/activate /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/ && \
								$(RSCRIPT) modules/variant_callers/updateSamples.R --patient $$*")
	
endef
$(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call combine-samples,$(sample))))
