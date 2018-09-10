include modules/Makefile.inc

LOGDIR ?= log/run_pyclone.$(NOW)
PHONY += pyclone

run_pyclone : $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/trace/alpha.tsv.bz2)

define run-pyclone
pyclone/%/trace/alpha.tsv.bz2 : pyclone/%/config.yaml
	$$(call RUN,-c -s 4G -m 6G,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
								PyClone run_analysis --config_file pyclone/$$*/config.yaml --seed 210783")
endef
$(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call run-pyclone,$(sample))))
