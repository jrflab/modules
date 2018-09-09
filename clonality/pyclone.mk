include modules/Makefile.inc

LOGDIR ?= log/pyclone.$(NOW)
PHONY += pyclone

pyclone : $(foreach sample,$(NORMAL_SAMPLES),sufam/$(sample)/$(sample).taskcomplete)

define make-input-pyclone
pyclone/%/%.taskcomplete : sufam/%.tsv
	$$(call RUN,-c -s 4G -m 6G,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate && \
								$(RSCRIPT) modules/clonality/pyclone_make_input.R --file_name $$< && \
								touch $$*")

endef
$(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call make-input-pyclone,$(sample))))
