include modules/Makefile.inc

LOGDIR ?= log/pyclone.$(NOW)
PHONY += pyclone

pyclone : $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample).taskcomplete)

define make-input-pyclone
pyclone/%.timestamp : sufam/%.tsv
	$$(call RUN,-c -s 4G -m 6G,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
								if [ ! -d pyclone/$$* ]; then mkdir pyclone/$$*; fi && \
								$(RSCRIPT) modules/clonality/pyclone_make_input.R --file_name $$< && \
								touch pyclone/$$*.timestamp")
endef
$(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call make-input-pyclone,$(sample))))

define build-mutations-file
pyclone/$2/$1.yaml : pyclone/$2/$1.tsv
	$$(call RUN,-c -s 4G -m 6G,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
								PyClone build_mutations_file --in_file $$< --out_file $$<< --prior parental_copy_number")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call build-mutations-file,$(tumor.$(pair)),$(normal.$(pair)))))
		
define complete-task
pyclone/%.taskcomplete : pyclone/%.timestamp
	$$(call RUN,-c -s 1G -m 1G,"touch pyclone/$$*.taskcomplete")
endef
$(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call complete-task,$(sample))))
