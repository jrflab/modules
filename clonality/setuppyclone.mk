include modules/Makefile.inc

LOGDIR ?= log/setup_pyclone.$(NOW)
PHONY += pyclone

setup_pyclone : $(foreach pair,$(SAMPLE_PAIRS),pyclone/$(normal.$(pair))/$(tumor.$(pair)).yaml)
configure_pyclone : $(foreach sample,$(NORMAL_SAMPLES),pyclone/$(sample)/config.yaml)

define make-input-pyclone
pyclone/$2/$1.tsv pyclone/$2/$1.yaml : $(wildcard $(foreach set,$(SAMPLE_SETS),sufam/$(set).tsv))
	$$(call RUN,-c -s 4G -m 6G -w 7200,"if [ ! -d pyclone/$2 ]; then mkdir pyclone/$2; fi && \
										$(RSCRIPT) modules/clonality/tsvforpyclone.R --sample_name $1")
	
	$$(call RUN,-c -s 4G -m 6G -w 7200,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
										while [ ! -f pyclone/$2/$1.tsv ]; do sleep 120; done && \
										PyClone build_mutations_file --in_file pyclone/$2/$1.tsv --out_file pyclone/$2/$1.yaml --prior total_copy_number")

endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call make-input-pyclone,$(tumor.$(pair)),$(normal.$(pair)))))

define make-config-yaml
pyclone/%/config.yaml : pyclone/%/
	$$(call RUN,-c -s 4G -m 6G -w 7200,"$(RSCRIPT) modules/clonality/pycloneconfig.R --sample_name $$*")
endef
$(foreach sample,$(NORMAL_SAMPLES),\
		$(eval $(call make-config-yaml,$(sample))))
