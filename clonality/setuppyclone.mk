include modules/Makefile.inc

LOGDIR ?= log/setup_pyclone.$(NOW)
PHONY += pyclone

setup_pyclone : $(foreach set,$(SAMPLE_SETS),pyclone/$(set)/config.yaml)

define make-input-pyclone
pyclone/%/config.yaml : sufam/%.tsv
	$$(call RUN, -s 4G -m 6G,"mkdir -p pyclone/$$(*) && \
							  $(RSCRIPT) modules/clonality/tsvforpyclone.R --sample_set $$(*) --normal_samples $(NORMAL_SAMPLES) && \
							  $(RSCRIPT) modules/clonality/pycloneconfig.R --sample_set $$(*) --normal_samples $(NORMAL_SAMPLES) && \
							  source /home/$(USER)/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/$(USER)/share/usr/anaconda-envs/PyClone-0.13.1 && \
							  for i in $(ls pyclone/$$(*)/*.tsv); do j=${i%.*}; PyClone build_mutations_file --in_file $i --out_file $j.yaml --prior total_copy_number; done;")

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call make-input-pyclone,$(set))))

#define make-config-yaml
#pyclone/%/config.yaml : $(wildcard $(foreach set,$(SAMPLE_SETS),sufam/$(set).tsv))
#	$$(call RUN,-c -s 4G -m 6G,"if [ ! -d pyclone/$$(*) ]; then mkdir pyclone/$$(*); fi && \
#								")
#endef
#$(foreach sample,$(NORMAL_SAMPLES),\
#		$(eval $(call make-config-yaml,$(sample))))
		
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
