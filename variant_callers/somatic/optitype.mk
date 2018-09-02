include modules/Makefile.inc

LOGDIR ?= log/hla_optitype.$(NOW)
PHONY += hla_optitype

hla_optitype : $(foreach sample,$(SAMPLES),hla_optitype/$(sample)/$(sample).bam)

define hla-optitype
hla_optitype/%/%fastq : bam/%.bam
	$$(call RUN,-n 8 -s 12G -m 24G, "source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate \
									 /home/${USER}/share/usr/anaconda-envs/optitype && \
								 	 if [ ! -d hla_optitype/$* ]; then mkdir hla_optitype/$*; fi && \
								 	 razers3 -i 95 -m 1 -dr 0 -o hla_optitype/$*/$*.bam /home/${USER}/share/usr/anaconda-envs/optitype/share/optitype-1.3.2-1/data/hla_reference_dna.fasta hla_optitype/$*/$*.fastq")
endef
$(foreach pair,$(SAMPLES),\
		$(eval $(call hla-optitype,$(sample),$(sample))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
