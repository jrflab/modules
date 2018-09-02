include modules/Makefile.inc

LOGDIR ?= log/hla_optitype.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: hla_optitype

VPATH = bam

hla_optitype : $(foreach sample,$(SAMPLES),hla_optitype/$(sample).bam)

define hla-optitype
hla_optitype/%.bam : %.bam
	$$(call RUN,-n 4 -s 12G -m 16G,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate \
								 	/home/${USER}/share/usr/anaconda-envs/optitype && \
								 	$(SAMTOOLS2) view -f 0x04 -h -@ 4 -b $$< -o hla_optitype/$$*.fastq")
endef
$(foreach pair,$(SAMPLES),\
		$(eval $(call hla-optitype,$sample)))


# hla_optitype/$$*.bam /home/${USER}/share/usr/anaconda-envs/optitype/share/optitype-1.3.2-1/data/hla_reference_dna.fasta hla_optitype/$$*.fastq")
