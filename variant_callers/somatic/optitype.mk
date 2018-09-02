include modules/Makefile.inc

LOGDIR ?= log/hla_optitype.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: hla_optitype

VPATH = bam

hla_optitype : $(foreach sample,$(SAMPLES),hla_optitype/$(sample).taskcomplete)

define hla-optitype
hla_optitype/%.bam : %.bam
	$$(call RUN,-n 4 -s 12G -m 16G,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate \
								 	/home/${USER}/share/usr/anaconda-envs/optitype && \
								 	$(SAMTOOLS2) view -@ 4 $$< 6 -b -o hla_optitype/$$*.bam")

hla_optitype/%_1.fastq hla_optitype/%_2.fastq : hla_optitype/%.bam
	$$(call RUN,-n 4 -s 12G -m 16G,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate \
								 	/home/${USER}/share/usr/anaconda-envs/optitype && \
								 	$(SAMTOOLS2) sort -T $$(<D)/$$* -O bam -n -@ 4 -m 6G $$< | $(SAMTOOLS2) fastq -f 1 -1 > hla_optitype/$$*.1.fastq -2 > hla_optitype/$$*.2.fastq -")

hla_optitype/%_1_razers3.bam hla_optitype/%_2_razers3.bam : hla_optitype/%_1.fastq hla_optitype/%_2.fastq
	$$(call RUN,-n 4 -s 12G -m 16G,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate \
								 	/home/${USER}/share/usr/anaconda-envs/optitype && \
								 	razers3 -i 95 -m 1 -dr 0 -o hla_optitype/$$*_1_razers3.bam /home/${USER}/share/usr/anaconda-envs/optitype/share/optitype-1.3.2-1/data/hla_reference_dna.fasta hla_optitype/$$*_1.fastq && \
								 	razers3 -i 95 -m 1 -dr 0 -o hla_optitype/$$*_2_razers3.bam /home/${USER}/share/usr/anaconda-envs/optitype/share/optitype-1.3.2-1/data/hla_reference_dna.fasta hla_optitype/$$*_2.fastq")

hla_optitype/%.taskcomplete : hla_optitype/%_1_razers3.bam hla_optitype/%_2_razers3.bam
	$$(call RUN,-s 1G -m 1G,"touch hla_optitype/$$*.taskcomplete")
endef
$(foreach pair,$(SAMPLES),\
		$(eval $(call hla-optitype,$sample)))
