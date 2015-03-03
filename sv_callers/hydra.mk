# run hydra

LOGDIR = log/hydra.$(NOW)

HYDRA = $(HOME)/share/usr/bin/hydra
override HYDRA_OPTS ?= -mld 500 -mn 1500
BAM_TO_FASTQ = $(HOME)/share/usr/bin/bamToFastq
BAM_TO_BED = /opt/common/bedtools/bedtools-2.17.0/bin/bamToBed
DEDUP_DISCORDANTS = $(HOME)/share/usr/bin/dedupDiscordants.py
PAIR_DISCORDANTS = $(HOME)/share/usr/bin/pairDiscordants.py

include modules/Makefile.inc
include modules/variant_callers/gatk.inc

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all

all : $(foreach sample,$(SAMPLES),hydra/breaks/$(sample).breaks)


#hydra/disc_fastq/%.disc.1.fastq.gz hydra/disc_fastq/%.disc.2.fastq.gz : bam/%.bam
#$(INIT) $(SAMTOOLS) view -uF 2 $< | $(BAM_TO_FASTQ) -bam stdin -fq1 >( gzip -c > hydra/disc_fastq/$*.disc.1.fastq.gz) -fq2  >( gzip -c > hydra/disc_fastq/$*.disc.2.fastq.gz)

hydra/bam/%.disc.bam : bam/%.bam
	$(call LSCRIPT,"$(SAMTOOLS) view -bF 2 $< > $@")

hydra/bed/%.disc.bedpe : hydra/bam/%.disc.bam
	$(call LSCRIPT,"$(BAM_TO_BED) -i $< -tag NM | $(PAIR_DISCORDANTS) -i stdin -m hydra -z 800 > $@")

hydra/bed/%.disc.dedup.bedpe : hydra/bed/%.disc.bedpe
	$(call LSCRIPT,"$(DEDUP_DISCORDANTS) -i $< -s 3 > $@")

hydra/breaks/%.breaks : hydra/bed/%.disc.dedup.bedpe
	$(call LSCRIPT,"$(HYDRA) -in $< -out $@ $(HYDRA_OPTS)")
