# This module gets the unmapped reads from the samples bam files and maps them
# against a virus database. The coverage of the virus genomes is determined #
# with BEDtools.
# Author: inodb

##### MAKE INCLUDES #####
include modules/Makefile.inc

LOGDIR = log/virus_detection_bowtie2.$(NOW)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: virus_detection_bowtie2

VIRUS_DB ?= /ifs/e63data/reis-filho/reference/viral.1.1.genomic.renamed.fna
VIRUS_DB_NAME ?= ncbi-viral

virus_detection_bowtie2 : $(foreach sample,$(SAMPLES),virus/bowtie2/default/$(VIRUS_DB_NAME)/$(sample).bam.bedcov.summary.txt)

VIRUS_MATE = <(grep '/$2$$' -A3 --no-group-separator $1)

PYTHON_ENV ?= /ifs/e63data/reis-filho/usr/anaconda-envs/pyenv27-chasm

virus/fastq/%.unmapped.il.fq virus/fastq/%.unmapped.single.fq: bam/%.bam
	$(call RUN,-s 2G -m 3G,"$(MKDIR) $(@D);  samtools view -uf 4 $< | samtools bam2fq - -s $(@D)/$*.unmapped.single.fq > $(@D)/$*.unmapped.il.fastq")

virus/bowtie2/default/$(VIRUS_DB_NAME)/%.bam: virus/fastq/%.unmapped.il.fq virus/fastq/%.unmapped.single.fq
	$(call RUN,-s 2G -m 3G,"$(MKDIR) $(@D); bowtie2 -x $(VIRUS_DB) -U $(<<) -1 $(call VIRUS_MATE,$<,1) -2 $(call VIRUS_MATE,$<,2) | samtools view -bS - | samtools sort - -O bam -T `mktemp` > $@")

virus/bowtie2/default/$(VIRUS_DB_NAME)/%.bam.bedcov.txt: virus/bowtie2/default/$(VIRUS_DB_NAME)/%.bam
	$(call RUN,-s 2G -m 3G,"genomeCoverageBed -ibam $< > $@")

virus/bowtie2/default/$(VIRUS_DB_NAME)/%.bam.bedcov.summary.txt: virus/bowtie2/default/$(VIRUS_DB_NAME)/%.bam.bedcov.txt
	$(INIT) unset PYTHONPATH && \
	source $(PYTHON_ENV)/bin/activate $(PYTHON_ENV) && \
	echo -e '\
import pandas as pd\n\
cov = pd.read_csv("$<", sep="\t", names=["ref","cov","count","len","ratio"])\n\
cov.groupby("ref").apply(lambda x: pd.Series({"mean_cov":(x["cov"] * x["count"] / x["len"]).sum(),"perc_cov":(100 - x[x["cov"] == 0]["ratio"].iloc[0]*100 if len(x[x["cov"]==0]) > 0 else 1)})).sort("perc_cov", ascending=False).to_csv("$@", sep="\t")' | python
