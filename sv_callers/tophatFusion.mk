include modules/Makefile.inc
include modules/variant_callers/gatk.inc

TOPHAT := $(HOME)/share/usr/bin/tophat2
TOPHAT_OPTS := --no-coverage-search --fusion-ignore-chromosomes MT --fusion-search --keep-fasta-order

ifeq ($(BAM_PHRED64),true)
	TOPHAT_OPTS += --solexa1.3-quals
endif

.SECONDARY:
.DELETEONERROR:
.PHONY: all

all : $(foreach sample,$(SAMPLES),tophat/$(sample)/fusions.out)

tophat/%/fusions.out : fastq/%.1.fastq.gz fastq/%.2.fastq.gz tophat/ins_size/%.insert_size.txt
	DIST_OPTS=`perl -e '$$text = do {local $$/ ; <>}; $$text =~ m/Read length: mean (\d+).*\nRead span: mean (\d+).*STD=(\d+)/; print "--mate-inner-dist " . ($$2 - $$1 * 2) . " --mate-std-dev $$3";' $(word 3,$^)`; \
	$(call LSCRIPT_NAMED_PARALLEL_MEM,$*_tophat,4,6G,10G,"$(TOPHAT) $(TOPHAT_OPTS) -p 4 $$DIST_OPTS -o $(@D) $(BOWTIE_REF) $(<) $(word 2,$(^))")

tophat/ins_size/%.insert_size.txt : bam/%.bam
	$(call LSCRIPT,"$(SAMTOOLS) view $< | $(GET_INSERT_SIZE) - > $@")

tophat/fusions/%.fusions.ft.txt : tophat/%/fusions.out
	awk '$5 > 100 { print }' $< 

