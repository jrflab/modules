# Run gatk depth of coverage on bam files

include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

SAMPLE_FILE ?= samples.txt
SAMPLES ?= $(shell cat $(SAMPLE_FILE))

LOGDIR = log/read_depth.$(NOW)

VPATH ?= bam

SPLIT_CHR ?= true

.PHONY: all

all: $(foreach sample,$(SAMPLES),gatk/read_depth/$(sample).grp)

ifeq ($(SPLIT_CHR),true)
define read-depth-chr
gatk/chr_read_depth/%.$1.read_depth : %.bam
	$$(call INIT_MEM,8G,12G) $$(call GATK_MEM,7G) -T DepthOfCoverage -R $$(REF_FASTA) -L $1 -o $$@ -I $$< &> $$(LOG)
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call read-depth-chr,$(chr))))
gatk/read_depth/% : $(foreach chr,$(CHROMOSOMES),gatk/chr_read_depth/%.$(chr).read_depth)
	$(INIT) head -1 $< > $@ && sed '1d' $^ >> $@
else
gatk/read_depth/%.read_depth : %.bam
	$(call INIT_MEM,8G,12G) $(call GATK_MEM,7G) -T DepthOfCoverage -R $(REF_FASTA) -o $@ -I $< &> $(LOG)
endif
