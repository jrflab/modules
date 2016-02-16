# run lumpy

LOGDIR = log/lumpy.$(NOW)

include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LUMPY_DIR = $(HOME)/share/usr/lumpy-sv
LUMPY_SCRIPTS_DIR = $(LUMPY_DIR)/scripts
EXTRACT_DISCORDANT = $(PYTHON) $(LUMPY_SCRIPTS_DIR)/extractSplitReads_BwaMem
LUMPY = $(LUMPY_DIR)/bin/lumpy
LUMPYEXPRESS = $(LUMPY_DIR)/bin/lumpyexpress
LUMPYEXPRESS_OPTS = -K $(LUMPY_DIR)/bin/lumpyexpress.config
LUMPY_HISTO = $(PERL) $(LUMPY_SCRIPTS_DIR)/pairend_distro.pl
LUMPY_UNMAPPED_TO_FASTQ = $(PERL) $(LUMPY_SCRIPTS_DIR)/split_unmapped_to_fasta.pl
LUMPY_UNMAPPED_TO_FASTQ_OPTS = -b 20

SAMBLASTER = $(HOME)/share/usr/bin/samblaster

# deprecated, using lumpyexpress instead
#LUMPY_OPTS = -g $(REF_FASTA) -tt 1e-3 -mw 4
#LUMPY_PE_PARAMS = min_non_overlap:150$(,)discordant_z:4$(,)back_distance:20$(,)weight:1$(,)id:1$(,)min_mapping_threshold:1
#LUMPY_SR_PARAMS = back_distance:20$(,)weight:1$(,)id:1$(,)min_mapping_threshold:1

BWASW_OPTS = -H

ANNOVAR_PROTOCOL = refGene$(,)cytoBand$(,)genomicSuperDups
ANNOVAR_OPERATION = g$(,)r$(,)r

LUMPY_SUFFIX = sv_som_ft.pass.$(ANNOVAR_REF)_multianno

.SECONDARY:
.DELETE_ON_ERROR:


PHONY += lumpy lumpy_vcf lumpy_tables
lumpy : lumpy_vcf lumpy_tables
ifdef SAMPLE_PAIRS
lumpy_vcf : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).lumpy.$(LUMPY_SUFFIX).vcf)
lumpy_tables : $(foreach pair,$(SAMPLE_PAIRS),tables/$(pair).lumpy.$(LUMPY_SUFFIX).tab.txt) alltables/allTN.lumpy.$(LUMPY_SUFFIX).tab.txt
else
lumpy_vcf : $(foreach sample,$(SAMPLES),vcf/$(sample).lumpy.$(LUMPY_SUFFIX).vcf)
lumpy_tables : $(foreach sample,$(SAMPLES),vcf/$(sample).lumpy.$(LUMPY_SUFFIX).tab.txt) alltables/all.lumpy.$(LUMPY_SUFFIX).tab.txt
endif

lumpy/bam/%.split.bam : bam/%.bam
	$(call LSCRIPT_MEM,7G,9G,"$(SAMTOOLS2) view -h $< | $(EXTRACT_DISCORDANT) -i stdin | $(SAMTOOLS2) view -b - > $@")

lumpy/bam/%.disc.bam : bam/%.bam
	$(call LSCRIPT_MEM,7G,9G,"$(SAMTOOLS2) view -b -F 1294 $< > $@")

lumpy/metrics/%.read_len : bam/%.bam
	$(INIT) $(SAMTOOLS) view $< | tail -n+100000 | head -1 | awk '{ print length($$10) }' > $@

lumpy/metrics/%.histo lumpy/metrics/%.histo.txt: bam/%.bam lumpy/metrics/%.read_len
	READ_LEN=`cat $(word 2,$^)`; \
	$(call LSCRIPT,"$(SAMTOOLS) view $< | tail -n+100000 | $(LUMPY_HISTO) -rl $$READ_LEN -X 4 -N 10000 -o lumpy/metrics/$*.histo > lumpy/metrics/$*.histo.txt")

lumpy/vcf/%.lumpy.vcf : bam/%.bam lumpy/bam/%.split.bam lumpy/bam/%.disc.bam
	$(call LSCRIPT_MEM,15G,20G,"$(LUMPYEXPRESS) $(LUMPYEXPRESS_OPTS) -B $<$ -S $(<<) -D $(<<<) -o $@")

ifdef SAMPLE_PAIRS
define lumpy-tumor-normal
lumpy/vcf/$1_$2.lumpy.vcf : bam/$1.bam bam/$2.bam lumpy/bam/$1.split.bam lumpy/bam/$2.split.bam lumpy/bam/$1.disc.bam lumpy/bam/$2.disc.bam
	$$(call LSCRIPT_MEM,30G,60G,"$$(LUMPYEXPRESS) $$(LUMPYEXPRESS_OPTS) -B $$<$$(,)$$(<<) -S $$(<<<)$$(,)$$(<<<<) -D $$(word 5,$$^)$$(,)$$(word 6,$$^) -o $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call lumpy-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif

SORT_VCF = $(PERL) $(HOME)/share/usr/bin/vcfsorter.pl
vcf/%.lumpy.vcf : lumpy/vcf/%.lumpy.vcf
	$(call LSCRIPT,"$(SORT_VCF) $(REF_DICT) $< > $@")

include modules/vcf_tools/vcftools.mk

.PHONY: $(PHONY)

