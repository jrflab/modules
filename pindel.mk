# Run pindel
# Detect indels
##### DEFAULTS ######
REF ?= hg19
LOGDIR = log/pindel.$(NOW)

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

GET_INSERT_SIZE = $(HOME)/share/usr/bin/getInsertSize.py
PINDEL = $(HOME)/share/usr/bin/pindel
PINDEL2VCF = $(HOME)/share/usr/bin/pindel2vcf

PINDEL2VCF_OPTS = -G -co 50 -ir 2 -il 3 -pr 2 -pl 3

VPATH ?= bam

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all vcfs tables

##### MAIN TARGETS ######
ANN_TYPES = eff # annotated
EFF_TYPES = silent missense nonsilent_cds

FILTER_SUFFIX := dbsnp.nsfp.eff
TABLE_SUFFIXES = $(foreach eff,$(EFF_TYPES),pindel.$(FILTER_SUFFIX).$(eff).pass.novel)

VCFS = $(foreach set,$(SAMPLE_SETS),vcf/$(set).pindel.$(FILTER_SUFFIX).vcf)
TABLES = $(foreach set,$(SAMPLE_SETS),$(foreach suff,$(TABLE_SUFFIXES),tables/$(set).$(suff).txt))
#TABLES += $(foreach suff,$(TABLE_SUFFIXES),tables/all.$(suff).txt)

all : vcfs tables
vcfs : $(VCFS)
tables : $(TABLES)

pindel/ins_size/%.insert_size.txt : bam/%.bam
	$(call LSCRIPT,"$(SAMTOOLS) view $< | $(GET_INSERT_SIZE) - > $@")

define pindel-config
pindel/config/$$(subst $$( ),_,$1).pindel_config.txt : $$(foreach sample,$1,bam/$$(sample).bam pindel/ins_size/$$(sample).insert_size.txt )
	$$(INIT) rm -f $$@; for sample in $1; do echo "bam/$$$$sample.bam " `perl -ne 'if (/^Read span: mean ([0-9]+)/) { print "$$$$1"}' pindel/ins_size/$$$$sample.insert_size.txt` " $$$$sample" >> $$@; done
endef
$(foreach i,$(shell seq 1 $(NUM_SETS)),$(eval $(call pindel-config,$(set.$i))))

define pindel-chr
pindel/%.$1.pindel_timestamp : pindel/config/%.pindel_config.txt
	$$(call LSCRIPT_PARALLEL_MEM,4,3G,6G,"$$(MKDIR) pindel/$$* && $$(PINDEL) -T 4 -f $$(REF_FASTA) -i $$< -o pindel/$$*/$1 -c $1 && touch $$@")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call pindel-chr,$(chr))))

pindel/%.pindel_timestamp : $(foreach chr,$(CHROMOSOMES),pindel/%.$(chr).pindel_timestamp)
	$(INIT) touch $@

define pindel-chr-vcf
pindel/chr_vcf/%.pindel_$1.vcf : pindel/%.$1.pindel_timestamp
	$$(INIT) $$(PINDEL2VCF) -P pindel/$$*/$1 -v $$@ -c $1 -r $$(REF_FASTA) -R $$(REF_NAME) -d $$(REF_DATE) $$(PINDEL2VCF_OPTS) &> $$(LOG)
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call pindel-chr-vcf,$(chr))))

define merge-pindel-chr
vcf/$1.pindel.vcf : $$(foreach chr,$$(CHROMOSOMES),pindel/chr_vcf/$1.pindel_$$(chr).vcf)
	$$(call LSCRIPT_MEM,4G,6G,"$$(call GATK_MEM,3G) -T CombineVariants --assumeIdenticalSamples $$(foreach i,$$^, --variant $$i) -R $$(REF_FASTA) -o $$@ &> $$(LOG)")
endef
$(foreach set,$(SAMPLE_SETS),$(eval $(call merge-pindel-chr,$(set))))

include ~/share/modules/vcftools.mk
