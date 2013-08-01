# Run VarScan on tumour-normal matched pairs
# Detect point mutations
##### DEFAULTS ######


REF ?= hg19
LOGDIR = log/varscan.$(NOW)
SAMPLE_PAIR_FILE ?= sample_pairs.txt
SAMPLE_FILE ?= samples.txt
TUMOR_SAMPLES ?= $(shell cut -f 1 $(SAMPLE_PAIR_FILE))
NORMAL_SAMPLES ?= $(shell cut -f 2 $(SAMPLE_PAIR_FILE))
SAMPLES ?= $(shell cat $(SAMPLE_FILE))

$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval normal_lookup.$(word $i,$(TUMOR_SAMPLES)) := $(word $i,$(NORMAL_SAMPLES))))
$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval tumor_lookup.$(word $i,$(NORMAL_SAMPLES)) := $(word $i,$(TUMOR_SAMPLES))))

SPLIT_CHR ?= true

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc

VARSCAN_JAR = $(JARDIR)/VarScan.v2.3.5.jar
VARSCAN = $(JAVA) -Xmx8G -jar $(VARSCAN_JAR)
SEGMENTCNV = $(HOME)/share/scripts/segmentCNV.R

VPATH ?= bam
SPLIT_CHR ?= true

.DELETE_ON_ERROR:

.SECONDARY: 

.PHONY: all vcfs copycalls segments cnv

all : vcfs cnv

cnv : copycalls segments
vcfs : $(foreach tumor,$(TUMOR_SAMPLES),vcf/$(tumor)_$(normal_lookup.$(tumor)).varscan_indels.vcf) $(foreach tumor,$(TUMOR_SAMPLES),vcf/$(tumor)_$(normal_lookup.$(tumor)).varscan_snps.vcf)
copycalls : $(foreach tumor,$(TUMOR_SAMPLES),varscan/copycall/$(tumor)_$(normal_lookup.$(tumor)).copycall)
segments : $(foreach tumor,$(TUMOR_SAMPLES),varscan/segment/$(tumor)_$(normal_lookup.$(tumor)).segment.txt)
	

ifeq ($(SPLIT_CHR),true)
#varscan/chr_vcf/$1_$2.$3.indels.vcf varscan/chr_vcf/$1_$2.$3.snps.vcf : $1.bam $2.bam
define varscan-somatic-tumor-normal-chr
varscan/chr_vcf/$1_$2.$3.varscan_timestamp : $1.bam $2.bam
	$$(call LSCRIPT_MEM,9G,12G,"tmpfile=`mktemp`; \
	mkfifo $$$$tmpfile.$1; mkfifo $$$$tmpfile.$2; \
	$$(SAMTOOLS) mpileup -l $3 -q 1 -f $$(REF_FASTA) $$< > $$$$tmpfile.$1 & \
	$$(SAMTOOLS) mpileup -l $3 -q 1 -f $$(REF_FASTA) $$(word 2,$$^) > $$$$tmpfile.$2 & \
	$$(VARSCAN) somatic $$$$tmpfile.$1 $$$$tmpfile.$2 --output-vcf --output-indel varscan/chr_vcf/$1_$2.$3.indels.vcf --output-snp varscan/chr_vcf/$1_$2.$3.snps.vcf &> $$(LOG) && touch $$@)")

varscan/chr_vcf/$1_$2.$3.indels.vcf : varscan/chr_vcf/$1_$2.$3.varscan_timestamp

varscan/chr_vcf/$1_$2.$3.snps.vcf : varscan/chr_vcf/$1_$2.$3.varscan_timestamp

endef
$(foreach chr,$(CHROMOSOMES),$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call varscan-somatic-tumor-normal-chr,$(tumor),$(normal_lookup.$(tumor)),$(chr)))))


define varscan-ped-tumor-normal
vcf/$1_$2.varscan_%.vcf : $$(foreach chr,$$(CHROMOSOMES),varscan/chr_vcf/$1_$2.$$(chr).%.vcf)
	$$(INIT) echo "##PEDIGREE=<Derived=$1,Original=$2>" > $$@ && \
	grep '^#' $$< >> $$@ && \
	cat $$^ | grep -v '^#' | $$(VCF_SORT) $$(REF_DICT) -  >> $$@
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call varscan-ped-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))

else # no splitting by chr

define varscan-somatic-tumor-normal
varscan/vcf/$1_$2.varscan_timestamp : $1.bam $2.bam
	$$(call LSCRIPT_MEM,9G,12G,"tmpfile=`mktemp`; \
	mkfifo $$$$tmpfile.$1; mkfifo $$$$tmpfile.$2; \
	$$(SAMTOOLS) mpileup -q 1 -f $$(REF_FASTA) $$< > $$$$tmpfile.$1 2> /dev/null & \
	$$(SAMTOOLS) mpileup -q 1 -f $$(REF_FASTA) $$(word 2,$$^) > $$$$tmpfile.$2 2> /dev/null & \
	$$(VARSCAN) somatic $$$$tmpfile.$2 $$$$tmpfile.$1 --output-vcf --output-indel varscan/vcf/$1_$2.indels.vcf --output-snp varscan/vcf/$1_$2.snps.vcf &> $$(LOG)")

varscan/chr_vcf/$1_$2.indels.vcf : varscan/chr_vcf/$1_$2.varscan_timestamp

varscan/chr_vcf/$1_$2.snps.vcf : varscan/chr_vcf/$1_$2.varscan_timestamp

vcf/$1_$2.varscan_%.vcf : varscan/vcf/$1_$2.%.vcf
	$$(INIT) echo "##PEDIGREE=<Derived=$1,Original=$2>" > $$@ && cat $$< >> $$@
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call varscan-somatic-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))
endif

define varscan-copynum-tumor-normal
varscan/copynum/$1_$2.copynumber :  $1.bam $2.bam
	$$(call LSCRIPT_MEM,9G,12G,"$$(SAMTOOLS) mpileup -q 1 -f $$(REF_FASTA) $$(word 2,$$^) $$< | awk 'NF == 9 { print }' |  $$(VARSCAN) copynumber - $$(basename $$@) --mpileup 1 &> $$(LOG)")
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call varscan-copynum-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))

varscan/copycall/%.copycall : varscan/copynum/%.copynumber
	$(call LSCRIPT_MEM,9G,12G,"$(VARSCAN) copyCaller $< --output-file $@ &> $(LOG)")

varscan/segment/%.segment.txt : varscan/copycall/%.copycall
	$(call LSCRIPT_MEM,4G,6G,"$(RSCRIPT) $(SEGMENTCNV) --prefix=varscan/segment/$* $< &> $(LOG)")

include ~/share/modules/gatk.mk
