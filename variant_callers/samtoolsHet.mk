# Run samtools to detect heterozygous positions
##### DEFAULTS ######
include ~/share/modules/Makefile.inc
include ~/share/modules/variant_callers/gatk.inc

LOGDIR ?= log/samtools_het.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all

all : $(foreach s,$(SAMPLES),vcf/$s.het_snp.vcf vcf/$s.het_snp.vcf.idx)

include ~/share/modules/vcf_tools/vcftools.mk

define hetsnp-chr
chr_vcf/%.$1.het_snp.vcf : bam/%.bam
	$$(call LSCRIPT_MEM,6G,8G,"$$(SAMTOOLS2) mpileup -r $1 -f $$(REF_FASTA) -g -I $$< | $$(BCFTOOLS2) call -c | $$(BCFTOOLS2) view -g het | $$(VCFUTILS) varFilter -d 10 -a 5 - > $$@")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call hetsnp-chr,$(chr))))

vcf/%.het_snp.vcf : $(foreach chr,$(CHROMOSOMES),chr_vcf/%.$(chr).het_snp.vcf)
	$(INIT) grep -P '^#' $< > $@ && sed '/^#/d' $^ >> $@

