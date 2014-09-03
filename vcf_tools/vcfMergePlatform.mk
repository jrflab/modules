# Merge vcf files

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc
include ~/share/modules/variant_callers/gatk.inc

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all 

LOGDIR = log/vcf_merge.$(NOW)
PLATFORMS = rnaseq exonseq
platform_suffix.rnaseq := gatk_snps
platform_suffix.exonseq := mutect museq

MERGE_SUFFIXES := $(foreach platform,$(PLATFORMS),$(foreach suff,$(platform_suffix.$(platform)),$(platform)_$(suff)))
MERGE_SUFFIX := $(subst $( ),.,$(MERGE_SUFFIXES))

VCF_SUFFIX := dp_ft.dbsnp.nsfp.eff
TABLE_SUFFIX := $(VCF_SUFFIX).ft
EFF_TYPES = silent missense nonsilent_cds nonsilent

VCF_SAMPLES := $(shell seq 0 $(shell expr $(words $(MERGE_SUFFIXES)) - 1))
VCF_GEN_IDS = GT

#all : $(foreach tumor,$(TUMOR_SAMPLES),merged_vcf/$(tumor).$(MERGE_SUFFIX).$(VCF_SUFFIX).vcf)
FULL_TABLES = $(foreach eff,$(EFF_TYPES),merged_tables/all.$(MERGE_SUFFIX).$(TABLE_SUFFIX).$(eff).pass.novel.txt)
all : $(foreach eff,$(EFF_TYPES),$(foreach tumor,$(TUMOR_SAMPLES),merged_tables/$(tumor).$(MERGE_SUFFIX).$(TABLE_SUFFIX).$(eff).pass.novel.txt)) $(FULL_TABLES)
#all : $(foreach tumor,$(TUMOR_SAMPLES),merged_tables/$(tumor)_$(normal_lookup.$(tumor)).$(MERGE_SUFFIX).$(VCF_SUFFIX).txt)
#all : merged_tables/all.$(MERGE_SUFFIX).$(TABLE_SUFFIX).txt

define select-tumor-variants-platform
merged_vcf/$1.$3_%.vcf : $3/vcf/$1_$2.%.vcf
	$$(call LSCRIPT_MEM,4G,6G,"$$(call GATK_MEM,3G) -T SelectVariants \
	-R $$(REF_FASTA)  --variant $$<  -o $$@ -sn $1 &> $$(LOG)")
endef
$(foreach platform,$(PLATFORMS),$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call select-tumor-variants-platform,$(tumor),$(normal_lookup.$(tumor)),$(platform)))))

define variants-platform
merged_vcf/$1.$2_%.vcf : $2/vcf/$1.%.vcf
	cp $$< $$@
endef
$(foreach platform,$(PLATFORMS),$(foreach sample,$(SAMPLES),$(eval $(call variants-platform,$(sample),$(platform)))))

merged_vcf/%.$(MERGE_SUFFIX).vcf : $(foreach suff,$(MERGE_SUFFIXES),merged_vcf/%.$(suff).vcf)
	$(call LSCRIPT_MEM,6G,8G,"$(call GATK_MEM,5G) -T CombineVariants -R $(REF_FASTA) $(foreach suff,$(MERGE_SUFFIXES),--variant:$(suff) merged_vcf/$*.$(suff).vcf ) -o $@ -genotypeMergeOptions UNIQUIFY")

merged_tables/%.$(MERGE_SUFFIX).$(VCF_SUFFIX).txt : merged_vcf/%.$(MERGE_SUFFIX).$(VCF_SUFFIX).vcf
	$(INIT) S1=`grep '^#CHROM' $< | cut -f 10`; S2=`grep '^#CHROM' $< | cut -f 11`; \
		S3=`grep '^#CHROM' $< | cut -f 12`; \
		S4=`grep '^#CHROM' $< | cut -f 13`; \
	$(VCF_EFF_ONE_PER_LINE) < $< | $(call SNP_SIFT_MEM,2G) extractFields - $(ALL_VCF_EFF_FIELDS) set | sed "1s/GEN\[0\]/$$S1/g; 1s/GEN\[1\]/$$S2/g; 1s/GEN\[2\]/$$S3/g; 1s/GEN\[3\]/$$S4/g " | $(PERL) $(VCF_JOIN_EFF) > $@
	
merged_tables/%.ft.txt : merged_tables/%.txt
	grep -v 'FilteredInAll' $< > $@

merged_tables/all.%.txt : $(foreach tumor,$(TUMOR_SAMPLES),merged_tables/$(tumor).%.txt)
	$(INIT) $(RBIND) --sampleName $< $^ > $@

include ~/share/modules/vcf_tools/vcftools.mk
