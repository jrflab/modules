# Merge vcf files

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all 

LOGDIR = log/vcf_merge.$(NOW)
VCF_SUFFIX := dp_ft.dbsnp.nsfp.eff
TABLE_SUFFIX := $(VCF_SUFFIX).ft
VARIANT_TYPES ?= mutect museq
MERGE_SUFFIX = $(subst $( ),_,$(VARIANT_TYPES))
VCF_SAMPLES = 0 1
ANN_VCF_FIELDS = ReadPosRankSum MQRankSum MQ MQ0 HRun FS OND BaseQRankSum HaplotypeScore

#all : $(foreach tumor,$(TUMOR_SAMPLES),merged_tables/$(tumor)_$(normal_lookup.$(tumor)).$(MERGE_SUFFIX).$(VCF_SUFFIX).txt)
all : merged_tables/all.$(MERGE_SUFFIX).$(TABLE_SUFFIX).txt

vcf/%.$(MERGE_SUFFIX).vcf : $(foreach type,$(VARIANT_TYPES),vcf/%.$(type).vcf)
	$(call LSCRIPT_MEM,6G,8G,"$(call GATK_MEM,5G) -T CombineVariants -R $(REF_FASTA) $(foreach type,$(VARIANT_TYPES),--variant:$(type) vcf/$*.$(type).vcf ) -o $@ -genotypeMergeOptions UNSORTED")

tables/%.$(MERGE_SUFFIX).$(VCF_SUFFIX).txt : vcf/%.$(MERGE_SUFFIX).$(VCF_SUFFIX).vcf
	S1=`grep '^#CHROM' $< | cut -f 10`; S2=`grep '^#CHROM' $< | cut -f 11`; \
	$(VCF_EFF_ONE_PER_LINE) < $< | $(call SNP_SIFT_MEM,2G) extractFields - $(ALL_VCF_EFF_FIELDS) set | sed "1s/GEN\[0\]/$$S1/g; 1s/GEN\[1\]/$$S2/g" | $(PERL) $(VCF_JOIN_EFF) > $@
	
tables/%.ft.txt : tables/%.txt
	grep -v 'FilteredInAll' $< > $@

tables/all.$(MERGE_SUFFIX).$(TABLE_SUFFIX).txt : $(foreach tumor,$(TUMOR_SAMPLES),tables/$(tumor)_$(normal_lookup.$(tumor)).$(MERGE_SUFFIX).$(TABLE_SUFFIX).txt)
	$(INIT) $(RSCRIPT) $(RBIND) --tumorNormal $^ > $@

VCF_ANNOTATIONS = RMSMappingQuality AlleleBalance ClippingRankSumTest QualByDepth TandemRepeatAnnotator ReadPosRankSumTest MappingQualityZero MappingQualityRankSumTest HaplotypeScore FisherStrand Coverage BaseQualityRankSumTest HomopolymerRun
define variant-annotate-tumor-normal
merged_vcf/$1_$2.%.ann.vcf : merged_vcf/$1_$2.%.vcf bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_PARALLEL_MEM,4,2G,3G,"$$(call GATK_MEM,8G) -T VariantAnnotator -nt 4 -R $$(REF_FASTA) $$(foreach ann,$$(VCF_ANNOTATIONS),-A $$(ann) ) --dbsnp $$(DBSNP) $$(foreach i,$$(filter %.bam,$$^),-I $$i) -V $$< -o $$@")
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call variant-annotate-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))


#$(call LSCRIPT_MEM,2G,5G,"S1=`grep '^#CHROM' $< | cut -f 10`; S2=`grep '^#CHROM' $< | cut -f 11`; \
#$(VCF_EFF_ONE_PER_LINE) < $< | $(call SNP_SIFT_MEM,2G) extractFields $< $(ALL_VCF_EFF_FIELDS) set | sed \"1s/GEN\[0\]/\$$S1/g; 1s/GEN\[1\]/\$$S2/g\" | $(PERL) $(VCF_JOIN_EFF) > $@")

include ~/share/modules/vcftools.mk
