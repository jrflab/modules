# Merge vcf files

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all 

LOGDIR = log/vcf_merge.$(NOW)
VCF_SUFFIX := dp_ft.som_ft.dbsnp.nsfp.chasm.fathmm.eff
TABLE_SUFFIX := $(VCF_SUFFIX).tab.ft
EFF_TYPES = silent missense nonsilent_cds nonsilent
TABLE_SUFFIXES := $(foreach eff,$(EFF_TYPES),$(TABLE_SUFFIX).$(eff))
VARIANT_TYPES ?= mutect museq gatk_snps
MERGE_SUFFIX = $(subst $( ),_,$(VARIANT_TYPES))
VCF_FIELDS += set


#all : $(foreach tumor,$(TUMOR_SAMPLES),merged_tables/$(tumor)_$(normal_lookup.$(tumor)).$(MERGE_SUFFIX).$(VCF_SUFFIX).txt)
all : $(foreach suff,$(TABLE_SUFFIXES),tables/all.$(MERGE_SUFFIX).$(suff).txt)

vcf/%.$(MERGE_SUFFIX).vcf : $(foreach type,$(VARIANT_TYPES),vcf/%.$(type).vcf vcf/%.$(type).vcf.idx )
	$(call LSCRIPT_MEM,6G,8G,"$(call GATK_MEM,5G) -T CombineVariants -R $(REF_FASTA) $(foreach type,$(VARIANT_TYPES),--variant:$(type) vcf/$*.$(type).vcf ) -o $@ -genotypeMergeOptions UNSORTED")

#tables/%.$(MERGE_SUFFIX).$(VCF_SUFFIX).txt : vcf/%.$(MERGE_SUFFIX).$(VCF_SUFFIX).vcf
#S1=`grep '^#CHROM' $< | cut -f 10`; S2=`grep '^#CHROM' $< | cut -f 11`; \
#$(VCF_EFF_ONE_PER_LINE) < $< | $(call SNP_SIFT_MEM,2G) extractFields - $(ALL_VCF_EFF_FIELDS) set | sed "1s/GEN\[0\]/$$S1/g; 1s/GEN\[1\]/$$S2/g" | $(PERL) $(VCF_JOIN_EFF) > $@
	
tables/%.ft.txt : tables/%.txt
	grep -v 'FilteredInAll' $< > $@

tables/all.$(MERGE_SUFFIX).$(TABLE_SUFFIX).txt : $(foreach tumor,$(TUMOR_SAMPLES),tables/$(tumor)_$(normal_lookup.$(tumor)).$(MERGE_SUFFIX).$(TABLE_SUFFIX).txt)
	$(INIT) $(RSCRIPT) $(RBIND) --tumorNormal $^ > $@


#$(call LSCRIPT_MEM,2G,5G,"S1=`grep '^#CHROM' $< | cut -f 10`; S2=`grep '^#CHROM' $< | cut -f 11`; \
#$(VCF_EFF_ONE_PER_LINE) < $< | $(call SNP_SIFT_MEM,2G) extractFields $< $(ALL_VCF_EFF_FIELDS) set | sed \"1s/GEN\[0\]/\$$S1/g; 1s/GEN\[1\]/\$$S2/g\" | $(PERL) $(VCF_JOIN_EFF) > $@")

include ~/share/modules/vcftools.mk

