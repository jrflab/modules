# dbsnp annotations
vcf/%.dbsnp.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call LSCRIPT_CHECK_MEM,20G,25G,"$(call SNP_SIFT_MEM,10G) annotate \
		$(SNP_SIFT_OPTS) $(DBSNP) $< > $@.tmp && if grep -q '^#CHROM' $@.tmp; then mv $@.tmp $@; fi"))

vcf/%.hotspot_ann.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call LSCRIPT_CHECK_MEM,20G,23G,"$(call SNP_SIFT_MEM,10G) annotate $(SNP_SIFT_OPTS) \
		$(HOTSPOT_UNMERGED_VCF) $< > $@.tmp && if grep -q '^#CHROM' $@.tmp; then mv $@.tmp $@; fi"))

# mouse genome project dbsnp
vcf/%.mgp_dbsnp.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call LSCRIPT_CHECK_MEM,13G,15G,"$(call SNP_SIFT_MEM,10G) annotate \
		-tabix $(SNP_SIFT_OPTS) $(MGP_SNP_DBSNP) $< | $(call SNP_SIFT_MEM,10G) annotate \
		-tabix $(SNP_SIFT_OPTS) $(MGP_INDEL_DBSNP) > $@.tmp && if grep -q '^#CHROM' $@.tmp; then mv $@.tmp $@; fi"))

vcf/%.cosmic.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call LSCRIPT_CHECK_MEM,20G,24G,"$(call SNP_SIFT_MEM,10G) annotate $(SNP_SIFT_OPTS) \
		$(COSMIC) $< > $@.tmp && if grep -q '^#CHROM' $@.tmp; then mv $@.tmp $@; fi"))

TRANSFIC = $(RSCRIPT) modules/vcf_tools/transficVcf.R
TRANSFIC_PERL_SCRIPT = $(HOME)/share/usr/transfic/bin/transf_scores.pl
vcf/%.transfic.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call LSCRIPT_MEM,9G,12G,"$(TRANSFIC) --genome $(REF) --transfic $(TRANSFIC_PERL_SCRIPT) --outFile $@ $<")

# add exon distance
vcf/%.exondist.vcf : vcf/%.vcf
	$(call LSCRIPT_CHECK_MEM,2G,3G,"$(INTRON_POSN_LOOKUP) $< > $@")

# run snp eff
SNP_EFF_FLAGS ?= -canon # -ud 0  -no-intron -no-intergenic -no-utr
SNP_EFF_OPTS = -c $(SNP_EFF_CONFIG) -i vcf -o vcf $(SNP_EFF_FLAGS)
vcf/%.eff.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call LSCRIPT_CHECK_MEM,18G,20G,"$(call SNP_EFF_MEM,11G) ann $(SNP_EFF_OPTS) $(SNP_EFF_GENOME) -s $(@D)/$*.eff_summary.html $< > $@.tmp && if grep -q '^#CHROM' $@.tmp; then mv $@.tmp $@; fi"))


vcf/%.clinvar.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call LSCRIPT_CHECK_MEM,18G,23G,"$(call SNP_SIFT_MEM,11G) annotate $(SNP_SIFT_OPTS) \
		$(CLINVAR) $< > $@.tmp && if grep -q '^#CHROM' $@.tmp; then mv $@.tmp $@; fi"))

vcf/%.exac_nontcga.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call LSCRIPT_CHECK_MEM,18G,20G,"$(call SNP_SIFT_MEM,11G) annotate $(SNP_SIFT_OPTS) \
		-info ExAC_AF $(EXAC_NONTCGA) $< > $@.tmp && if grep -q '^#CHROM' $@.tmp; then mv $@.tmp $@; fi"))

HAPLOTYPE_INSUF_BED = $(HOME)/share/reference/haplo_insuff_genes.bed
CANCER_GENE_CENSUS_BED = $(HOME)/share/reference/annotation_gene_lists/cancer_gene_census_genes_v20150303.bed
KANDOTH_BED = $(HOME)/share/reference/annotation_gene_lists/Kandoth_127genes.bed
LAWRENCE_BED = $(HOME)/share/reference/annotation_gene_lists/Lawrence_cancer5000-S.bed

ADD_GENE_LIST_ANNOTATION = $(RSCRIPT) modules/vcf_tools/addGeneListAnnotationToVcf.R
vcf/%.gene_ann.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call LSCRIPT_CHECK_MEM,8G,12G,"$(ADD_GENE_LIST_ANNOTATION) --genome $(REF) --geneBed $(HAPLOTYPE_INSUF_BED)$(,)$(CANCER_GENE_CENSUS_BED)$(,)$(KANDOTH_BED)$(,)$(LAWRENCE_BED) --name hap_insuf$(,)cancer_gene_census$(,)kandoth$(,)lawrence --outFile $@ $<"))

# Copy number regulated genes annotated per subtype
# FYI Endometrioid_MSI-L has no copy number regulated genes
CN_ENDOMETRIAL_SUBTYPES = CN_high CN_low Endometrioid_MSI_H Endometrioid_MSS Endometrioid MSI POLE Serous
CN_BREAST_SUBTYPES = ER_negative ER_positive HER2_postitive Pam50_Basal Pam50_Her2 Pam50_LumA Pam50_LumB Pam50_Normal Triple_negative
CN_ENDOMETRIAL_BED = $(foreach set,$(CN_ENDOMETRIAL_SUBTYPES), $(HOME)/share/reference/annotation_gene_lists/cn_reg/endometrial/copy_number_regulated_genes_subtype_$(set)_spearmanrsquare0.4_fdrbh_adjp_lt0.05.HUGO.bed)
CN_BREAST_BED = $(foreach set,$(CN_BREAST_SUBTYPES), $(HOME)/share/reference/annotation_gene_lists/cn_reg/breast/metabric_subtype_$(set)_copy_number_regulated_genes_std0.5_spearmanrsquare0.4_fdrbh_adjp_lt0.05.HUGO.bed)
vcf/%.cn_reg.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call LSCRIPT_MEM,8G,12G,"$(ADD_GENE_LIST_ANNOTATION) --genome $(REF) --geneBed $(subst $(space),$(,),$(strip $(CN_ENDOMETRIAL_BED)) $(strip $(CN_BREAST_BED))) --name $(subst $(space),$(,),$(foreach set,$(strip $(CN_ENDOMETRIAL_SUBTYPES)),endometrial_$(set)) $(foreach set,$(strip $(CN_BREAST_SUBTYPES)),breast_$(set))) --outFile $@ $<"))


define ad-tumor-normal
vcf/$1_$2.%.ad.vcf : vcf/$1_$2.%.vcf bam/$1.bam bam/$2.bam bam/$1.bai bam/$2.bai
	$$(call LSCRIPT_CHECK_PARALLEL_MEM,4,2G,3G,"$$(call GATK_MEM,8G) -T VariantAnnotator -nt 4 -R $$(REF_FASTA) -A DepthPerAlleleBySample --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -o $$@ -L $$<")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call ad-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

ANNOVAR_PROTOCOL ?= refGene$(,)cytoBand$(,)genomicSuperDups$(,)esp6500siv2_all$(,)1000g2014oct_all$(,)1000g2014oct_afr$(,)1000g2014oct_eas$(,)1000g2014oct_eur$(,)snp138$(,)ljb26_all
ANNOVAR_OPERATION ?= g$(,)r$(,)r$(,)f$(,)f$(,)f$(,)f$(,)f$(,)f$(,)f
ANNOVAR_OPTS = --dot2underline -remove -protocol $(ANNOVAR_PROTOCOL) -operation $(ANNOVAR_OPERATION) -nastring . -vcfinput -buildver $(ANNOVAR_REF)
vcf/%.$(ANNOVAR_REF)_multianno.vcf : vcf/%.vcf
	$(call LSCRIPT_CHECK_MEM,7G,9G,"$(ANNOVAR) -out $(@D)/$* $(ANNOVAR_OPTS) $< $(ANNOVAR_DB)")
