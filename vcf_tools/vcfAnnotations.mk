# dbsnp annotations
vcf/%.dbsnp.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call RUN,-c -s 20G -m 25G,"$(call SNP_SIFT_MEM,10G) annotate \
		$(SNP_SIFT_OPTS) $(DBSNP) $< > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)"))

vcf/%.hotspot_int_ann.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call RUN,-c -s 20G -m 23G,"$(call SNP_SIFT_MEM,10G) annotate $(SNP_SIFT_OPTS) \
		$(HOTSPOT_VCF.int) $< > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)"))

vcf/%.hotspot_ext_ann.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call RUN,-c -s 20G -m 23G,"$(call SNP_SIFT_MEM,10G) annotate $(SNP_SIFT_OPTS) \
		$(HOTSPOT_VCF.ext) $< > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)"))

# mouse genome project dbsnp
vcf/%.mgp_dbsnp.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call RUN,-c -s 13G -m 15G,"$(call SNP_SIFT_MEM,10G) annotate \
		-tabix $(SNP_SIFT_OPTS) $(MGP_SNP_DBSNP) $< | $(call SNP_SIFT_MEM,10G) annotate \
		-tabix $(SNP_SIFT_OPTS) $(MGP_INDEL_DBSNP) > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)"))

vcf/%.cosmic.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call RUN,-c -s 20G -m 24G,"$(call SNP_SIFT_MEM,10G) annotate $(SNP_SIFT_OPTS) \
		$(COSMIC) $< > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)"))

vcf/%.cosmic_nc.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call RUN,-c -s 20G -m 24G,"$(call SNP_SIFT_MEM,10G) annotate $(SNP_SIFT_OPTS) \
		$(COSMIC_NONCODING) $< > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)"))

TRANSFIC = $(RSCRIPT) modules/vcf_tools/transficVcf.R
TRANSFIC_PERL_SCRIPT = $(HOME)/share/usr/transfic/bin/transf_scores.pl
vcf/%.transfic.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call RUN,-s 9G -m 12G,"$(TRANSFIC) --genome $(REF) --transfic $(TRANSFIC_PERL_SCRIPT) --outFile $@.tmp $< && $(call VERIFY_VCF,$@.tmp,$@)")

# add exon distance
vcf/%.exondist.vcf : vcf/%.vcf
	$(call RUN,-c -s 2G -m 3G,"$(INTRON_POSN_LOOKUP) $< > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)")

# run snp eff
SNP_EFF_FLAGS ?= -canon # -ud 0  -no-intron -no-intergenic -no-utr
SNP_EFF_OPTS = -c $(SNP_EFF_CONFIG) -i vcf -o vcf $(SNP_EFF_FLAGS)
vcf/%.eff.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call RUN,-c -s 18G -m 20G,"$(call SNP_EFF_MEM,11G) ann $(SNP_EFF_OPTS) $(SNP_EFF_GENOME) -s $(@D)/$*.eff_summary.html $< > $@.tmp \
		&& $(call VERIFY_VCF,$@.tmp,$@)"))

vcf/%.clinvar.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call RUN,-c -s 18G -m 23G,"$(call SNP_SIFT_MEM,11G) annotate $(SNP_SIFT_OPTS) \
		$(CLINVAR) $< > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)"))

vcf/%.exac_nontcga.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call RUN,-c -s 18G -m 20G,"$(call SNP_SIFT_MEM,11G) annotate $(SNP_SIFT_OPTS) \
		-info ExAC_AF $(EXAC_NONTCGA) $< > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)"))

GNOMAD_INFO ?= AF
vcf/%.gnomad.vcf : $(foreach chr,$(GNOMAD_CHROMOSOMES),chr_vcf/%.gnomad.$(chr).vcf)
	$(call RUN, -c -s 6G -m 7G,"$(MERGE_VCF) --out_file $@.tmp $^ && $(call VERIFY_VCF,$@.tmp,$@)")

define gnomad-chr
chr_vcf/%.gnomad.$1.vcf : vcf/%.vcf
	$$(call CHECK_VCF,$$(call RUN,-c -s 18G -m 20G,"$$(call SNP_SIFT_MEM,11G) annotate $$(SNP_SIFT_OPTS) \
		-info $$(GNOMAD_INFO) $$(GNOMAD_DB_DIR)/$$(GNOMAD_PREFIX).$1.vcf.gz $$< > $$@.tmp && $$(call VERIFY_VCF,$$@.tmp,$$@)"))
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call gnomad-chr,$(chr))))
 
ADD_GENE_LIST_ANNOTATION = $(RSCRIPT) modules/vcf_tools/addGeneListAnnotationToVcf.R
vcf/%.gene_ann.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call RUN,-c -s 8G -m 12G,"$(ADD_GENE_LIST_ANNOTATION) --genome $(REF) \
		--geneBed $(HAPLOTYPE_INSUF_BED)$(,)$(CANCER_GENE_CENSUS_BED)$(,)$(KANDOTH_BED)$(,)$(LAWRENCE_BED) \
		--name hap_insuf$(,)cancer_gene_census$(,)kandoth$(,)lawrence --outFile $@.tmp $< && \
		$(call VERIFY_VCF,$@.tmp,$@)"))

# Copy number regulated genes annotated per subtype
# FYI Endometrioid_MSI-L has no copy number regulated genes
CN_ENDOMETRIAL_SUBTYPES = CN_high CN_low Endometrioid_MSI_H Endometrioid_MSS Endometrioid MSI POLE Serous
CN_BREAST_SUBTYPES = ER_negative ER_positive HER2_postitive Pam50_Basal Pam50_Her2 Pam50_LumA Pam50_LumB Pam50_Normal Triple_negative
CN_ENDOMETRIAL_BED = $(foreach set,$(CN_ENDOMETRIAL_SUBTYPES), $(HOME)/share/reference/annotation_gene_lists/cn_reg/endometrial/copy_number_regulated_genes_subtype_$(set)_spearmanrsquare0.4_fdrbh_adjp_lt0.05.HUGO.bed)
CN_BREAST_BED = $(foreach set,$(CN_BREAST_SUBTYPES), $(HOME)/share/reference/annotation_gene_lists/cn_reg/breast/metabric_subtype_$(set)_copy_number_regulated_genes_std0.5_spearmanrsquare0.4_fdrbh_adjp_lt0.05.HUGO.bed)
vcf/%.cn_reg.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call RUN,-c -s 8G -m 12G,"$(ADD_GENE_LIST_ANNOTATION) --genome $(REF) \
		--geneBed $(subst $(space),$(,),$(strip $(CN_ENDOMETRIAL_BED)) $(strip $(CN_BREAST_BED))) \
		--name $(subst $(space),$(,),$(foreach set,$(strip $(CN_ENDOMETRIAL_SUBTYPES)),endometrial_$(set)) \
		$(foreach set,$(strip $(CN_BREAST_SUBTYPES)),breast_$(set))) --outFile $@.tmp $< && \
		$(call VERIFY_VCF,$@.tmp,$@)"))


define ad-tumor-normal
vcf/$1_$2.%.ad.vcf : vcf/$1_$2.%.vcf bam/$1.bam bam/$2.bam bam/$1.bai bam/$2.bai
	$$(call RUN,-c -n 4 -s 2G -m 3G,"$$(call GATK_MEM,8G) -T VariantAnnotator -nt 4 -R $$(REF_FASTA) \
		-A DepthPerAlleleBySample --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -o $$@ -L $$<")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(eval $(call ad-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

ANNOVAR_PROTOCOL ?= refGene$(,)cytoBand$(,)genomicSuperDups$(,)esp6500siv2_all$(,)1000g2014oct_all$(,)1000g2014oct_afr$(,)1000g2014oct_eas$(,)1000g2014oct_eur$(,)snp138$(,)ljb26_all
ANNOVAR_OPERATION ?= g$(,)r$(,)r$(,)f$(,)f$(,)f$(,)f$(,)f$(,)f$(,)f
ANNOVAR_OPTS = --dot2underline -remove -protocol $(ANNOVAR_PROTOCOL) -operation $(ANNOVAR_OPERATION) -nastring . -vcfinput -buildver $(ANNOVAR_REF)
vcf/%.$(ANNOVAR_REF)_multianno.vcf : vcf/%.vcf
	$(call RUN,-c -s 7G -m 9G,"$(ANNOVAR) -out $(@D)/$* $(ANNOVAR_OPTS) $< $(ANNOVAR_DB)")

ONCOTATOR = oncotator
ONCOTATOR_OPTS = -v --db-dir $(ONCOTATOR_DB) $(if $(ONCOTATOR_TX_OVERRIDES),-c $(ONCOTATOR_TX_OVERRIDES))
vcf/%.oncotator.vcf : vcf/%.vcf
	$(call RUN,-c -v $(ONCOTATOR_ENV) -s 8G -m 12G,"$(ONCOTATOR) $(ONCOTATOR_OPTS) -i VCF -o VCF $< $@.tmp $(ONCOTATOR_REF) && \
		sed -i 's/^##INFO=<ID=HapScore$(,)Number=.$(,)Type=Integer/##INFO=<ID=HapScore$(,)Number=.$(,)Type=String/' $@.tmp && \
		perl -lane 'if (/^#/) { print; } else { for \$$i (7 .. \$$#F) { \$$F[\$$i] =~ s/\|/$(,)/g; } print join \"\t\"$(,) @F;}' $@.tmp > $@ && \
		rm $@.tmp")

CMO_ANN = python modules/vcf_tools/annotate_vcf2maf.py
CMO_ANN_OPTS = --vcf2maf '$(VCF2MAF)' \
			   --vcf2maf_opts '--vep-path $(VEP_PATH) --vep-data $(VEP_DATA) --ncbi-build $(VEP_REF) \
			   --maf-center mskcc.org --tmp-dir $(shell mktemp -d) \
			   --custom-enst $(VEP_OVERRIDES) --species $(VEP_SPECIES)' \
			   --filter_vcf $(EXAC_NONTCGA) --ref_fasta $(REF_FASTA) \
			   $(if $(CMO_HOTSPOT_FILE), --hotspot_list $(CMO_HOTSPOT_FILE))

vcf/%.cmo_ann.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call RUN,-c -v $(VEP_ENV) -n 4 -s 3G -m 3G,"$(CMO_ANN) $(CMO_ANN_OPTS) \
		--vep_forks 4 $< > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)"))

VEP = variant_effect_predictor.pl
VEP_OPTS = --species $(VEP_SPECIES) --everything --cache --no_progress --format vcf --dir $(VEP_DATA) --force_overwrite --vcf --offline --assembly $(VEP_REF)
vcf/%.vep.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call RUN,-c -v $(VEP_ENV) -s 6G -m 8G,"$(VEP) $(VEP_OPTS) -i $< -o $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)"))

BED_ANNOTATE_VCF = python modules/vcf_tools/bed_annotate_vcf.py
vcf/%.fuentes.vcf : vcf/%.vcf
	$(call RUN,-c -s 4G -m 6G,"$(BED_ANNOTATE_VCF) --info_tag fuentes $(FUENTES_BED) $< > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)")

vcf/%.dgd.vcf : vcf/%.vcf
	$(call RUN,-c -s 4G -m 6G,"$(BED_ANNOTATE_VCF) --info_tag dgd $(DGD_BED) $< > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)")

ONCOKB_VCF = python modules/vcf_tools/oncokb_vcf.py
vcf/%.oncokb.vcf : vcf/%.vcf
	$(call RUN,-c -s 4G -m 6G,"$(ONCOKB_VCF) --oncokb $(ONCOKB) $< > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)")

