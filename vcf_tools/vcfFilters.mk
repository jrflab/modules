DEPTH_FILTER ?= 5

FALSE_POSITIVE_BED = $(HOME)/share/reference/fuentes_blacklist.include_cosmic.hg19.bed
vcf/ft/fp_ft/%.fp_ft.vcf : vcf/%.vcf
	$(call LSCRIPT_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ --maskName 'FuentesFalsePositive' --mask $(FALSE_POSITIVE_BED)")

DGD_BED = $(HOME)/share/reference/dgd.include_cosmic.hg19.bed
vcf/ft/dgd_ft/%.dgd_ft.vcf : vcf/dgd_ft/%.vcf
	$(call LSCRIPT_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ --maskName 'DuplicateGenesDB' --mask $(DGD_BED)")

ENCODE_BED = $(HOME)/share/reference/wgEncodeDacMapabilityConsensusExcludable.include_cosmic.bed
vcf/ft/encode_ft/%.encode_ft.vcf : vcf/%.vcf
	$(call LSCRIPT_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ --maskName 'encode' --mask $(ENCODE_BED)")

%.het_ft.vcf : %.vcf
	$(call LSCRIPT_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ \
		--genotypeFilterExpression 'isHet == 1' --genotypeFilterName 'Heterozygous positions'")

vcf/ft/dbsnp_ft/%.dbsnp_ft.vcf : vcf/%.vcf
	$(INIT) awk '/^#/ || $$3 ~ /^rs/ {print}' $< > $@

# apply overall depth filter
vcf/ft/dp_ft/%.dp_ft.vcf : vcf/%.vcf
	$(call LSCRIPT_CHECK_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ \
		--filterExpression 'DP < $(DEPTH_FILTER)' --filterName Depth")

%.sdp_ft.vcf : %.vcf
	$(call CHECK_VCF,$(call LSCRIPT_CHECK_MEM,2G,5G,"$(call SNP_SIFT_MEM,2G) filter $(SNP_SIFT_OPTS) \
		-f $< '(exists GEN[*].DP) & (GEN[*].DP > 20)' > $@"))

# apply HRun filter
vcf/ft/hrun_ft/%.hrun_ft.vcf : vcf/%.vcf
	$(call LSCRIPT_CHECK_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ \
		--filterExpression 'HRun > $(HRUN_FILTER)' --filterName HRun")

vcf/vaf_ft/%.vaf_ft.vcf : vcf/%.vcf
	$(call LSCRIPT_CHECK_MEM,2G,5G,"$(call SNP_SIFT_MEM,2G) filter $(SNP_SIFT_OPTS) '(exists GEN[*].VAF) & (GEN[0].VAF > 0.05) & (GEN[1].VAF < 0.05)' < $< > $@")

# target region filter
vcf/ft/target_ft/%.target_ft.vcf : vcf/%.vcf
	$(call LSCRIPT_CHECK_MEM,8G,12G,"$(call GATK_MEM2,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ \
		--mask $(TARGETS_FILE) --maskName targetInterval --filterNotInMask")

# varscan depth filter (b/c varscan is dumb and only gives variant depth)
vcf/ft/vdp_ft/%.vdp_ft.vcf : vcf/%.vcf
	$(call LSCRIPT_CHECK_MEM,2G,5G,"cat $< | $(call SNP_SIFT_MEM,2G) filter $(SNP_SIFT_OPTS) '(exists GEN[*].AD) & (GEN[*].AD > $(DEPTH_FILTER))' > $@")

ifdef SAMPLE_SET_PAIRS
define somatic-filter-vcf-set
vcf/som_ft/$1.%.som_ft.vcf : vcf/$1.%.vcf
	$$(INIT) $$(SOMATIC_FILTER_VCF) -n $(normal.$1) -f 0.03 $$< > $$@ 2> $$(LOG)
endef
$(foreach set,$(SAMPLE_SET_PAIRS),$(eval $(call somatic-filter-vcf-set,$(set))))
endif

ifdef SAMPLE_PAIRS
# ff normal filter :
# filter if normal depth > 20 and normal VAF > 1/5 * tumor VAF
# or normal variant depth greater than 1
# ffpe normal filter :
# filter if normal depth > 20 and normal variant depth > 1/3 * tumor variant depth
# or normal variant depth greater than 1
define som-ad-ft-tumor-normal
vcf/ft/som_ad_ft/$1_$2.%.som_ad_ft.vcf : vcf/$1_$2.%.vcf
	$$(call LSCRIPT_CHECK_MEM,8G,12G,"$$(call GATK_MEM,8G) -T VariantFiltration -R $$(REF_FASTA) -V $$< -o $$@ \
		--filterExpression 'vc.getGenotype(\"$1\").getAD().1 < $(DEPTH_FILTER)' \
		--filterName tumorVarAlleleDepth \
		--filterExpression 'if (vc.getGenotype(\"$2\").getDP() > 20) { ( vc.getGenotype(\"$2\").getAD().1 * 1.0 / vc.getGenotype(\"$2\").getDP()) > ( vc.getGenotype(\"$1\").getAD().1 * 1.0 / vc.getGenotype(\"$1\").getDP()) / 5.0 } else { vc.getGenotype(\"$2\").getAD().1 > 1 }' \
		--filterName somaticAlleleDepth \
		--filterExpression 'vc.getGenotype(\"$1\").getDP() <= $$(DEPTH_FILTER) || vc.getGenotype(\"$2\").getDP() <= $$(DEPTH_FILTER)' \
		--filterName depthFilter && sed -i 's/getGenotype(\"\([^\"]*\)\")/getGenotype(\1)/g' $$@")

vcf/ft/ffpe_som_ad_ft/$1_$2.%.ffpe_som_ad_ft.vcf : vcf/$1_$2.%.vcf
	$$(call LSCRIPT_CHECK_MEM,8G,12G,"$$(call GATK_MEM,8G) -T VariantFiltration -R $$(REF_FASTA) -V $$< -o $$@ \
		--filterExpression 'vc.getGenotype(\"$1\").getAD().1 < $(DEPTH_FILTER)' \
		--filterName tumorVarAlleleDepth \
		--filterExpression 'if (vc.getGenotype(\"$2\").getDP() > 20) { ( vc.getGenotype(\"$2\").getAD().1 * 1.0 / vc.getGenotype(\"$2\").getDP()) > ( vc.getGenotype(\"$1\").getAD().1 * 1.0 / vc.getGenotype(\"$1\").getDP()) / 3.0 } else { vc.getGenotype(\"$2\").getAD().1 > 1 }' \
		--filterName somaticAlleleDepth \
		--filterExpression 'vc.getGenotype(\"$1\").getDP() <= $$(DEPTH_FILTER) || vc.getGenotype(\"$2\").getDP() <= $$(DEPTH_FILTER)' \
		--filterName depthFilter && sed -i 's/getGenotype(\"\([^\"]*\)\")/getGenotype(\1)/g' $$@")

# somatic filter for structural variants
vcf/ft/sv_som_ft/$1_$2.%.sv_som_ft.vcf : vcf/$1_$2.%.vcf
	$$(call LSCRIPT_CHECK_MEM,8G,12G,"$$(call GATK_MEM,8G) -T VariantFiltration -R $$(REF_FASTA) -V $$< -o $$@ \
		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"SU\") <= $$(DEPTH_FILTER)' \
		--filterName svSupport \
		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"SU\") < 5 * vc.getGenotype(\"$2\").getAnyAttribute(\"SU\")' \
		--filterName somaticSvSupport && sed -i 's/getGenotype(\"\([^\"]*\)\")/getGenotype(\1)/g' $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call som-ad-ft-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif

NORMAL_VCF ?= $(HOME)/share/reference/spowellnormal.gatk_variants.vcf
ifdef NORMAL_VCF
vcf/nft/%.nft.vcf : vcf/%.vcf
	$(call LSCRIPT_CHECK_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ --maskName 'normal' --mask $(NORMAL_VCF)")

# allele count filtering for hotspots: any alt allele count > 0
vcf/ft/ac_ft/%.ac_ft.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call LSCRIPT_CHECK_MEM,9G,12G,"$(call SNP_SIFT_MEM,8G) filter \
		$(SNP_SIFT_OPTS) \" ( AC[*] > 0 ) \" $< > $@"))
endif
