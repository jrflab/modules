DEPTH_FILTER ?= 5

FALSE_POSITIVE_BED = $(HOME)/share/reference/fuentes_blacklist.include_cosmic.hg19.bed
vcf/%.fp_ft.vcf : vcf/%.vcf
	$(call RUN,-c -s 8G -m 12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@.tmp --maskName 'FuentesFalsePositive' --mask $(FALSE_POSITIVE_BED) && $(call VERIFY_VCF,$@.tmp,$@)")

DGD_BED = $(HOME)/share/reference/dgd.include_cosmic.hg19.bed
vcf/%.dgd_ft.vcf : vcf/%.vcf
	$(call RUN,-c -s 8G -m 12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@.tmp --maskName 'DuplicateGenesDB' --mask $(DGD_BED) && $(call VERIFY_VCF,$@.tmp,$@)")

ENCODE_BED = $(HOME)/share/reference/wgEncodeDacMapabilityConsensusExcludable.include_cosmic.bed
vcf/%.encode_ft.vcf : vcf/%.vcf
	$(call RUN,-c -s 8G -m 12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@.tmp --maskName 'encode' --mask $(ENCODE_BED) && $(call VERIFY_VCF,$@.tmp,$@)")

%.het_ft.vcf : %.vcf
	$(call RUN,-c -s 8G -m 12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@.tmp \
		--genotypeFilterExpression 'isHet == 1' --genotypeFilterName 'Heterozygous positions' && $(call VERIFY_VCF,$@.tmp,$@)")

vcf/%.dbsnp_ft.vcf : vcf/%.vcf
	$(INIT) awk '/^#/ || $$3 ~ /^rs/ {print}' $< > $@

# apply overall depth filter
vcf/%.dp_ft.vcf : vcf/%.vcf
	$(call RUN,-c -s 8G -m 12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@.tmp \
		--filterExpression 'DP < $(DEPTH_FILTER)' --filterName Depth && $(call VERIFY_VCF,$@.tmp,$@)")

vcf/%.sdp_ft.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call RUN,-c -s 2G -m 5G,"$(call SNP_SIFT_MEM,2G) filter $(SNP_SIFT_OPTS) \
		-f $< '(exists GEN[*].DP) & (GEN[*].DP > 20)' > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)"))

# apply HRun filter
vcf/%.hrun_ft.vcf : vcf/%.vcf
	$(call RUN,-c -s 8G -m 12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@.tmp \
		--filterExpression 'HRun > $(HRUN_FILTER)' --filterName HRun && $(call VERIFY_VCF,$@.tmp,$@)")

vcf/%.vaf_ft.vcf : vcf/%.vcf
	$(call RUN,-c -s 2G -m 5G,"$(call SNP_SIFT_MEM,2G) filter $(SNP_SIFT_OPTS) '(exists GEN[*].VAF) & (GEN[0].VAF > 0.05) & (GEN[1].VAF < 0.05)' < $< > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)")

# target region filter
INTERVAL_FILTER_VCF = python modules/vcf_tools/interval_filter_vcf.py
%.target_ft.vcf : %.vcf
	$(call RUN,-c -s 4G -m 6G,"$(INTERVAL_FILTER_VCF) $(TARGETS_FILE) $< > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)")

INTERVAL_DEPTH_FILTER_VCF = python modules/vcf_tools/interval_depth_filter_vcf.py
INTERVAL_DEPTH_THRESHOLD ?= 50
%.target_dp_ft.vcf : %.vcf
	$(call RUN,-c -s 4G -m 6G,"$(INTERVAL_DEPTH_FILTER_VCF) --depth_threshold $(INTERVAL_DEPTH_THRESHOLD) \
		$(TARGETS_FILE) $< > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)")

# varscan depth filter (b/c varscan is dumb and only gives variant depth)
vcf/%.vdp_ft.vcf : vcf/%.vcf
	$(call RUN,-c -s 2G -m 5G,"cat $< | $(call SNP_SIFT_MEM,2G) filter $(SNP_SIFT_OPTS) '(exists GEN[*].AD) & (GEN[*].AD > $(DEPTH_FILTER))' > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)")

ifdef SAMPLE_SET_PAIRS
define somatic-filter-vcf-set
vcf/$1.%.som_ft.vcf : vcf/$1.%.vcf
	$$(INIT) $$(SOMATIC_FILTER_VCF) -n $(normal.$1) -f 0.03 $$< > $$@ 2> $$(LOG)
endef
$(foreach set,$(SAMPLE_SET_PAIRS),$(eval $(call somatic-filter-vcf-set,$(set))))
endif

MUTECT2_SOMATIC_AD_FILTER_VCF = python modules/vcf_tools/somatic_ad_filter_vcf.py
ifdef SAMPLE_PAIRS
# ff normal filter :
# filter if normal depth > 20 and normal VAF > 1/5 * tumor VAF
# or normal variant depth greater than 1
# ffpe normal filter :
# filter if normal depth > 20 and normal variant depth > 1/3 * tumor variant depth
# or normal variant depth greater than 1
define som-ad-ft-tumor-normal
vcf/$1_$2.%.som_ad_ft.vcf : vcf/$1_$2.%.vcf
	$$(call RUN,-c -s 4G -m 7G,"$$(MUTECT2_SOMATIC_AD_FILTER_VCF) --tumor $1 --normal $2 $$< > $$@.tmp && $$(call VERIFY_VCF,$$@.tmp,$$@)")

vcf/$1_$2.%.ffpe_som_ad_ft.vcf : vcf/$1_$2.%.vcf
	$$(call RUN,-c -s 8G -m 12G,"$$(call GATK_MEM,8G) -T VariantFiltration -R $$(REF_FASTA) -V $$< -o $$@ \
		--filterExpression 'if (vc.getGenotype(\"$2\").getDP() > 20) { ( vc.getGenotype(\"$2\").getAD().1 * 1.0 / vc.getGenotype(\"$2\").getDP()) > ( vc.getGenotype(\"$1\").getAD().1 * 1.0 / vc.getGenotype(\"$1\").getDP()) / 3.0 } else { vc.getGenotype(\"$2\").getAD().1 > 1 }' \
		--filterName somaticAlleleDepth \
		--filterExpression 'vc.getGenotype(\"$1\").getDP() <= $$(DEPTH_FILTER) || vc.getGenotype(\"$2\").getDP() <= $$(DEPTH_FILTER)' \
		--filterName depthFilter && sed -i 's/getGenotype(\"\([^\"]*\)\")/getGenotype(\1)/g' $$@")

# somatic filter for structural variants
vcf/$1_$2.%.sv_som_ft.vcf : vcf/$1_$2.%.vcf
	$$(call RUN,-c -s 8G -m 12G,"$$(call GATK_MEM,8G) -T VariantFiltration -R $$(REF_FASTA) -V $$< -o $$@ \
		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"SU\") <= $$(DEPTH_FILTER)' \
		--filterName svSupport \
		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"SU\") < 5 * vc.getGenotype(\"$2\").getAnyAttribute(\"SU\")' \
		--filterName somaticSvSupport && sed -i 's/getGenotype(\"\([^\"]*\)\")/getGenotype(\1)/g' $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call som-ad-ft-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif

NORMAL_VCF ?= $(HOME)/share/reference/spowellnormal.gatk_variants.vcf
ifdef NORMAL_VCF
vcf/%.nft.vcf : vcf/%.vcf
	$(call RUN,-c -s 8G -m 12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@.tmp --maskName 'normal' --mask $(NORMAL_VCF) && $(call VERIFY_VCF,$@.tmp,$@)")

# allele count filtering for hotspots: any alt allele count > 0
vcf/%.ac_ft.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call RUN,-c -s 12G -m 14G,"$(call SNP_SIFT_MEM,8G) filter \
		$(SNP_SIFT_OPTS) \" ( AC[*] > 0 ) \" $< > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)"))
endif
