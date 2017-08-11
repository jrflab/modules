# post-annotation filters
VCF_POST_ANN_FILTER_EXPRESSION ?= ExAC_AF > 0.1
vcf/%.cft.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call RUN,-c -s 8G -m 12G,"$(call GATK_MEM,8G) -T VariantFiltration -U LENIENT_VCF_PROCESSING -R $(REF_FASTA) -V $< -o $@.tmp \
		--filterExpression '$(VCF_POST_ANN_FILTER_EXPRESSION)' --filterName customFilter && $(call VERIFY_VCF,$@.tmp,$@)"))

COMMON_FILTER_VCF = $(PYTHON) modules/vcf_tools/common_filter_vcf.py
vcf/%.common_ft.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call RUN,-c -s 4G -m 5G,"$(COMMON_FILTER_VCF) $< > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)"))
