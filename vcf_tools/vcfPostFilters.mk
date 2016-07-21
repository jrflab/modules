# post-annotation filters
VCF_POST_ANN_FILTER_EXPRESSION ?= ExAC_AF > 0.1
vcf/%.cft.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call LSCRIPT_CHECK_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ \
		--filterExpression '$(VCF_POST_ANN_FILTER_EXPRESSION)' --filterName customFilter"))

COMMON_FILTER_VCF = $(PYTHON) modules/vcf_tools/common_filter_vcf.py
vcf/%.common_ft.vcf : vcf/%.vcf
	$(call CHECK_VCF,$(call LSCRIPT_MEM,4G,5G,"$(COMMON_FILTER_VCF) $< > $@"))
