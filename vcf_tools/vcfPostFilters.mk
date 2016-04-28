# post-annotation filters
VCF_POST_ANN_FILTER_EXPRESSION ?= ExAC_AF > 0.1
vcf/ft2/cft/%.cft.vcf : vcf/ann/%.vcf
	$(call LSCRIPT_CHECK_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ \
		--filterExpression '$(VCF_POST_ANN_FILTER_EXPRESSION)' --filterName customFilter")

COMMON_FILTER_VCF = $(PYTHON) modules/vcf_tools/common_filter_vcf.py
vcf/ft2/common_ft/%.common_ft.vcf : vcf/ann/%.vcf
	$(call LSCRIPT_MEM,4G,5G,"$(COMMON_FILTER_VCF) $< $@")
