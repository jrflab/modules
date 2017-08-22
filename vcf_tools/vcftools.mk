# vim: set ft=make :
# sub module containing vcf related tools

ifndef VCFTOOLS_MK

include modules/Makefile.inc
#include modules/variant_callers/gatk.inc

LOGDIR ?= log/vcf.$(NOW)


..DUMMY := $(shell mkdir -p version; echo "$(SNP_EFF) &> version/snp_eff.txt")

# flags for non-gatk snp eff
SNP_SIFT_OPTS = -c $(SNP_EFF_CONFIG)

MUT_ASS = $(RSCRIPT) modules/vcf_tools/mutAssVcf.R

%.vcf.idx : %.vcf
	$(call RUN,-c -s 4G -m 8G,"$(IGVTOOLS) index $< && sleep 10")

%.pass.vcf : %.vcf
	$(call CHECK_VCF,$(call RUN,-c -s 8G -m 12G,"$(call SNP_SIFT_MEM,7G) filter $(SNP_SIFT_OPTS) \
		-f $< \"( na FILTER ) | (FILTER = 'PASS')\" > $@.tmp && $(call VERIFY_VCF,$@.tmp,$@)"))

define rename-samples-tumor-normal
vcf/$1_$2.%.rn.vcf : vcf/$1_$2.%.vcf
	$$(INIT) perl -ne 'if (/^#CHROM/) { s/NORMAL/$2/; s/TUMOR/$1/; } print;' $$< > $$@
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call rename-samples-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

# VariantEval: generate vcf report
reports/%.grp : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).%.vcf)
	$(call RUN,-s 2G -m 5G,"$(call GATK_MEM,2G) -T VariantEval $(foreach sm,$(REPORT_STRATIFICATION), --stratificationModule $(sm)) -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach eval,$(filter %.vcf,$^), --eval:$(call strip-suffix,$(notdir $(eval))) $(eval)) -o $@")

%.norm.vcf.gz : %.vcf
	$(call RUN,-s 9G -m 12G,"sed '/^##GATKCommandLine/d;/^##MuTect/d;' $< | \
		$(VT) view -h -f PASS - | \
		$(VT) decompose -s - | \
		$(VT) normalize -r $(REF_FASTA) - | \
		$(call SNP_EFF_MEM,8G) ann -c $(SNP_EFF_CONFIG) $(SNP_EFF_GENOME) -formatEff -classic | \
		bgzip -c > $@.tmp && if zgrep -q '^#CHROM' $@.tmp; then mv $@.tmp $@; else false; fi")

%.vcf.gz.tbi : %.vcf.gz
	$(call RUN,-c -w 1:00:00 -s 3G -m 5G,"$(BCFTOOLS2) index -t -f $<")

%.vcf.gz : %.vcf
	$(call RUN,-c -s 2G -m 3G,"bgzip -c $< > $@")

define vcf2maf-tumor-normal
maf/$1_$2.%.maf : vcf_ann/$1_$2.%.vcf
	$$(call RUN,-c -v $$(VEP_ENV) -s 9G -m 12G,"$$(VCF2MAF) --input-vcf $$< --tumor-id $1 --normal-id $2 --filter-vcf $$(EXAC_NONTCGA) --ref-fasta $$(REF_FASTA) --tmp-dir `mktemp -d` --vep-path $$(VEP_PATH) --vep-data $$(VEP_DATA) --output-maf $$@.tmp && sed '/^#/d' $$@.tmp > $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call vcf2maf-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

maf/allTN.%.maf : $(foreach pair,$(SAMPLE_PAIRS),maf/$(pair).%.maf)
	$(INIT) csvstack -t $^ | csvformat -T > $@

define vcf2maf-sample
maf/$1.%.maf : vcf_ann/$1.%.vcf
	$$(call RUN,-s 9G -m 12G,"$$(VCF2MAF) --input-vcf $$< --tumor-id $1 $$(if $$(EXAC_NONTCGA),--filter-vcf $$(EXAC_NONTCGA)) --ref-fasta $$(REF_FASTA) --vep-path $$(VEP_PATH) --vep-data $$(VEP_DATA) --tmp-dir `mktemp -d` --output-maf $$@")
endef
$(foreach sample,$(SAMPLES),$(eval $(call vcf2maf-sample,$(sample))))

maf/all.%.maf : $(foreach sample,$(SAMPLES),maf/$(sample).%.maf)
	$(INIT) csvstack -t $^ | csvformat -T > $@

ifdef SAMPLE_PAIRS
alltables/allTN.%.txt : $(foreach pair,$(SAMPLE_PAIRS),tables/$(pair).%.txt)
	$(call RUN,-s 5G -m 12G,"$(RSCRIPT) $(RBIND) --tumorNormal $^ > $@")
endif

# extract vcf to table
VCF_FIELDS = CHROM POS ID REF ALT FILTER
ANN_FIELDS = $(addprefix ANN[*].,ALLELE EFFECT IMPACT GENE GENEID FEATURE FEATUREID BIOTYPE RANK HGVS_C HGVS_P CDNA_POS CDNA_LEN CDS_POS CDS_LEN AA_POS AA_LEN DISTANCE ERRORS)
tables/%.opl_tab.txt : vcf_ann/%.vcf
	$(call RUN,-c -s 9G -m 15G,"format_fields=\$$(grep '^##FORMAT=<ID=' $< | sed 's/CGC_Other_Syndrome\/Disease/CGC_Other_Syndrome_or_Disease/; s/GERP++/GERPpp/; s/(/_/; s/)//; s/.*ID=//; s/,.*//;' | tr '\n' ' '); \
	N=\$$(expr \$$(grep '^#CHROM' $< | wc -w) - 10); \
	fields='$(VCF_FIELDS)'; \
	for f in \$$format_fields; do \
		for i in \$$(seq 0 \$$N); do \
			fields+=' 'GEN[\$$i].\$$f; \
		done; \
	done; \
	fields+=' '\$$(grep '^##INFO=<ID=' $< | grep -v '=REF,' | sed 's/CGC_Other_Syndrome\/Disease/CGC_Other_Syndrome_or_Disease/; s/GERP++/GERPpp/; s/(/_/; s/)//; s/.*ID=//; s/,.*//; s/\bANN\b/$(ANN_FIELDS)/; ' | tr '\n' ' '); \
	sed 's/CGC_Other_Syndrome\/Disease/CGC_Other_Syndrome_or_Disease/; s/GERP++/GERPpp/; s/(/_/; s/)// ' $< | $(VCF_EFF_ONE_PER_LINE) | $(call SNP_SIFT_MEM,6G) extractFields - \$$fields > $@.tmp && \
	if grep -q '^CHROM' $@.tmp; then mv $@.tmp $@; else false; fi && \
	for i in \`seq 0 \$$N\`; do \
	S=\$$(grep '^#CHROM' $< | cut -f \$$((\$$i + 10))); \
	sed -i \"1s/GEN\[\$$i\]/\$$S/g;\" $@; \
	done")

%.tab.txt : %.opl_tab.txt
	$(call RUN,-s 8G -m 20G,"$(PERL) $(VCF_JOIN_EFF) < $< > $@")

%.high_moderate.txt : %.txt
	col=$$(head -1 $< | tr '\t' '\n' | grep -n "IMPACT" | sed 's/:.*//'); \
	$(INIT) head -1 $< > $@ && awk -v col=$$col 'match($$col, /MODERATE/) || match($$col, /HIGH/)' $< >> $@

%.low_modifier.txt : %.txt
	col=$$(head -1 $< | tr '\t' '\n' | grep -n "IMPACT" | sed 's/:.*//'); \
	$(INIT) head -1 $< > $@ && awk -v col=$$col '! (match($$col, /MODERATE/) || match($$col, /HIGH/)) && (match($$col, /LOW/) || match($$col,/MODIFIER/))' $< >> $@

%.synonymous.txt : %.txt
	col_imp=$$(head -1 $< | tr '\t' '\n' | grep -n "IMPACT" | sed 's/:.*//'); \
	col_eff=$$(head -1 $< | tr '\t' '\n' | grep -n "EFFECT" | sed 's/:.*//'); \
	$(INIT) head -1 $< > $@ && awk -v col_imp=$$col_imp -v col_eff=$$col_eff '! (match($$col_imp, /MODERATE/) || match($$col_imp, /HIGH/)) && (match($$col_imp, /LOW/) && (match($$col_eff, /synonymous_variant/)))' $< >> $@

%.nonsynonymous.txt : %.txt
	col=$$(head -1 $< | tr '\t' '\n' | grep -n "IMPACT" | sed 's/:.*//'); \
	$(INIT) head -1 $< > $@ && awk -v col=$$col 'match($$col, /MODERATE/) || match($$col, /HIGH/)' $< >> $@

include modules/vcf_tools/vcfPostFilters.mk
include modules/vcf_tools/vcfFilters.mk
include modules/vcf_tools/vcfAnnotations.mk
include modules/vcf_tools/vcfPostAnnotations.mk

endif
VCFTOOLS_MK = true

