include modules/Makefile.inc

LOGDIR ?= log/varscanTN.$(NOW)

IGNORE_FP_FILTER ?= true
VALIDATION ?= false
FP_FILTER = $(PERL) $(HOME)/share/usr/bin/fpfilter.pl
BAM_READCOUNT = $(HOME)/share/usr/bin/bam-readcount
VARSCAN_TO_VCF = $(PERL) modules/variant_callers/somatic/varscanTNtoVcf.pl
MIN_MAP_QUAL ?= 1
MIN_VAR_FREQ ?= $(if $(findstring false,$(VALIDATION)),0.05,0.000001)

VARSCAN_MEM = $(JAVA7) -Xmx$1 -jar $(VARSCAN_JAR)
VARSCAN = $(call VARSCAN_MEM,8G)
VARSCAN_OPTS = $(if $(findstring true,$(VALIDATION)),--validation 1 --strand-filter 0) --min-var-freq $(MIN_VAR_FREQ)
VARSCAN_SOURCE_ANN_VCF = python modules/vcf_tools/annotate_source_vcf.py --source varscan
VPATH ?= bam
VARSCAN_VARIANT_TYPES = varscan_indels varscan_snps

varscan : $(foreach chr,$(CHROMOSOMES),$(foreach pair,$(SAMPLE_PAIRS),varscan/chr_tables/$(pair).$(chr).varscan_timestamp)) \
	  $(foreach chr,$(CHROMOSOMES),$(foreach pair,$(SAMPLE_PAIRS),varscan/chr_tables/$(pair).$(chr).snp.txt)) \
	  $(foreach chr,$(CHROMOSOMES),$(foreach pair,$(SAMPLE_PAIRS),varscan/chr_tables/$(pair).$(chr).indel.txt)) \
	  $(foreach pair,$(SAMPLE_PAIRS),varscan/tables/$(pair).snp.txt) \
	  $(foreach pair,$(SAMPLE_PAIRS),varscan/tables/$(pair).indel.txt) \
	  $(foreach pair,$(SAMPLE_PAIRS),varscan/tables/$(pair).snp.Somatic.txt) \
	  $(foreach pair,$(SAMPLE_PAIRS),varscan/tables/$(pair).indel.Somatic.txt) \
	  $(foreach pair,$(SAMPLE_PAIRS),varscan/vcf/$(pair).snp.Somatic.vcf) \
	  $(foreach pair,$(SAMPLE_PAIRS),varscan/vcf/$(pair).indel.Somatic.vcf) \
	  $(foreach type,$(VARSCAN_VARIANT_TYPES),$(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).$(type).vcf))

define varscan-somatic-tumor-normal-chr
varscan/chr_tables/$1_$2.$3.varscan_timestamp : bam/$1.bam bam/$2.bam bam/$1.bam.bai bam/$2.bam.bai
	if [[ $$$$($$(SAMTOOLS) view $$< $3 | head -1 | wc -l) -gt 0 ]]; then \
		$$(call RUN,-s 9G -m 12G -w 72:00:00,"set -o pipefail && \
		rm -rf varscan/chr_tables/$1_$2.$3.snp.txt && \
		rm -rf varscan/chr_tables/$1_$2.$3.indel.txt && \
		$$(VARSCAN) somatic \
		<($$(SAMTOOLS) mpileup -A -r $3 -q $$(MIN_MAP_QUAL) -f $$(REF_FASTA) $$(word 2,$$^)) \
		<($$(SAMTOOLS) mpileup -A -r $3 -q $$(MIN_MAP_QUAL) -f $$(REF_FASTA) $$<) \
		$$(VARSCAN_OPTS) \
		--output-indel varscan/chr_tables/$1_$2.$3.indel.txt --output-snp varscan/chr_tables/$1_$2.$3.snp.txt && touch $$@"); \
	else \
		$$(INIT) \
		echo 'chrom	position	ref	var	normal_reads1	normal_reads2	normal_var_freq	normal_gt	tumor_reads1	tumor_reads2	tumor_var_freq	tumor_gt	somatic_status	variant_p_value	somatic_p_value	tumor_reads1_plus	tumor_reads1_minus	tumor_reads2_plus	tumor_reads2_minus	normal_reads1_plus	normal_reads1_minus	normal_reads2_plus	normal_reads2_minus' > varscan/chr_tables/$1_$2.$3.indel.txt; \
		echo 'chrom	position	ref	var	normal_reads1	normal_reads2	normal_var_freq	normal_gt	tumor_reads1	tumor_reads2	tumor_var_freq	tumor_gt	somatic_status	variant_p_value	somatic_p_value	tumor_reads1_plus	tumor_reads1_minus	tumor_reads2_plus	tumor_reads2_minus	normal_reads1_plus	normal_reads1_minus	normal_reads2_plus	normal_reads2_minus' > varscan/chr_tables/$1_$2.$3.snp.txt; \
		touch $$@; \
	fi

varscan/chr_tables/$1_$2.$3.indel.txt : varscan/chr_tables/$1_$2.$3.varscan_timestamp

varscan/chr_tables/$1_$2.$3.snp.txt : varscan/chr_tables/$1_$2.$3.varscan_timestamp

varscan/chr_tables/$1_$2.$3.%.fp_pass.txt : varscan/chr_tables/$1_$2.$3.%.txt bamrc/$1.$3.bamrc.gz
	$$(call RUN,-s 8G -m 55G,"$$(VARSCAN) fpfilter $$< <(zcat $$(<<)) --output-file $$@")

endef
$(foreach chr,$(CHROMOSOMES), \
	$(foreach pair,$(SAMPLE_PAIRS), \
		$(eval $(call varscan-somatic-tumor-normal-chr,$(tumor.$(pair)),$(normal.$(pair)),$(chr)))))

define merge-varscan-pair-type
varscan/tables/$1.$2.txt : $$(foreach chr,$$(CHROMOSOMES),\
	$$(if $$(findstring true,$$(VALIDATION) $$(IGNORE_FP_FILTER)),\
	varscan/chr_tables/$1.$$(chr).$2.txt,\
	varscan/chr_tables/$1.$$(chr).$2.fp_pass.txt))
	$$(INIT) head -1 $$< > $$@ && for x in $$^; do sed 1d $$$$x >> $$@; done

endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(foreach type,snp indel, \
		$(eval $(call merge-varscan-pair-type,$(pair),$(type)))))
	
define filter-varscan-pair-type
varscan/tables/$1.$2.Somatic.txt : varscan/tables/$1.$2.txt
	$$(call RUN,-s 5G -m 8G,"set -o pipefail && \
				$$(VARSCAN) somaticFilter $$(<) && \
				$$(VARSCAN) processSomatic $$(<) && \
				cp varscan/tables/$1.$2.txt.Somatic varscan/tables/$1.$2.Somatic.txt && \
				rm varscan/tables/$1.$2.txt.Somatic && \
				cp varscan/tables/$1.$2.txt.Germline varscan/tables/$1.$2.Germline.txt && \
				rm varscan/tables/$1.$2.txt.Germline && \
				cp varscan/tables/$1.$2.txt.LOH varscan/tables/$1.$2.LOH.txt && \
				rm varscan/tables/$1.$2.txt.LOH && \
				cp varscan/tables/$1.$2.txt.hc varscan/tables/$1.$2.hc.txt && \
				rm varscan/tables/$1.$2.txt.hc")

endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(foreach type,snp indel, \
		$(eval $(call filter-varscan-pair-type,$(pair),$(type)))))

define convert-varscan-tumor-normal
varscan/vcf/$1_$2.$3.Somatic.vcf : varscan/tables/$1_$2.$3.Somatic.txt
	$$(call RUN,-s 4G -m 8G,"set -o pipefail && \
				 $$(VARSCAN_TO_VCF) -f $$(REF_FASTA) -t $1 -n $2 $$(<) | $$(VCF_SORT) $$(REF_DICT) - > $$(@)")


endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(foreach type,snp indel, \
		$(eval $(call convert-varscan-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)),$(type)))))

vcf/%.varscan_indels.vcf : varscan/vcf/%.indel.Somatic.vcf
	$(INIT) $(VARSCAN_SOURCE_ANN_VCF) < $< > $@

vcf/%.varscan_snps.vcf : varscan/vcf/%.snp.Somatic.vcf
	$(INIT) $(VARSCAN_SOURCE_ANN_VCF) < $< > $@

include modules/variant_callers/gatk.mk

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY: varscan
