# vim: set ft=make :
# sub module containing vcf related tools

ifndef VCFTOOLS_MK

include modules/Makefile.inc
#include modules/variant_callers/gatk.inc

LOGDIR ?= log/vcf.$(NOW)

..DUMMY := $(shell mkdir -p version; echo "$(SNP_EFF) &> version/snp_eff.txt")

ANNOVAR = $(PERL) $(HOME)/share/usr/annovar/table_annovar.pl
ANNOVAR_PROTOCOL ?= refGene$(,)cytoBand$(,)genomicSuperDups$(,)esp6500siv2_all$(,)1000g2014oct_all$(,)1000g2014oct_afr$(,)1000g2014oct_eas$(,)1000g2014oct_eur$(,)snp138$(,)ljb26_all
ANNOVAR_OPERATION ?= g$(,)r$(,)r$(,)f$(,)f$(,)f$(,)f$(,)f$(,)f$(,)f
ANNOVAR_OPTS = --dot2underline -remove -protocol $(ANNOVAR_PROTOCOL) -operation $(ANNOVAR_OPERATION) -nastring . -vcfinput -buildver $(ANNOVAR_REF)

# flags for non-gatk snp eff
VCF_FIELDS = CHROM POS ID REF ALT FILTER
ANN_FIELDS = $(addprefix ANN[*].,ALLELE EFFECT IMPACT GENE GENEID FEATURE FEATUREID BIOTYPE RANK HGVS_C HGVS_P CDNA_POS CDNA_LEN CDS_POS CDS_LEN AA_POS AA_LEN DISTANCE ERRORS)
SNP_EFF_FLAGS ?= -canon # -ud 0  -no-intron -no-intergenic -no-utr
SNP_EFF_OPTS = -c $(SNP_EFF_CONFIG) -i vcf -o vcf $(SNP_EFF_FLAGS)
DEPTH_FILTER ?= 5
SNP_SIFT_OPTS = -c $(SNP_EFF_CONFIG)

CHASM = $(RSCRIPT) modules/vcf_tools/chasmVcf.R 
#CHASM_DIR = /ifs/opt/common/CHASM/CHASMDL.1.0.7
CHASM_DIR = $(HOME)/share/usr/CHASM
CHASM_PYTHON_ENV = $(HOME)/share/usr/anaconda-envs/pyenv27-chasm
CHASM_CLASSIFIER ?= Breast

FATHMM = $(MY_RSCRIPT) modules/vcf_tools/fathmmVcf.R 
FATHMM_DIR = $(HOME)/share/usr/fathmm
FATHMM_PYTHON = $(HOME)/share/usr/bin/python
FATHMM_PYTHONPATH = $(HOME)/share/usr/lib/python:$(HOME)/share/usr/lib/python2.7
FATHMM_OPTS = --genome $(REF) --ensemblTxdb $(ENSEMBL_TXDB) --ref $(REF_FASTA) --fathmmDir $(FATHMM_DIR) --python $(FATHMM_PYTHON) \
			  --mysqlHost $(EMBL_MYSQLDB_HOST) --mysqlPort $(EMBL_MYSQLDB_PORT) --mysqlUser $(EMBL_MYSQLDB_USER) --mysqlPassword $(EMBL_MYSQLDB_PW) --mysqlDb $(EMBL_MYSQLDB_DB)

TRANSFIC = $(RSCRIPT) modules/vcf_tools/transficVcf.R
TRANSFIC_PERL_SCRIPT = $(HOME)/share/usr/transfic/bin/transf_scores.pl

MUT_ASS = $(RSCRIPT) modules/vcf_tools/mutAssVcf.R

ADD_GENE_LIST_ANNOTATION = $(RSCRIPT) modules/vcf_tools/addGeneListAnnotationToVcf.R

NON_SILENT_EFF = START_GAINED SPLICE_SITE_ACCEPTOR SPLICE_SITE_DONOR START_LOST NON_SYNONYMOUS_CODING FRAME_SHIFT CODON_CHANGE CODON_INSERTION CODON_CHANGE_PLUS_CODON_INSERTION CODON_DELETION CODON_CHANGE_PLUS_CODON_DELETION STOP_GAINED STOP_LOST NON_SYNONYMOUS_START
NON_SILENT_CODING_EFF = START_GAINED START_LOST NON_SYNONYMOUS_CODING FRAME_SHIFT CODON_CHANGE CODON_INSERTION CODON_CHANGE_PLUS_CODON_INSERTION CODON_DELETION CODON_CHANGE_PLUS_CODON_DELETION STOP_GAINED STOP_LOST NON_SYNONYMOUS_START

NORMAL_VCF ?= $(HOME)/share/reference/spowellnormal.gatk_variants.vcf
ifdef NORMAL_VCF
%.nft.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ --maskName 'normal' --mask $(NORMAL_VCF) && $(RM) $< $<.idx")
endif

# run snp eff
%.eff.vcf : %.vcf %.vcf.idx
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,9G,14G,"$(call SNP_EFF_MEM,8G) ann $(SNP_EFF_OPTS) $(SNP_EFF_GENOME) -s $*.eff_summary.html $< > $@ && $(RM) $^"))

# run snp sift to annotated with dbnsfp
%.nsfp.vcf : %.vcf %.vcf.idx
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,9G,12G,"$(call SNP_SIFT_MEM,8G) dbnsfp $(SNP_SIFT_OPTS) -f $(subst $( ),$(,),$(NSFP_FIELDS)) -db $(DB_NSFP) $< | sed '/^##INFO=<ID=dbNSFP/ s/Character/String/' > $@ && $(RM) $^"))

# run gatk snp eff
%.gatk_eff.vcf : %.vcf %.vcf.idx
	$(call LSCRIPT_MEM,5G,8G,"$(call SNP_EFF_MEM,4G) eff -i vcf -o gatk $(SNP_EFF_GENOME) $< > $@")

# process snp eff output with gatk %=sample.indels/snps
%.annotated.vcf : %.vcf %.gatk_eff.vcf %.gatk_eff.vcf.idx %.vcf.idx 
	$(call LSCRIPT_PARALLEL_MEM,5,2G,3G,"$(call GATK_MEM,8G) -T VariantAnnotator \
	-R $(REF_FASTA) -nt 5 -A SnpEff  --variant $<  --snpEffFile $(word 2,$^) -o $@ &> $(LOGDIR)/$@.log")

# snp sift using GMAF > 1% filtered dbsnp 
%.dbsnp.vcf : %.vcf %.vcf.idx 
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,9G,12G,"$(call SNP_SIFT_MEM,8G) annotate $(SNP_SIFT_OPTS) $(DBSNP1PC) $< > $@ && $(RM) $^"))

# mouse genome project dbsnp
%.mgp_dbsnp.vcf : %.vcf %.vcf.idx 
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,33G,65G,"$(call SNP_SIFT_MEM,45G) annotate -tabix $(SNP_SIFT_OPTS) $(MGP_SNP_DBSNP) $< | $(call SNP_SIFT_MEM,10G) annotate -tabix $(SNP_SIFT_OPTS) $(MGP_INDEL_DBSNP) > $@ && $(RM) $^"))

%.cosmic.vcf : %.vcf %.vcf.idx 
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,9G,12G,"$(call SNP_SIFT_MEM,8G) annotate $(SNP_SIFT_OPTS) $(COSMIC) $< > $@ && $(RM) $^"))

%.clinvar.vcf : %.vcf %.vcf.idx 
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,9G,12G,"$(call SNP_SIFT_MEM,8G) annotate $(SNP_SIFT_OPTS) $(CLINVAR) $< > $@ && $(RM) $^"))

%.exac_nontcga.vcf : %.vcf %.vcf.idx 
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,9G,12G,"$(call SNP_SIFT_MEM,8G) annotate $(SNP_SIFT_OPTS) -info ExAC_AF $(EXAC_NONTCGA) $< > $@ && $(RM) $^"))


# apply overall depth filter
%.dp_ft.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ --filterExpression 'DP < $(DEPTH_FILTER)' --filterName Depth && $(RM) $<")

%.sdp_ft.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,2G,5G,"$(call SNP_SIFT_MEM,2G) filter $(SNP_SIFT_OPTS) -f $< '(exists GEN[*].DP) & (GEN[*].DP > 20)' > $@"))

# apply HRun filter
%.hrun_ft.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ --filterExpression 'HRun > $(HRUN_FILTER)' --filterName HRun && $(RM) $< $<.idx")

%.pass.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,2G,5G,"$(call SNP_SIFT_MEM,2G) filter $(SNP_SIFT_OPTS) -f $< \"( na FILTER ) | (FILTER = 'PASS')\" > $@"))

# apply dp filter for somatic sniper
%.ss_dp_ft.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ --filterExpression 'vc.getGenotype(\"TUMOR\").getDP() < $(DEPTH_FILTER) || vc.getGenotype(\"NORMAL\").getDP() < $(DEPTH_FILTER)' --filterName depthFilter && $(RM) $< $<.idx")

# somatic sniper somatic flag filter
%.ss_ft.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ --filterExpression 'vc.getGenotype(\"TUMOR\").getAttributeAsInt(\"SS\", 0) != 2'  --filterName nonSomatic && $(RM) $< $<.idx")

# target region filter
%.target_ft.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,8G,12G,"$(call GATK_MEM2,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ --mask $(TARGETS_FILE) --maskName targetInterval --filterNotInMask && $(RM) $< $<.idx")


# varscan TN variant allele frequency: min tumor freq > 5% ; max normal freq < 5%
%.freq_ft.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,2G,5G,"sed '/##FORMAT=<ID=FREQ/ s/String/Float/; /^#/! s/%//g' $< | $(call SNP_SIFT_MEM,2G) filter $(SNP_SIFT_OPTS) '(exists GEN[*].FREQ) & (GEN[0].FREQ < 5) & (GEN[1].FREQ[0] > 5)' > $@ && $(RM) $< $<.idx")

%.vaf_ft.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,2G,5G,"$(call SNP_SIFT_MEM,2G) filter $(SNP_SIFT_OPTS) '(exists GEN[*].VAF) & (GEN[0].VAF > 0.05) & (GEN[1].VAF < 0.05)' < $< > $@ && $(RM) $< $<.idx")


# varscan depth filter (b/c varscan is dumb and only gives variant depth)
%.vdp_ft.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,2G,5G,"cat $< | $(call SNP_SIFT_MEM,2G) filter $(SNP_SIFT_OPTS) '(exists GEN[*].AD) & (GEN[*].AD > $(DEPTH_FILTER))' > $@ && $(RM) $< $<.idx")

# add exon distance
%.exondist.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,2G,3G,"$(INTRON_POSN_LOOKUP) $< > $@")

%.vcf.idx : %.vcf
	$(call LSCRIPT_CHECK_MEM,4G,8G,"$(IGVTOOLS) index $< && sleep 10")

%.chasm.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,8G,17G,"unset PYTHONPATH && source $(CHASM_PYTHON_ENV)/bin/activate $(CHASM_PYTHON_ENV) && $(CHASM) --genome $(REF) --classifier $(subst $( ),$(,),$(CHASM_CLASSIFIER)) --chasmDir $(CHASM_DIR) --python $(shell which python) --outFile $@ $< && $(RM) $< $<.idx"))

%.fathmm.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,8G,10G,"PYTHONPATH=$(FATHMM_PYTHONPATH) $(FATHMM) $(FATHMM_OPTS) --outFile $@ $< && $(RM) $< $<.idx"))

%.mutass.vcf : %.vcf
	$(call LSCRIPT_MEM,12G,15G,$(MUT_ASS) --outFile $@ --maData $(MUT_ASS_RDATA) $<)

%.transfic.vcf : %.vcf
	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_MEM,9G,12G,"$(TRANSFIC) --genome $(REF) --transfic $(TRANSFIC_PERL_SCRIPT) --outFile $@ $< && $(RM) $< $<.idx"))

ifdef SAMPLE_SET_PAIRS
define somatic-filter-vcf-set
vcf/$1.%.som_ft.vcf : vcf/$1.%.vcf
	$$(INIT) $$(SOMATIC_FILTER_VCF) -n $(normal.$1) -f 0.03 $$< > $$@ 2> $$(LOG) && $$(RM) $$< $$<.idx
endef
$(foreach set,$(SAMPLE_SET_PAIRS),$(eval $(call somatic-filter-vcf-set,$(set))))
endif

ifdef SAMPLE_PAIRS
# ff normal filter :
# filter if normal depth > 20 and normal variant depth > 1/5 * tumor variant depth
# or normal variant depth greater than 1
# ffpe normal filter :
# filter if normal depth > 20 and normal variant depth > 1/3 * tumor variant depth
# or normal variant depth greater than 1
define som-ad-ft-tumor-normal
vcf/$1_$2.%.som_ad_ft.vcf : vcf/$1_$2.%.vcf
	$$(call LSCRIPT_CHECK_MEM,8G,12G,"$$(call GATK_MEM,8G) -T VariantFiltration -R $$(REF_FASTA) -V $$< -o $$@ \
		--filterExpression 'vc.getGenotype(\"$1\").getAD().1 < $(DEPTH_FILTER)' \
		--filterName tumorVarAlleleDepth \
		--filterExpression 'if (vc.getGenotype(\"$2\").getDP() > 20) { ( vc.getGenotype(\"$2\").getAD().1 * 1.0 / vc.getGenotype(\"$2\").getDP()) > ( vc.getGenotype(\"$1\").getAD().1 * 1.0 / vc.getGenotype(\"$1\").getDP()) / 5.0 } else { vc.getGenotype(\"$2\").getAD().1 > 1 }' \
		--filterName somaticAlleleDepth \
		--filterExpression 'vc.getGenotype(\"$1\").getDP() <= $$(DEPTH_FILTER) || vc.getGenotype(\"$2\").getDP() <= $$(DEPTH_FILTER)' \
		--filterName depthFilter && $$(RM) $$< $$<.idx")

vcf/$1_$2.%.ffpe_som_ad_ft.vcf : vcf/$1_$2.%.vcf
	$$(call LSCRIPT_CHECK_MEM,8G,12G,"$$(call GATK_MEM,8G) -T VariantFiltration -R $$(REF_FASTA) -V $$< -o $$@ \
		--filterExpression 'vc.getGenotype(\"$1\").getAD().1 < $(DEPTH_FILTER)' \
		--filterName tumorVarAlleleDepth \
		--filterExpression 'if (vc.getGenotype(\"$2\").getDP() > 20) { ( vc.getGenotype(\"$2\").getAD().1 * 1.0 / vc.getGenotype(\"$2\").getDP()) > ( vc.getGenotype(\"$1\").getAD().1 * 1.0 / vc.getGenotype(\"$1\").getDP()) / 3.0 } else { vc.getGenotype(\"$2\").getAD().1 > 1 }' \
		--filterName somaticAlleleDepth \
		--filterExpression 'vc.getGenotype(\"$1\").getDP() <= $$(DEPTH_FILTER) || vc.getGenotype(\"$2\").getDP() <= $$(DEPTH_FILTER)' \
		--filterName depthFilter && $$(RM) $$< $$<.idx")

# somatic filter for structural variants
vcf/$1_$2.%.sv_som_ft.vcf : vcf/$1_$2.%.vcf
	$$(call LSCRIPT_CHECK_MEM,8G,12G,"$$(call GATK_MEM,8G) -T VariantFiltration -R $$(REF_FASTA) -V $$< -o $$@ \
		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"SU\") <= $$(DEPTH_FILTER)' \
		--filterName svSupport \
		--filterExpression 'vc.getGenotype(\"$1\").getAnyAttribute(\"SU\") < 5 * vc.getGenotype(\"$2\").getAnyAttribute(\"SU\")' \
		--filterName somaticSvSupport \
		&& $$(RM) $$< $$<.idx")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call som-ad-ft-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))


define rename-samples-tumor-normal
vcf/$1_$2.%.rn.vcf : vcf/$1_$2.%.vcf
	$$(INIT) perl -ne 'if (/^#CHROM/) { s/NORMAL/$2/; s/TUMOR/$1/; } print;' $$< > $$@ && $$(RM) $$< $$<.idx
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call rename-samples-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))


define ad-tumor-normal
vcf/$1_$2.%.ad.vcf : vcf/$1_$2.%.vcf bam/$1.bam bam/$2.bam bam/$1.bai bam/$2.bai
	$$(call LSCRIPT_CHECK_PARALLEL_MEM,4,2G,3G,"$$(call GATK_MEM,8G) -T VariantAnnotator -nt 4 -R $$(REF_FASTA) -A DepthPerAlleleBySample --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -o $$@ -L $$<")
endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call ad-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))

define annotate-tumor-normal
vcf/$1_$2.%.ann.vcf : vcf/$1_$2.%.vcf bam/$1.bam bam/$2.bam bam/$1.bai bam/$2.bai
	$$(call LSCRIPT_CHECK_PARALLEL_MEM,4,2G,3G,"$$(call GATK_MEM,8G) -T VariantAnnotator -nt 4 -R $$(REF_FASTA) $$(foreach ann,$$(VCF_ANNOTATIONS),-A $$(ann) ) --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -o $$@ -L $$< && $$(RM) $$< $$<.idx")
endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call annotate-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))

define hrun-tumor-normal
vcf/$1_$2.%.hrun.vcf : vcf/$1_$2.%.vcf bam/$1.bam bam/$2.bam bam/$1.bai bam/$2.bai
	$$(call LSCRIPT_CHECK_PARALLEL_MEM,4,2G,3G,"$$(call GATK_MEM,8G) -T VariantAnnotator -nt 4 -R $$(REF_FASTA) -A HomopolymerRun --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -L $$< -o $$@ && $$(RM) $$< $$<.idx")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call hrun-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))
endif

define annotate-sample
vcf/$1.%.ann.vcf : vcf/$1.%.vcf bam/$1.bam bam/$1.bai
	$$(call LSCRIPT_PARALLEL_MEM,4,2G,3G,"$$(call GATK_MEM,8G) -T VariantAnnotator -nt 4 -R $$(REF_FASTA) $$(foreach ann,$$(VCF_ANNOTATIONS),-A $$(ann) ) --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -L $$< -V $$< -o $$@ && $$(RM) $$< $$<.idx")
endef
$(foreach sample,$(SAMPLES),$(eval $(call annotate-sample,$(sample))))

define hrun-sample
vcf/$1.%.hrun.vcf : vcf/$1.%.vcf bam/$1.bam bam/$1.bai
	$$(call LSCRIPT_CHECK_PARALLEL_MEM,4,2G,3G,"$$(call GATK_MEM,8G) -T VariantAnnotator -nt 4 -R $$(REF_FASTA) -L $$< -A HomopolymerRun --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -o $$@ && $$(RM) $$< $$<.idx")
endef
$(foreach sample,$(SAMPLES),$(eval $(call hrun-sample,$(sample))))

# VariantEval: generate vcf report
reports/%.grp : $(foreach sample,$(SAMPLES),vcf/$(sample).%.vcf) $(foreach sample,$(SAMPLES),vcf/$(sample).%.vcf.idx)
	$(call LSCRIPT_MEM,2G,5G,"$(call GATK_MEM,2G) -T VariantEval $(foreach sm,$(REPORT_STRATIFICATION), --stratificationModule $(sm)) -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach eval,$(filter %.vcf,$^), --eval:$(call strip-suffix,$(notdir $(eval))) $(eval)) -o $@")
ifdef SAMPLE_PAIRS
reports/%.grp : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).%.vcf vcf/$(pair).%.vcf.idx)
	$(call LSCRIPT_MEM,2G,5G,"$(call GATK_MEM,2G) -T VariantEval $(foreach sm,$(REPORT_STRATIFICATION), --stratificationModule $(sm)) -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach eval,$(filter %.vcf,$^), --eval:$(call strip-suffix,$(notdir $(eval))) $(eval)) -o $@")
endif

# merge tables
alltables/all.%.txt : $(foreach sample,$(SAMPLES),tables/$(sample).%.txt)
	$(call LSCRIPT_MEM,5G,12G,"$(RBIND) --sampleName $< $^ > $@")
ifdef SAMPLE_SETS
alltables/allSS.%.txt : $(foreach set,$(SAMPLE_SETS),tables/$(set).%.txt)
	$(call LSCRIPT_MEM,5G,12G,"$(RSCRIPT) $(RBIND) --normalLast $^ > $@")
endif
ifdef SAMPLE_PAIRS
alltables/allTN.%.txt : $(foreach pair,$(SAMPLE_PAIRS),tables/$(pair).%.txt)
	$(call LSCRIPT_MEM,5G,12G,"$(RSCRIPT) $(RBIND) --tumorNormal $^ > $@")
endif

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

# DEPRECATED
#%.nonsilent.txt : %.txt
#	$(INIT) head -1 $< > $@ && sed '1d' $< | grep $(foreach eff,$(NON_SILENT_EFF), -e $(eff)) >> $@ || true
#
#%.nonsilent_cds.txt : %.txt
#	$(INIT) head -1 $< > $@ && sed '1d' $< | grep $(foreach eff,$(NON_SILENT_CODING_EFF), -e $(eff)) >> $@ || true
#
#%.missense.txt : %.txt
#	$(INIT) head -1 $< > $@ && sed '1d' $< | grep -e MISSENSE >> $@ || true
#
#%.silent.txt : %.txt
#	$(INIT) head -1 $< > $@ && sed '1d' $< | grep -e SILENT >> $@ || true

# extract vcf to table
tables/%.opl_tab.txt : vcf/%.vcf
	$(call LSCRIPT_CHECK_MEM,2G,5G,"format_fields=\$$(grep '^##FORMAT=<ID=' $< | sed 's/.*ID=//; s/,.*//;' | tr '\n' ' '); \
	N=\$$(expr \$$(grep '^#CHROM' $< | wc -w) - 10); \
	fields='$(VCF_FIELDS)'; \
	for f in \$$format_fields; do \
		for i in \$$(seq 0 \$$N); do \
			fields+=' 'GEN[\$$i].\$$f; \
		done; \
	done; \
	fields+=' '\$$(grep '^##INFO=<ID=' $< | grep -v '=REF,' | sed 's/.*ID=//; s/,.*//; s/\bANN\b/$(ANN_FIELDS)/; ' | tr '\n' ' '); \
	$(VCF_EFF_ONE_PER_LINE) < $< | $(call SNP_SIFT_MEM,2G) extractFields - \$$fields > $@; \
	for i in \`seq 0 \$$N\`; do \
	S=\$$(grep '^#CHROM' $< | cut -f \$$((\$$i + 10))); \
	sed -i \"1s/GEN\[\$$i\]/\$$S/g;\" $@; \
	done")

%.tab.txt : %.opl_tab.txt
	$(call LSCRIPT_MEM,8G,20G,"$(PERL) $(VCF_JOIN_EFF) < $< > $@")
	
%.pass.txt : %.txt
	$(INIT) head -1 $< > $@ && awk '$$6 == "PASS" { print }' $< >> $@ || true

%.novel.txt : %.txt
	$(INIT) id=`head -1 $< | tr '\t' '\n' | grep -n "^ID$$" | sed 's/:.*//'`; \
	gmaf=`head -1 $< | tr '\t' '\n' | grep -n "^GMAF$$" | sed 's/:.*//'`; \
	awk -v id=$$id -v gmaf=$$gmaf 'NR == 1 || length($$id) == 1 || $$gmaf < 0.01 { print }' $< > $@ || true

FALSE_POSITIVE_BED = $(HOME)/share/reference/fuentes_blacklist.include_cosmic.hg19.bed
%.fp_ft.vcf : %.vcf
	$(call LSCRIPT_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ --maskName 'FuentesFalsePositive' --mask $(FALSE_POSITIVE_BED) && $(RM) $< $<.idx")

DGD_BED = $(HOME)/share/reference/dgd.include_cosmic.hg19.bed
%.dgd_ft.vcf : %.vcf
	$(call LSCRIPT_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ --maskName 'DuplicateGenesDB' --mask $(DGD_BED) && $(RM) $< $<.idx")

ENCODE_BED = $(HOME)/share/reference/wgEncodeDacMapabilityConsensusExcludable.include_cosmic.bed
%.encode_ft.vcf : %.vcf
	$(call LSCRIPT_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ --maskName 'encode' --mask $(ENCODE_BED) && $(RM) $< $<.idx")

%.het_ft.vcf : %.vcf
	$(call LSCRIPT_MEM,8G,12G,"$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) -V $< -o $@ \
		--genotypeFilterExpression 'isHet == 1' --genotypeFilterName 'Heterozygous positions'")

%.dbsnp_ft.vcf : %.vcf
	$(INIT) awk '/^#/ || $$3 ~ /^rs/ {print}' $< > $@

HAPLOTYPE_INSUF_BED = $(HOME)/share/reference/haplo_insuff_genes.bed
CANCER_GENE_CENSUS_BED = $(HOME)/share/reference/annotation_gene_lists/cancer_gene_census_genes_v20150303.bed
KANDOTH_BED = $(HOME)/share/reference/annotation_gene_lists/Kandoth_127genes.bed
LAWRENCE_BED = $(HOME)/share/reference/annotation_gene_lists/Lawrence_cancer5000-S.bed
%.gene_ann.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,8G,12G,"$(ADD_GENE_LIST_ANNOTATION) --genome $(REF) --geneBed $(HAPLOTYPE_INSUF_BED)$(,)$(CANCER_GENE_CENSUS_BED)$(,)$(KANDOTH_BED)$(,)$(LAWRENCE_BED) --name hap_insuf$(,)cancer_gene_census$(,)kandoth$(,)lawrence --outFile $@ $< && $(RM) $< $<.idx")

# Copy number regulated genes annotated per subtype
# FYI Endometrioid_MSI-L has no copy number regulated genes
CN_ENDOMETRIAL_SUBTYPES = CN_high CN_low Endometrioid_MSI_H Endometrioid_MSS Endometrioid MSI POLE Serous
CN_BREAST_SUBTYPES = ER_negative ER_positive HER2_postitive Pam50_Basal Pam50_Her2 Pam50_LumA Pam50_LumB Pam50_Normal Triple_negative
CN_ENDOMETRIAL_BED = $(foreach set,$(CN_ENDOMETRIAL_SUBTYPES), $(HOME)/share/reference/annotation_gene_lists/cn_reg/endometrial/copy_number_regulated_genes_subtype_$(set)_spearmanrsquare0.4_fdrbh_adjp_lt0.05.HUGO.bed)
CN_BREAST_BED = $(foreach set,$(CN_BREAST_SUBTYPES), $(HOME)/share/reference/annotation_gene_lists/cn_reg/breast/metabric_subtype_$(set)_copy_number_regulated_genes_std0.5_spearmanrsquare0.4_fdrbh_adjp_lt0.05.HUGO.bed)
%.cn_reg.vcf : %.vcf
	$(call LSCRIPT_MEM,8G,12G,"$(ADD_GENE_LIST_ANNOTATION) --genome $(REF) --geneBed $(subst $(space),$(,),$(strip $(CN_ENDOMETRIAL_BED)) $(strip $(CN_BREAST_BED))) --name $(subst $(space),$(,),$(foreach set,$(strip $(CN_ENDOMETRIAL_SUBTYPES)),endometrial_$(set)) $(foreach set,$(strip $(CN_BREAST_SUBTYPES)),breast_$(set))) --outFile $@ $< && $(RM) $< $<.idx")

#-cancer does nothing
#%.som_eff.vcf : %.vcf %.vcf.pair
#	$(call CHECK_VCF,$<,$@,$(call LSCRIPT_CHECK_MEM,9G,14G,"$(call SNP_EFF_MEM,8G) ann -cancer -cancerSamples $(<<) $(SNP_EFF_OPTS) $(SNP_EFF_GENOME) $< > $@ && $(RM) $^"))

%.$(ANNOVAR_REF)_multianno.vcf : %.vcf
	$(call LSCRIPT_CHECK_MEM,7G,9G,"$(ANNOVAR) -out $* $(ANNOVAR_OPTS) $< $(ANNOVAR_DB) && $(RM) $< $<.idx")

%.norm.vcf.gz : %.vcf
	$(call LSCRIPT_MEM,9G,12G,"sed '/^##GATKCommandLine/d;/^##MuTect/d;' $< | \
		$(VT) view -h -f PASS - | \
		$(VT) decompose -s - | \
		$(VT) normalize -r $(REF_FASTA) - | \
		$(call SNP_EFF_MEM,8G) ann -c $(SNP_EFF_CONFIG) $(SNP_EFF_GENOME) -formatEff -classic | \
		bgzip -c > $@")

%.vcf.gz.tbi : %.vcf.gz
	$(call LSCRIPT_MEM,3G,5G,"$(VT) index $<")

endif
VCFTOOLS_MK = true

