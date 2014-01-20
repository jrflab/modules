# vim: set ft=make :
# sub module containing vcf related tools

#include ~/share/modules/Makefile.inc
#include ~/share/modules/gatk.inc

# flags for non-gatk snp eff
SNP_EFF_FLAGS ?= -canon -ud 0 -no-intron -no-intergenic -no-utr
DEPTH_FILTER ?= 5

CHASM = $(RSCRIPT) $(HOME)/share/scripts/chasmVcf.R 
#CHASM_DIR = /ifs/opt/common/CHASM/CHASMDL.1.0.7
CHASM_DIR = $(HOME)/share/usr/CHASM
CHASM_PYTHON = $(HOME)/share/usr/bin/python
CHASM_CLASSIFIER = Breast

FATHMM = $(RSCRIPT) $(HOME)/share/scripts/fathmmVcf.R 
FATHMM_DIR = $(HOME)/share/usr/fathmm
FATHMM_PYTHON = $(HOME)/share/usr/bin/python
FATHMM_PYTHONPATH = $(HOME)/share/usr/lib/python:$(HOME)/share/usr/lib/python2.7


ifdef NORMAL_VCF
%.nft.vcf : %.vcf
	$(call LSCRIPT_MEM,2G,2G,"$(NORMAL_FILTER) $< $(NORMAL_VCF) > $@")
endif

# run snp eff
%.eff.vcf : %.vcf %.vcf.idx
	$(call LSCRIPT_MEM,9G,14G,"$(call SNP_EFF_MEM,8G) -i vcf -o vcf $(SNP_EFF_FLAGS) $(SNP_EFF_GENOME) $< > $@  2> $(LOGDIR)/$(@F).log")

# run snp sift to annotated with dbnsfp
%.nsfp.vcf : %.vcf %.vcf.idx
	$(call LSCRIPT_MEM,9G,12G,"$(call SNP_SIFT_MEM,8G) dbnsfp -f $(subst $( ),$(,),$(NSFP_FIELDS)) -v $(DB_NSFP) $< > $@ 2> $(LOG)")

# run gatk snp eff
%.gatk_eff.vcf : %.vcf %.vcf.idx
	$(call LSCRIPT_MEM,5G,8G,"$(call SNP_EFF_MEM,4G) -i vcf -o gatk $(SNP_EFF_GENOME) $< > $@  2> $(LOGDIR)/$(@F).log")

# process snp eff output with gatk %=sample.indels/snps
%.annotated.vcf : %.vcf %.gatk_eff.vcf %.gatk_eff.vcf.idx %.vcf.idx 
	$(call LSCRIPT_PARALLEL_MEM,5,2G,3G,"$(call GATK_MEM,8G) -T VariantAnnotator \
	-R $(REF_FASTA) -nt 5 -A SnpEff  --variant $<  --snpEffFile $(word 2,$^) -o $@ &> $(LOGDIR)/$@.log")

%.dbsnp.vcf : %.vcf %.vcf.idx 
	$(call LSCRIPT_MEM,9G,12G,"$(call CHECK_VCF,$<,$@,$(call SNP_SIFT_MEM,8G) annotate $(DBSNP) $< > $@ 2> $(LOG))")
#$(call LSCRIPT_MEM,9G,12G,"$(call SNP_SIFT_MEM,8G) annotate $(DBSNP) $< > $@")

# apply sample depth filter
%.dp_ft.vcf : %.vcf
	$(call LSCRIPT_MEM,2G,5G,"$(call SNP_SIFT_MEM,2G) filter -f $< -p -a AD -r PASS -i AllelicDepth '(exists GEN[*].AD) & (GEN[*].AD[1] > $(DEPTH_FILTER))' > $@")

# apply dp filter for somatic sniper
%.ss_dp_ft.vcf : %.vcf
	$(call LSCRIPT_MEM,2G,5G,"$(call SNP_SIFT_MEM,2G) filter -p -a DP -i Depth -r PASS -f $< '(exists GEN[ALL].DP) & (GEN[ALL].DP >= $(DEPTH_FILTER))' > $@")

%.ss_ft.vcf : %.vcf
	$(call LSCRIPT_MEM,2G,5G,"$(call SNP_SIFT_MEM,2G) filter -p -a SS -i 'non-ref normal' -r PASS -f $< '(exists GEN[ALL].SS) & (GEN[0].SS = 0)' > $@")


# varscan TN variant allele frequency: min tumor freq > 5% ; max normal freq < 5%
%.freq_ft.vcf : %.vcf
	$(call LSCRIPT_MEM,2G,5G,"sed '/##FORMAT=<ID=FREQ/ s/String/Float/; /^#/! s/%//g' $< | $(call SNP_SIFT_MEM,2G) filter '(exists GEN[*].FREQ) & (GEN[0].FREQ < 5) & (GEN[1].FREQ[0] > 5)' > $@")

# varscan depth filter (b/c varscan is dumb and only gives variant depth)
%.vdp_ft.vcf : %.vcf
	$(call LSCRIPT_MEM,2G,5G,"cat $< | $(call SNP_SIFT_MEM,2G) filter '(exists GEN[*].AD) & (GEN[*].AD > $(DEPTH_FILTER))' > $@")

# add exon distance
%.exondist.vcf : %.vcf
	$(call LSCRIPT_MEM,2G,3G,"$(INTRON_POSN_LOOKUP) $< > $@")

# extract vcf to table
tables/%.opl_tab.txt : vcf/%.vcf
	$(call LSCRIPT_MEM,2G,5G,"NS=$(call COUNT_SAMPLES,$*); \
	$(VCF_EFF_ONE_PER_LINE) < $< | $(call SNP_SIFT_MEM,2G) extractFields - $(call SNP_VCF_EFF_FIELDS,$(call COUNT_SAMPLES,$*)) > $@; \
	for i in \`seq 0 \$$((\$$NS - 1))\`; do \
	S=\`grep '^#CHROM' $< | cut -f \$$((\$$i + 10))\`; \
	sed -i \"1s/GEN\[\$$i\]/\$$S/g;\" $@; \
	done")


%.tab.txt : %.opl_tab.txt
	$(INIT) $(PERL) $(VCF_JOIN_EFF) < $< > $@ 2> $(LOG)
	
tables/%.nsfp.annotated.txt : vcf/%.nsfp.annotated.vcf
	$(call LSCRIPT_MEM,2G,5G,"$(call SNP_SIFT_MEM,2G) extractFields $< $(ALL_VCF_GATK_FIELDS) > $@")

%.pass.txt : %.txt
	$(INIT) head -1 $< > $@ && awk '$$6 == "PASS" { print }' $< >> $@ || true

%.novel.txt : %.txt
	$(INIT) id=`head -1 $< | tr '\t' '\n' | grep -n "^ID$$" | sed 's/:.*//'`; \
	gmaf=`head -1 $< | tr '\t' '\n' | grep -n "^GMAF$$" | sed 's/:.*//'`; \
	awk -v id=$$id -v gmaf=$$gmaf 'NR == 1 || length($$id) == 1 || $$gmaf < 0.01 { print }' $< > $@ || true

# merge tables
tables/all.%.txt : $(foreach sample,$(SAMPLES),tables/$(sample).%.txt)
	$(INIT) $(RBIND) --sampleName $< $^ > $@

ifdef SAMPLE_SETS
define somatic-filter-vcf-set
vcf/$$(subst $$( ),_,$1).%.som_ft.vcf : vcf/$$(subst $$( ),_,$1).%.vcf
	$$(INIT) $$(SOMATIC_FILTER_VCF) -n $$(word $$(words $1),$1) -f 0.03 $$< > $$@ 2> $$(LOG)
endef
$(foreach i,$(shell seq 1 $(NUM_SETS)),$(eval $(call somatic-filter-vcf-set,$(set.$i))))

tables/allSS.%.txt : $(foreach set,$(SAMPLE_SETS),tables/$(set).%.txt)
	$(INIT) $(RSCRIPT) $(RBIND) --normalLast $^ > $@
endif

ifdef SAMPLE_PAIRS
#define somatic-filter-vcf
#vcf/$1_$2.%.som_ft.vcf : vcf/$1_$2.%.vcf
#$$(INIT) $$(SOMATIC_FILTER_VCF)  -n $2 -f 0.03 $$< > $$@ 2> $$(LOG)
#endef
#$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call somatic-filter-vcf,$(tumor),$(normal_lookup.$(tumor)))))

tables/allTN.%.txt : $(foreach pair,$(SAMPLE_PAIRS),tables/$(pair).%.txt)
	$(INIT) $(RSCRIPT) $(RBIND) --tumorNormal $^ > $@

define som-ad-ft-tumor-normal
vcf/$1_$2.%.som_ad_ft.vcf : vcf/$1_$2.%.vcf
	$$(call LSCRIPT_MEM,8G,12G,"$$(call GATK_MEM,8G) -T VariantFiltration -R $$(REF_FASTA) -V $$< -o $$@ --filterExpression 'if (DP > 20) { vc.getGenotype(\"$2\").getAD().1 > vc.getGenotype(\"$1\").getAD.1 / 5 } else { vc.getGenotype(\"$2\").getAD.1 > 1 }' --filterName somaticAlleleDepth")
endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call som-ad-ft-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))


define ad-tumor-normal
vcf/$1_$2.%.ad.vcf : vcf/$1_$2.%.vcf bam/$1.bam bam/$2.bam bam/$1.bai bam/$2.bai
	$$(call LSCRIPT_PARALLEL_MEM,4,2G,3G,"$$(call GATK_MEM,8G) -T VariantAnnotator -nt 4 -R $$(REF_FASTA) -A DepthPerAlleleBySample --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -o $$@ -L $$<")
endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call ad-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))

define annotate-tumor-normal
vcf/$1_$2.%.ann.vcf : vcf/$1_$2.%.vcf bam/$1.bam bam/$2.bam bam/$1.bai bam/$2.bai
	$$(call LSCRIPT_PARALLEL_MEM,4,2G,3G,"$$(call GATK_MEM,8G) -T VariantAnnotator -nt 4 -R $$(REF_FASTA) $$(foreach ann,$$(VCF_ANNOTATIONS),-A $$(ann) ) --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -o $$@ -L $$<")
endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call annotate-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))

define hrun-tumor-normal
vcf/$1_$2.%.hrun.vcf : vcf/$1_$2.%.vcf bam/$1.bam bam/$2.bam bam/$1.bai bam/$2.bai
	$$(call LSCRIPT_PARALLEL_MEM,4,2G,3G,"$$(call GATK_MEM,8G) -T VariantAnnotator -nt 4 -R $$(REF_FASTA) -A HomopolymerRun --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -L $$< -o $$@")
endef
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call hrun-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))
endif

define annotate-sample
vcf/$1.%.ann.vcf : vcf/$1.%.vcf bam/$1.bam bam/$1.bai
	$$(call LSCRIPT_PARALLEL_MEM,4,2G,3G,"$$(call GATK_MEM,8G) -T VariantAnnotator -nt 4 -R $$(REF_FASTA) $$(foreach ann,$$(VCF_ANNOTATIONS),-A $$(ann) ) --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -L $$< -V $$< -o $$@")
endef
$(foreach sample,$(SAMPLES),$(eval $(call annotate-sample,$(sample))))

define hrun-sample
vcf/$1.%.hrun.vcf : vcf/$1.%.vcf bam/$1.bam bam/$1.bai
	$$(call LSCRIPT_PARALLEL_MEM,4,2G,3G,"$$(call GATK_MEM,8G) -T VariantAnnotator -nt 4 -R $$(REF_FASTA) -L $$< -A HomopolymerRun --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -o $$@")
endef
$(foreach sample,$(SAMPLES),$(eval $(call hrun-sample-chr,$(sample))))

tables/all.%.txt : $(foreach sample,$(SAMPLES),tables/$(sample).%.txt)
	$(call INIT_MEM,2G,3G) $(RBIND) --sampleName $< $^ > $@

# VariantEval: generate vcf report
reports/%.grp : $(foreach sample,$(SAMPLES),vcf/$(sample).%.vcf) $(foreach sample,$(SAMPLES),vcf/$(sample).%.vcf.idx)
	$(call LSCRIPT_MEM,2G,5G,"$(call GATK_MEM,2G) -T VariantEval $(foreach sm,$(REPORT_STRATIFICATION), --stratificationModule $(sm)) -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach eval,$(filter %.vcf,$^), --eval:$(call strip-suffix,$(notdir $(eval))) $(eval)) -o $@ &> $(LOG)")
ifdef SAMPLE_PAIRS
reports/%.grp : $(foreach pair,$(SAMPLE_PAIRS),vcf/$(pair).%.vcf vcf/$(pair).%.vcf.idx)
	$(call LSCRIPT_MEM,2G,5G,"$(call GATK_MEM,2G) -T VariantEval $(foreach sm,$(REPORT_STRATIFICATION), --stratificationModule $(sm)) -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach eval,$(filter %.vcf,$^), --eval:$(call strip-suffix,$(notdir $(eval))) $(eval)) -o $@ &> $(LOG)")
endif

NON_SILENT_EFF = START_GAINED SPLICE_SITE_ACCEPTOR SPLICE_SITE_DONOR START_LOST NON_SYNONYMOUS_CODING FRAME_SHIFT CODON_CHANGE CODON_INSERTION CODON_CHANGE_PLUS_CODON_INSERTION CODON_DELETION CODON_CHANGE_PLUS_CODON_DELETION STOP_GAINED STOP_LOST NON_SYNONYMOUS_START
%.nonsilent.txt : %.txt
	$(INIT) head -1 $< > $@ && sed '1d' $< | grep $(foreach eff,$(NON_SILENT_EFF), -e $(eff)) >> $@ || true

NON_SILENT_CODING_EFF = START_GAINED START_LOST NON_SYNONYMOUS_CODING FRAME_SHIFT CODON_CHANGE CODON_INSERTION CODON_CHANGE_PLUS_CODON_INSERTION CODON_DELETION CODON_CHANGE_PLUS_CODON_DELETION STOP_GAINED STOP_LOST NON_SYNONYMOUS_START
%.nonsilent_cds.txt : %.txt
	$(INIT) head -1 $< > $@ && sed '1d' $< | grep $(foreach eff,$(NON_SILENT_CODING_EFF), -e $(eff)) >> $@ || true

%.missense.txt : %.txt
	$(INIT) head -1 $< > $@ && sed '1d' $< | grep -e MISSENSE >> $@ || true

%.silent.txt : %.txt
	$(INIT) head -1 $< > $@ && sed '1d' $< | grep -e SILENT >> $@ || true

%.vcf.idx : %.vcf
	$(call LSCRIPT_MEM,4G,8G,"$(IGVTOOLS) index $< &> $(LOG)")

%.chasm.vcf : %.vcf
	$(call LSCRIPT_MEM,8G,17G,"$(CHASM) --genome $(REF) --chasmDir $(CHASM_DIR) --python $(CHASM_PYTHON) --outFile $@ $<")

%.fathmm.vcf : %.vcf %.fathmmInput.Rdata
	$(call LSCRIPT_MEM,12G,22G,"PYTHONPATH=$(FATHMM_PYTHONPATH) $(FATHMM) --genome $(REF) --fathmmDir $(FATHMM_DIR) --outFile $@ --python $(FATHMM_PYTHON) $^")

PRED_CODING = $(RSCRIPT) $(HOME)/share/scripts/vcfPredictCoding.R
%.predCoding.Rdata : %.vcf
	$(call LSCRIPT_MEM,12G,15G,"$(PRED_CODING) --outFile $@ --genome $(REF) --ref $(REF_FASTA) --ensemblTxdb $(ENSEMBL_TXDB) $<")

FATHMM_INPUT = $(RSCRIPT) $(HOME)/share/scripts/fathmmVcfInput.R
%.fathmmInput.Rdata : %.predCoding.Rdata
	$(INIT) $(FATHMM_INPUT) --outFile $@ --genome $(REF) --outFile $@ --ensemblTxdb $(ENSEMBL_TXDB) $< &> $(LOG)

MUT_ASS = $(RSCRIPT) $(HOME)/share/scripts/mutAssVcf.R
%.mutass.vcf : %.vcf
	$(call LSCRIPT_MEM,12G,15G,$(MUT_ASS) --outFile $@ --maData $(MUT_ASS_RDATA) $<)

TRANSFIC = $(RSCRIPT) $(HOME)/share/scripts/transficVcf.R
TRANSFIC_PERL_SCRIPT = $(HOME)/share/usr/transfic/bin/transf_scores.pl
%.transfic.vcf : %.vcf
	$(call LSCRIPT_MEM,3G,4G,"$(TRANSFIC) --genome $(REF) --transfic $(TRANSFIC_PERL_SCRIPT) --outFile $@ $<")

