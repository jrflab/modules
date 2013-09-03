# vim: set ft=make :
# sub module containing vcf related tools

include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

# flags for non-gatk snp eff
SNP_EFF_FLAGS ?= -canon -ud 0 -no-intron -no-intergenic -no-utr
DEPTH_FILTER ?= 5


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
	$(call LSCRIPT_MEM,9G,12G,"$(call SNP_SIFT_MEM,8G) annotate $(DBSNP) $< > $@ 2> $(LOG)")

# apply sample depth filter
%.dp_ft.vcf : %.vcf
	$(call LSCRIPT_MEM,2G,5G,"cat $< | $(call SNP_SIFT_MEM,2G) filter '(exists GEN[*].AD) & (GEN[*].AD[1] > $(DEPTH_FILTER))' > $@ 2> $(LOG)")

# add exon distance
%.exondist.vcf : %.vcf
	$(call INIT_MEM,2G,3G) $(INTRON_POSN_LOOKUP) $< > $@ 2> $(LOGDIR)/$@.log

# extract vcf to table
tables/%.nsfp.ann.opl_eff.txt : vcf/%.nsfp.ann.eff.vcf
	$(call LSCRIPT_MEM,2G,5G,"S1=`grep '^#CHROM' $< | cut -f 10`; \
		S2=`grep '^#CHROM' $< | cut -f 11`; \
		S3=`grep '^#CHROM' $< | cut -f 12`; \
		S4=`grep '^#CHROM' $< | cut -f 13`; \
	$(VCF_EFF_ONE_PER_LINE) < $< | $(call SNP_SIFT_MEM,2G) extractFields - $(ALL_VCF_EFF_FIELDS) | sed \"1s/GEN\[0\]/\$$S1/g; 1s/GEN\[1\]/\$$S2/g; 1s/GEN\[2\]/\$$S3/g; 1s/GEN\[3\]/\$$S4/g \" > $@")

%.eff.txt : %.opl_eff.txt
	$(INIT) $(PERL) $(VCF_JOIN_EFF) < $< > $@ 2> $(LOG)
	
tables/%.nsfp.annotated.txt : vcf/%.nsfp.annotated.vcf
	$(call LSCRIPT_MEM,2G,5G,"$(call SNP_SIFT_MEM,2G) extractFields $< $(ALL_VCF_GATK_FIELDS) > $@")

%.pass.txt : %.txt
	$(INIT) head -1 $< > $@ && grep -e "PASS" $< >> $@ || true

%.novel.txt : %.txt
	$(INIT) id=`head -1 $< | tr '\t' '\n' | grep -n "^ID$$" | sed 's/:.*//'`; \
	gmaf=`head -1 $< | tr '\t' '\n' | grep -n "^GMAF$$" | sed 's/:.*//'`; \
	awk -v id=$$id -v gmaf=$$gmaf 'NR == 1 || length($$id) == 1 || $$gmaf < 0.01 { print }' $< > $@ || true


# merge tables
tables/all.%.txt : $(foreach sample,$(SAMPLES),tables/$(sample).%.txt)
	$(call INIT_MEM,2G,3G) $(RBIND) --sampleName $< $^ > $@

ifdef TUMOR_SAMPLES
tables/allTN.%.txt : $(foreach tumor,$(TUMOR_SAMPLES),tables/$(tumor)_$(normal_lookup.$(tumor)).%.txt)
	$(INIT) $(RSCRIPT) $(RBIND) --tumorNormal $^ > $@

define annotate-tumor-normal
vcf/$1_$2.%.ann.vcf : vcf/$1_$2.%.vcf bam/$1.bam bam/$2.bam bam/$1.bai bam/$2.bai
	$$(call LSCRIPT_PARALLEL_MEM,4,2G,3G,"$$(call GATK_MEM,8G) -T VariantAnnotator -nt 4 -R $$(REF_FASTA) $$(foreach ann,$$(VCF_ANNOTATIONS),-A $$(ann) ) --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -o $$@")
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call annotate-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))
endif

define annotate-sample
vcf/$1.%.ann.vcf : vcf/$1.%.vcf bam/$1.bam bam/$1.bai
	$$(call LSCRIPT_PARALLEL_MEM,4,2G,3G,"$$(call GATK_MEM,8G) -T VariantAnnotator -nt 4 -R $$(REF_FASTA) $$(foreach ann,$$(VCF_ANNOTATIONS),-A $$(ann) ) --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -o $$@")
endef
$(foreach sample,$(SAMPLES),$(eval $(call annotate-sample,$(sample))))

tables/all.%.txt : $(foreach sample,$(SAMPLES),tables/$(sample).%.txt)
	$(call INIT_MEM,2G,3G) $(RBIND) --sampleName $< $^ > $@

# VariantEval: generate vcf report
reports/%.grp : $(foreach sample,$(SAMPLES),vcf/$(sample).%.vcf) $(foreach sample,$(SAMPLES),vcf/$(sample).%.vcf.idx)
	$(call LSCRIPT_MEM,2G,5G,"$(call GATK_MEM,2G) -T VariantEval $(foreach sm,$(REPORT_STRATIFICATION), --stratificationModule $(sm)) -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach eval,$(filter %.vcf,$^), --eval:$(call strip-suffix,$(notdir $(eval))) $(eval)) -o $@ &> $(LOG)")

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

CHASM = ssh unagi 


#%.txt : %.vcf
#	$(call INIT_MEM,2G,3G) $(VCF_TO_TABLE) $< > $@
#
