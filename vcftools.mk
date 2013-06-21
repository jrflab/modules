# vim: set ft=make :
# sub module containing vcf related tools

include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

#SNP_EFF_FLAGS ?= -ud 0 -no-intron -no-intergenic
DEPTH_FILTER ?= 10

%.vcf.idx : %.vcf
	$(call INIT_MEM,4G,8G) $(IGVTOOLS) index $< &> $(LOGDIR)/$(@F).log

# run snp eff
%.eff.vcf : %.vcf %.vcf.idx
	$(call INIT_MEM,9G,12G) $(call SNP_EFF_MEM,8G) -i vcf -o vcf $(SNP_EFF_GENOME) $< > $@  2> $(LOGDIR)/$(@F).log

# run snp sift to annotated with dbnsfp
%.nsfp.vcf : %.vcf %.vcf.idx
	$(call INIT_MEM,9G,12G) $(call SNP_SIFT_MEM,8G) dbnsfp -v $(DB_NSFP) $< > $@ 2> $(LOG)

# run gatk snp eff
%.gatk_eff.vcf : %.vcf %.vcf.idx
	$(call INIT_MEM,5G,8G) $(call SNP_EFF_MEM,4G) -i vcf -o gatk $(SNP_EFF_GENOME) $< > $@  2> $(LOGDIR)/$(@F).log

# process snp eff output with gatk %=sample.indels/snps
%.annotated.vcf : %.vcf %.gatk_eff.vcf %.gatk_eff.vcf.idx %.vcf.idx 
	$(call INIT_PARALLEL_MEM,5,2G,3G) $(call GATK_MEM,8G) -T VariantAnnotator \
	-R $(REF_FASTA) -nt 5 -A SnpEff  --variant $<  --snpEffFile $(word 2,$^) -o $@ &> $(LOGDIR)/$@.log

# apply sample depth filter
%.dp_ft.vcf : %.vcf
	$(call INIT_MEM,2G,5G) cat $< | $(call SNP_SIFT_MEM,2G) filter "GEN[ALL].DP > $(DEPTH_FILTER)" > $@ 2> $(LOG)

# add exon distance
%.exondist.vcf : %.vcf
	$(call INIT_MEM,2G,3G) $(INTRON_POSN_LOOKUP) $< > $@ 2> $(LOGDIR)/$@.log

# extract vcf to table
tables/%.eff.nsfp.txt : vcf/%.eff.nsfp.vcf
	$(call INIT_MEM,2G,5G) $(VCF_EFF_ONE_PER_LINE) < $< | $(call SNP_SIFT_MEM,2G) extractFields - $(ALL_VCF_EFF_FIELDS) > $@
	
tables/%.annotated.nsfp.txt : vcf/%.annotated.nsfp.vcf
	$(call INIT_MEM,2G,5G) $(call SNP_SIFT_MEM,2G) extractFields $< $(ALL_VCF_GATK_FIELDS) > $@


# extract vcf to table
tables/%.eff.nsfp.txt : vcf/%.eff.nsfp.vcf
	$(call INIT_MEM,2G,5G) $(VCF_EFF_ONE_PER_LINE) < $< | $(call SNP_SIFT_MEM,2G) extractFields - $(VCF_EFF_FIELDS) > $@

%.pass.txt : %.txt
	$(INIT) grep -v "REJECT" $< > $@

%.novel.txt : %.txt
	$(INIT) awk 'NR == 1 || length($$3) == 1 { print }' $< > $@

# merge tables
tables/all.%.txt : $(foreach sample,$(SAMPLES),tables/$(sample).%.txt)
	$(call INIT_MEM,2G,3G) $(RBIND) $^ > $@


# VariantEval: generate vcf report
reports/%.grp : $(foreach sample,$(SAMPLES),vcf/$(sample).%.vcf) $(foreach sample,$(SAMPLES),vcf/$(sample).%.vcf.idx)
	$(call INIT_MEM,2G,5G) $(call GATK_MEM,2G) -T VariantEval $(foreach sm,$(REPORT_STRATIFICATION), --stratificationModule $(sm)) -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach eval,$(filter %.vcf,$^), --eval:$(call strip-suffix,$(notdir $(eval))) $(eval)) -o $@ &> $(LOG)


#%.txt : %.vcf
#	$(call INIT_MEM,2G,3G) $(VCF_TO_TABLE) $< > $@
