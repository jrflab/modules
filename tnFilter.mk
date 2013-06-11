# Naive tumour-normal filter
TN_FILTER = $(PERL) $(HOME)/share/scripts/normalFilterVCF.pl

# $(eval $(call tn-filter,tumor.vcf,normal.vcf))
define tn-filter
$1.%.tnFiltered.vcf : $1.%.vcf $2.%.vcf
	SGE_RREQ="$$(SGE_RREQ) $$(call MEM_FREE,4G,8G)" \
	$$(TN_FILTER) $$^ > $$@
endef
$(foreach i,$(shell seq 1 $(NSAMPLES)),$(eval $(call tn-filter,$(word $i,$(TUMOR_SAMPLES)),$(word $i,$(NORMAL_SAMPLES)))))
