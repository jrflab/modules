# mutation signature report
include modules/Makefile.inc
include modules/variant_callers/gatk.inc
include modules/variant_callers/variantCaller.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR = log/mutsig_report.$(NOW)

VCF2VRANGES = $(RSCRIPT) modules/mut_sigs/vcf2VRanges.R
KNIT = $(RSCRIPT) modules/scripts/knit.R
ALEXANDROV_DATA = $(HOME)/share/reference/sanger_30_mutsig_prob.txt
MUTSIG_REPORT = modules/mut_sigs/mutSigReport2.Rmd
MUTSIG_REPORT_OPTS = --name $(PROJECT_NAME) \
					 --alexandrovData $(ALEXANDROV_DATA) \
					 $(if $(TARGETS_FILE),--targetBed $(TARGETS_FILE))


.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: mutect_mutsig_reports

mutect_mutsig_reports : mutsig_report/mutect/mutsig_report.timestamp

mutsig_report/mutect/mutsig_report.timestamp : $(foreach pair,$(SAMPLE_PAIRS),mutsig_report/vrange/$(pair).$(call SOMATIC_VCF_SUFFIXES,mutect).VRanges.Rdata)
	$(call LSCRIPT_NAMED_PARALLEL_MEM,mutect_mutsig_report,4,3G,5G,"$(KNIT) $(MUTSIG_REPORT) $(@D) --ncores 4 --outDir $(@D) $(MUTSIG_REPORT_OPTS) $^ && touch $@")

mutsig_report/vrange/%.VRanges.Rdata : vcf/%.vcf
	$(call LSCRIPT_MEM,7G,10G,"$(VCF2VRANGES) --genome $(REF) --outFile $@ $<")

#MAKE_TRINUC_MAF = $(PYTHON) $(HOME)/share/usr/mutation-signatures/make_trinuc_maf.py
#MUT_SIG_MAIN = $(PYTHON) $(HOME)/share/usr/mutation-signatures/main.py
#MUT_SIGS = $(HOME)/share/usr/mutation-signatures/Stratton_signatures30.txt
#mutsig_report2/mutect_mutsig_report.timestamp : allmaf/allTN.mutect.$(call SOMATIC_FILTER_SUFFIX,mutect).pass.trinuc.maf
	#$(call LSCRIPT_NAMED_MEM,$*_mutect_mutsig_report2,8G,12G,"$(MUT_SIG_MAIN) $(MUT_SIGS) $< $(@D)/mutect && touch $@")

include modules/vcf_tools/vcftools.mk
