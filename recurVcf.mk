# intersect vcf files
##### DEFAULTS ######

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all

VCF_SUFFIXES = gatk_snps.dp_ft.dbsnp.nsfp.chasm.fathmm.eff som_sniper.ss_dp_ft.ss_ft.pass.dbsnp.nsfp.chasm.fathmm.eff mutect.som_ad_ft.pass.dbsnp.nsfp.chasm.fathmm.eff

RECUR_VCF = $(RSCRIPT) $(HOME)/share/scripts/recurVcf.R

recur_pos/%.recur.txt : $(foreach suffix,$(VCF_SUFFIXES),vcf/%.$(suffix).vcf)
	$(INIT) $(RECUR_VCF) --outFile $@ $^
