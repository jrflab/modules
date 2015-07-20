# integrate only on rnaseq
# pre-req: tophat

##### MAKE INCLUDES #####
include modules/Makefile.inc

LOGDIR = log/integrate_rnaseq.$(NOW)

..DUMMY := $(shell mkdir -p version; echo "$(INTEGRATE) &> version/integrate.txt")

INTEGRATE = $(HOME)/share/usr/bin/Integrate
INTEGRATE_ONCOFUSE = $(RSCRIPT) modules/sv_callers/integrateOncofuse.R
INTEGRATE_ONCOFUSE_OPTS = --oncofuseJar $(ONCOFUSE_JAR) --oncofuseTissueType $(ONCOFUSE_TISSUE_TYPE) --java $(JAVA) 
ONCOFUSE_TISSUE_TYPE ?= EPI
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: integrate_rnaseq

integrate_rnaseq : $(foreach sample,$(SAMPLES),integrate/oncofuse/$(sample).oncofuse.txt)

integrate/reads/%.reads.txt integrate/sum/%.sum.tsv integrate/exons/%.exons.tsv integrate/breakpoints/%.breakpoints.tsv : bam/%.bam.md5 bam/%.bam.bai
	$(call LSCRIPT_MEM,8G,80G,"mkdir -p integrate/reads integrate/sum integrate/exons integrate/breakpoints; $(INTEGRATE) fusion -reads integrate/reads/$*.reads.txt -sum integrate/sum/$*.sum.tsv -ex integrate/exons/$*.exons.tsv -bk integrate/breakpoints/$*.breakpoints.tsv $(REF_FASTA) $(INTEGRATE_ANN) $(INTEGRATE_BWTS) $(<M) $(<M)")

integrate/oncofuse/%.oncofuse.txt : integrate/sum/%.sum.tsv integrate/exons/%.exons.tsv integrate/breakpoints/%.breakpoints.tsv
	$(call LSCRIPT_MEM,7G,10G,"$(INTEGRATE_ONCOFUSE) $(INTEGRATE_ONCOFUSE_OPTS) \
		--sumFile $< \
		--exonsFile $(<<) \
		--breakPointsFile $(<<<) \
		--outPrefix $(@D)/$*")


include modules/bam_tools/processBam.mk
