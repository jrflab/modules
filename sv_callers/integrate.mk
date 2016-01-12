# integrate on wgseq, rnaseq
# pre-req: tophat

##### MAKE INCLUDES #####
include modules/Makefile.inc

LOGDIR = log/integrate_rnaseq.$(NOW)

..DUMMY := $(shell mkdir -p version; echo "$(INTEGRATE) &> version/integrate.txt")

INTEGRATE_ONCOFUSE = $(RSCRIPT) modules/sv_callers/integrateOncofuse.R
INTEGRATE_ONCOFUSE_OPTS = --oncofuseJar $(ONCOFUSE_JAR) --oncofuseTissueType $(ONCOFUSE_TISSUE_TYPE) --java $(JAVA_BIN) 
ONCOFUSE_TISSUE_TYPE ?= EPI

WGSEQ_DIR = ../wgseq

RNA_TUMOR_NORMAL_FILE ?= rnaseq_tumor_normal.txt
SAMPLES = $(shell cut -f1 $(RNA_TUMOR_NORMAL_FILE))
T = $(shell cut -f2 $(RNA_TUMOR_NORMAL_FILE))
N = $(shell cut -f3 $(RNA_TUMOR_NORMAL_FILE))
$(foreach i,$(shell seq 1 $(words $(SAMPLES))),$(eval tumor.$(word $i,$(SAMPLES)) = $(word $i,$(T))))
$(foreach i,$(shell seq 1 $(words $(SAMPLES))),$(eval normal.$(word $i,$(SAMPLES)) = $(word $i,$(N))))

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: integrate

integrate : $(foreach sample,$(SAMPLES),integrate/oncofuse/$(sample).oncofuse.txt)

define integrate-rna-tumor-normal
integrate/sum/$1.sum.tsv : bam/$1.bam $$(WGSEQ_DIR)/bam/$2.bam $$(WGSEQ_DIR)/bam/$3.bam bam/$1.bam.bai $$(WGSEQ_DIR)/bam/$2.bam.bai $$(WGSEQ_DIR)/bam/$3.bam.bai
	$$(call LSCRIPT_MEM,8G,80G,"mkdir -p integrate/reads integrate/sum integrate/exons integrate/breakpoints; $$(INTEGRATE) fusion -reads integrate/reads/$1.reads.txt -sum integrate/sum/$1.sum.tsv -ex integrate/exons/$1.exons.tsv -bk integrate/breakpoints/$1.breakpoints.tsv $$(REF_FASTA) $$(INTEGRATE_ANN) $$(INTEGRATE_BWTS) $$(<) $$(<) $$(<<) $$(<<<)")
endef
$(foreach sample,$(SAMPLES),$(eval $(call integrate-rna-tumor-normal,$(sample),$(tumor.$(sample)),$(normal.$(sample)))))

integrate/oncofuse/%.oncofuse.txt : integrate/sum/%.sum.tsv
	$(call LSCRIPT_MEM,7G,10G,"$(INTEGRATE_ONCOFUSE) $(INTEGRATE_ONCOFUSE_OPTS) \
		--sumFile $< \
		--exonsFile integrate/exons/$*.exons.tsv \
		--breakpointsFile integrate/breakpoints/$*.breakpoints.tsv \
		--outPrefix $(@D)/$*")



