# integrate only on rnaseq
# pre-req: tophat

##### MAKE INCLUDES #####
include modules/Makefile.inc

LOGDIR = log/integrate_rnaseq.$(NOW)

..DUMMY := $(shell mkdir -p version; echo "$(INTEGRATE) &> version/integrate.txt")

INTEGRATE_MINW ?= 2.0
INTEGRATE_LARGENUM ?= 4
INTEGRATE_OPTS = -minW $(INTEGRATE_MINW) -largeNum $(INTEGRATE_LARGENUM)

INTEGRATE_ONCOFUSE = $(RSCRIPT) modules/sv_callers/integrateOncofuse.R
INTEGRATE_ONCOFUSE_OPTS = --oncofuseJar $(ONCOFUSE_JAR) --oncofuseTissueType $(ONCOFUSE_TISSUE_TYPE) --java $(JAVA_BIN) \
						  --mysqlHost $(EMBL_MYSQLDB_HOST) --mysqlPort $(EMBL_MYSQLDB_PORT) --mysqlUser $(EMBL_MYSQLDB_USER) --mysqlPassword $(EMBL_MYSQLDB_PW) --mysqlDb $(EMBL_MYSQLDB_DB)
ONCOFUSE_TISSUE_TYPE ?= EPI

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: integrate_rnaseq

integrate_rnaseq : $(foreach sample,$(SAMPLES),integrate_rnaseq/oncofuse/$(sample).oncofuse.txt)

integrate_rnaseq/reads/%.reads.txt integrate_rnaseq/sum/%.sum.tsv integrate_rnaseq/exons/%.exons.tsv integrate_rnaseq/breakpoints/%.breakpoints.tsv : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,8G,40G,"mkdir -p integrate_rnaseq/reads integrate_rnaseq/sum integrate_rnaseq/exons integrate_rnaseq/breakpoints; $(INTEGRATE) fusion $(INTEGRATE_OPTS) -reads integrate_rnaseq/reads/$*.reads.txt -sum integrate_rnaseq/sum/$*.sum.tsv -ex integrate_rnaseq/exons/$*.exons.tsv -bk integrate_rnaseq/breakpoints/$*.breakpoints.tsv $(REF_FASTA) $(INTEGRATE_ANN) $(INTEGRATE_BWTS) $(<) $(<)")

integrate_rnaseq/oncofuse/%.oncofuse.txt : integrate_rnaseq/sum/%.sum.tsv integrate_rnaseq/exons/%.exons.tsv integrate_rnaseq/breakpoints/%.breakpoints.tsv
	$(call LSCRIPT_MEM,7G,10G,"$(INTEGRATE_ONCOFUSE) $(INTEGRATE_ONCOFUSE_OPTS) \
		--ref $(REF) \
		--sumFile $< \
		--exonsFile $(<<) \
		--breakpointsFile $(<<<) \
		--outPrefix $(@D)/$*")


include modules/bam_tools/processBam.mk
