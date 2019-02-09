include modules/Makefile.inc
include modules/bam_tools/processBam.mk

LOGDIR = log/integrate_rnaseq.$(NOW)
.PHONY: integrate_rnaseq integrate_rnaseq/oncofuse integrate_rnaseq/reads integrate_rnaseq/sum integrate_rnaseq/exons integrate_rnaseq/breakpoints integrate_rnaseq/usv

INTEGRATE_MINW ?= 2.0
INTEGRATE_LARGENUM ?= 4
INTEGRATE_OPTS = -minW $(INTEGRATE_MINW) -largeNum $(INTEGRATE_LARGENUM)
INTEGRATE_ONCOFUSE = $(RSCRIPT) modules/sv_callers/integrateOncofuse.R
INTEGRATE_ONCOFUSE_OPTS = --oncofuseJar $(ONCOFUSE_JAR) --oncofuseTissueType $(ONCOFUSE_TISSUE_TYPE) --java $(JAVA_BIN) \
						  --mysqlHost $(EMBL_MYSQLDB_HOST) --mysqlPort $(EMBL_MYSQLDB_PORT) --mysqlUser $(EMBL_MYSQLDB_USER) \
						  $(if $(EMBL_MYSQLDB_PW),--mysqlPassword $(EMBL_MYSQLDB_PW)) --mysqlDb $(EMBL_MYSQLDB_DB)
ONCOFUSE_TISSUE_TYPE ?= EPI
INTEGRATE_TO_USV = python modules/sv_callers/integrate2usv.py


integrate_rnaseq: integrate_rnaseq/all.integrate.oncofuse.txt #$(foreach sample,$(TUMOR_SAMPLES),integrate_rnaseq/usv/$(sample).integrate_rnaseq.tsv)

define init-integrate
integrate_rnaseq/reads/%.reads.txt integrate_rnaseq/sum/%.sum.tsv integrate_rnaseq/exons/%.exons.tsv integrate_rnaseq/breakpoints/%.breakpoints.tsv : bam/%.bam bam/%.bam.bai
	$$(call RUN,-s 8G -m 40G,"mkdir -p integrate_rnaseq/reads integrate_rnaseq/sum integrate_rnaseq/exons integrate_rnaseq/breakpoints && \
							  $$(INTEGRATE) fusion $$(INTEGRATE_OPTS) -reads integrate_rnaseq/reads/$$(*).reads.txt -sum integrate_rnaseq/sum/$$(*).sum.tsv -ex integrate_rnaseq/exons/$$(*).exons.tsv -bk integrate_rnaseq/breakpoints/$$(*).breakpoints.tsv $$(REF_FASTA) $$(INTEGRATE_ANN) $$(INTEGRATE_BWTS) $$(<) $$(<)")

endef
$(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call init-integrate,$(sample))))

define init-oncofuse
integrate_rnaseq/oncofuse/%.oncofuse.txt : integrate_rnaseq/sum/%.sum.tsv integrate_rnaseq/exons/%.exons.tsv integrate_rnaseq/breakpoints/%.breakpoints.tsv
	$$(call RUN,-s 7G -m 10G,"$$(INTEGRATE_ONCOFUSE) $$(INTEGRATE_ONCOFUSE_OPTS) \
							  --ref $$(REF) \
					  		  --sumFile $$(<) \
							  --exonsFile $$(<<) \
							  --breakpointsFile $$(<<<) \
							  --outPrefix $$(@D)/$$(*)")
		
endef
$(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call init-oncofuse,$(sample))))

define integrate-usv
integrate_rnaseq/usv/%.integrate_rnaseq.tsv : integrate_rnaseq/breakpoints/%.breakpoints.tsv integrate_rnaseq/sum/%.sum.tsv integrate_rnaseq/exons/%.exons.tsv
	$$(call RUN, -s 8G -m 24G,"$$(INTEGRATE_TO_USV) --breakpoints_file $$(<) --sum_file $$(<<) --exons_file $$(<<<) > $$(@)")
	
endef
$(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call integrate-usv,$(sample))))


integrate_rnaseq/all.integrate.oncofuse.txt : $(foreach sample,$(TUMOR_SAMPLES),integrate_rnaseq/oncofuse/$(sample).oncofuse.txt)
	$(INIT) (head -1 $< | sed 's/^/sample\t/'; for x in $^; do sed "1d;s/^/$$(basename $${x%%.oncofuse.txt})\t/" $$x; done) > $@


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
