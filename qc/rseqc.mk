LOGDIR = log/rseqc.$(NOW)
## includes
include modules/Makefile.inc

INFER_EXPT = python $(HOME)/share/usr/anaconda/bin/infer_experiment.py
INNER_DIST = python $(HOME)/share/usr/anaconda/bin/inner_distance.py
JUNC_ANNOT = python $(HOME)/share/usr/anaconda/bin/junction_annotation.py
BAM_STAT = python $(HOME)/share/usr/anaconda/bin/junction_annotation.py
CLIP_PROFILE = python $(HOME)/share/usr/anaconda/bin/clipping_profile.py
MISMATCH_PROFILE = python $(HOME)/share/usr/anaconda/bin/mismatch_profile.py
INSERTION_PROFILE = python $(HOME)/share/usr/anaconda/bin/insertion_profile.py
DELETION_PROFILE = python $(HOME)/share/usr/anaconda/bin/deletion_profile.py
GENEBODY_COV = python $(HOME)/share/usr/anaconda/bin/geneBody_coverage.py
READ_HEXAMER = python $(HOME)/share/usr/anaconda/bin/read_hexamer.py
READ_QUALITY = python $(HOME)/share/usr/anaconda/bin/read_quality.py
READ_NVC = python $(HOME)/share/usr/anaconda/bin/read_NVC.py
READ_GC = python $(HOME)/share/usr/anaconda/bin/read_GC.py
READ_DUP = python $(HOME)/share/usr/anaconda/bin/read_duplication.py
READ_DIST = python $(HOME)/share/usr/anaconda/bin/read_distribution.py
RPKM_SAT = python $(HOME)/share/usr/anaconda/bin/RPKM_saturation.py
RPKM_COUNT = python $(HOME)/share/usr/anaconda/bin/RPKM_count.py


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: rseqc

rseqc : $(foreach sample,$(SAMPLES), \
	rseqc/infer/$(sample).infer \
	rseqc/gene_body_cov/$(sample).geneBodyCoverage.txt \
	rseqc/inner_dist/$(sample).innerDistance.txt)

rseqc/infer/%.infer : bam/%.bam bam/%.bam.bai
	$(call LSCRIPT_MEM,7G,8G,"$(INFER_EXPT) -i $< -r $(REF_HOUSEKEEPING_GENE_BED) > $@")

rseqc/gene_body_cov/%.geneBodyCoverage.txt : bam/%.bam
	$(call LSCRIPT_MEM,7G,8G,"$(GENEBODY_COV) -i $< -r $(REF_HOUSEKEEPING_GENE_BED) -o rseqc/gene_body_cov/$*")

rseqc/inner_dist/%.innerDistance.txt : bam/%.bam
	$(call LSCRIPT_MEM,7G,8G,"$(INNER_DIST) -i $< -r $(REF_HOUSEKEEPING_GENE_BED) -o rseqc/inner_dist/$*")

