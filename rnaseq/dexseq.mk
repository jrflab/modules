include modules/Makefile.inc

LOGDIR ?= log/exon_counts.$(NOW)
PHONY += dexseq

dexseq : $(foreach sample,$(TUMOR_SAMPLES),dexseq/$(sample).txt)

define exon-count
dexseq/%.txt : star/bam/%.star.sorted.filtered.bam
	$$(call RUN,-c -s 8G -m 12G -w 1440,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/dexseq && \
								 /home/${USER}/share/usr/anaconda-envs/dexseq/lib/R/library/DEXSeq/python_scripts/dexseq_count.py -f bam -p yes -r pos /home/${USER}/share/reference/Ensembl/Homo_sapiens.GRCh37.75.gff $$< dexseq/$$*.txt")
endef
$(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call exon-count,$sample)))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
