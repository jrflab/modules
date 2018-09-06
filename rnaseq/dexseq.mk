include modules/Makefile.inc

LOGDIR ?= log/exon_counts.$(NOW)
PHONY += dex_seq

dex_seq : $(foreach sample,$(TUMOR_SAMPLES),dex_seq/$(sample).taskcomplete)

define exon-count
dex_seq/$1.txt : star/bam/$1.star.sorted.filtered.bam
	$$(call RUN,-c -s 8G -m 12G,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules.0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/dexseq && \
								 /home/${USER}/share/usr/anaconda-envs/dexseq/lib/R/library/DEXSeq/python_scripts/dexseq_count.py -f bam /home/${USER}/share/reference/Ensembl/Homo_sapiens.GRCh37.75.gff star/bam/$$<.star.sorted.filtered.bam dex_seq/$$<.txt")
dex_seq/%.taskcomplete : dex_seq/%.txt
	$$(call RUN,-c -s 1G -m 1G,"echo $$<")
endef
$(foreach sample,$(TUMOR_SAMPLES),\
		$(eval $(call exon-count,$sample)))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
