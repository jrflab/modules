# run lumpy

LOGDIR = log/lumpy.$(NOW)

include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LUMPY_SCRIPTS_DIR = $(HOME)/share/usr/lumpy-sv/scripts
LUMPY = $(HOME)/share/usr/lumpy-sv/bin/lumpy
LUMPY_HISTO = $(PERL) $(LUMPY_SCRIPTS_DIR)/pairend_distro.pl
LUMPY_UNMAPPED_TO_FASTQ = $(PERL) $(LUMPY_SCRIPTS_DIR)/split_unmapped_to_fasta.pl
LUMPY_UNMAPPED_TO_FASTQ_OPTS = -b 20

LUMPY_OPTS = -tt 1e-3 -mw 4
LUMPY_PE_PARAMS = min_non_overlap:150$(,)discordant_z:4$(,)back_distance:20$(,)weight:1$(,)id:1$(,)min_mapping_threshold:1
LUMPY_SR_PARAMS = back_distance:20$(,)weight:1$(,)id:1$(,)min_mapping_threshold:1

BWASW_OPTS = -H

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all

all : $(foreach sample,$(SAMPLES),lumpy/bed/$(sample).sr.bedpe lumpy/bed/$(sample).pe.bedpe)

lumpy/fastq/%.um.fastq.gz : bam/%.bam
	$(call LSCRIPT,"$(SAMTOOLS) view $< | $(LUMPY_UNMAPPED_TO_FASTQ) $(LUMPY_UNMAPPED_TO_FASTQ_OPTS) | gzip -c > $@")

lumpy/sr_bam/%.sr.bam : lumpy/fastq/%.um.fastq.gz
	$(call LSCRIPT_PARALLEL_MEM,4,1G,2G,"$(BWA) bwasw $(BWASW_OPTS) -t 4 $(REF_FASTA) $< | samtools view -Sb - > $@")

lumpy/metrics/%.read_len : bam/%.bam
	$(INIT) $(MKDIR) $(@D); $(SAMTOOLS) view $< | tail -n+100000 | head -1 | awk '{ print length($$10) }' > $@

lumpy/metrics/%.histo lumpy/metrics/%.histo.txt: bam/%.bam lumpy/metrics/%.read_len
	READ_LEN=`cat $(word 2,$^)`; \
	$(call LSCRIPT,"$(SAMTOOLS) view $< | tail -n+100000 | $(LUMPY_HISTO) -rl $$READ_LEN -X 4 -N 10000 -o lumpy/metrics/$*.histo > lumpy/metrics/$*.histo.txt")

lumpy/bed/%.pe.bedpe : bam/%.bam lumpy/metrics/%.read_len lumpy/metrics/%.histo lumpy/metrics/%.histo.txt
	READ_LEN=`cat $(word 2,$^)`; \
	MEAN=`cut -f1 $(word 4,$^)`; \
	STDEV=`cut -f2 $(word 4,$^)`; \
	$(call LSCRIPT_MEM,4G,8G,"$(LUMPY) $(LUMPY_OPTS) -pe bam_file:$<$(,)histo_file:$(word 3,$^)$(,)$$MEAN$(,)$$STDEV$(,)read_length:$$READ_LEN$(,)$(LUMPY_PE_PARAMS) > $@")

lumpy/bed/%.sr.bedpe : lumpy/sr_bam/%.sr.bam lumpy/metrics/%.read_len lumpy/metrics/%.histo lumpy/metrics/%.histo.txt
	READ_LEN=`cat $(word 2,$^)`; \
	MEAN=`cut -f1 $(word 4,$^)`; \
	STDEV=`cut -f2 $(word 4,$^)`; \
	$(call LSCRIPT_MEM,4G,8G,"$(LUMPY) $(LUMPY_OPTS) -sr bam_file:$<$(,)$(LUMPY_SR_PARAMS) > $@")

ifdef SAMPLE_PAIRS
define lumpy-tumor-normal
lumpy/bed/$1_$2.pe.bedpe : bam/$1.bam bam/$2.bam lumpy/metrics/$1.read_len lumpy/metrics/$2.read_len lumpy/metrics/$1.histo lumpy/metrics/$2.histo lumpy/metrics/$1.histo.txt lumpy/metrics/$2.hist.txt
	READ_LEN1=`cat $$(word 3,$$^)`; \
	READ_LEN2=`cat $$(word 4,$$^)`; \
	MEAN1=`cut -f1 $$(word 6,$$^)`; \
	MEAN2=`cut -f1 $$(word 7,$$^)`; \
	STDEV1=`cut -f2 $$(word 6,$$^)`; \
	STDEV2=`cut -f2 $$(word 7,$$^)`; \
	$$(call LSCRIPT_MEM,4G,8G,"$$(LUMPY) $$(LUMPY_OPTS) \
	    -pe bam_file:$$<,histo_file:$$(word 5,$$^)$$(,)$$$$MEAN1$$(,)$$$$STDEV1$$(,)read_length:$$$$READ_LEN1$$(,)$(LUMPY_PE_PARAMS) \
	    -pe bam_file:$$(word 2,$$^)$$(,)histo_file:$$(word 6,$$^)$$(,)$$$$MEAN2$$(,)$$$$STDEV2$$(,)read_length:$$$$READ_LEN2$$(,)$(LUMPY_PE_PARAMS) > $@")
endef
$(foreach i,$(SETS_SEQ), \
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call lumpy-tumor-normal,$(tumor),$(call get_normal,$(set.$i)),$(chr)))))
endif

