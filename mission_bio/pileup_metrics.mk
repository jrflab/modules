include modules/Makefile.inc
include modules/config/waltz.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/pileup_metrics.$(NOW)

WALTZ_MIN_MAPQ ?= 30
WALTZ_BED_FILE ?= $(HOME)/share/lib/bed_files/MSK-ACCESS-v1_0-probe-A.sorted.bed

pileup_metrics : $(foreach sample,$(SAMPLES),waltz/$(sample)_cl_aln_srt_MD_IR_FX_BR-pileup.txt.gz) \
		 $(foreach sample,$(SAMPLES),waltz/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX-pileup.txt.gz) \
		 $(foreach sample,$(SAMPLES),waltz/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX-pileup.txt.gz) \
		 $(foreach sample,$(SAMPLES),waltz/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX-pileup.txt.gz) \
		 summary/intervals_summary.txt \
		 summary/pileup_summary_standard.txt

define waltz-genotype
waltz/$1_cl_aln_srt_MD_IR_FX_BR-pileup.txt.gz : bam/$1_cl_aln_srt_MD_IR_FX_BR.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"set -o pipefail && \
					 mkdir -p waltz && \
					 cd waltz && \
					 ln -sf ../bam/$1_cl_aln_srt_MD_IR_FX_BR.bam $1_cl_aln_srt_MD_IR_FX_BR.bam && \
					 ln -sf ../bam/$1_cl_aln_srt_MD_IR_FX_BR.bai $1_cl_aln_srt_MD_IR_FX_BR.bai && \
					 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1_cl_aln_srt_MD_IR_FX_BR.bam $$(REF_FASTA) $$(WALTZ_BED_FILE) && \
					 gzip $1_cl_aln_srt_MD_IR_FX_BR-pileup.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX_BR-pileup-without-duplicates.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX_BR-intervals.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX_BR-intervals-without-duplicates.txt && \
					 unlink $1_cl_aln_srt_MD_IR_FX_BR.bam && \
					 unlink $1_cl_aln_srt_MD_IR_FX_BR.bai && \
					 cd ..")

waltz/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX-pileup.txt.gz : bam/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"set -o pipefail && \
					 mkdir -p waltz && \
					 cd waltz && \
					 ln -sf ../bam/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bam $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bam && \
					 ln -sf ../bam/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bai $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bai && \
					 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bam $$(REF_FASTA) $$(WALTZ_BED_FILE) && \
					 gzip $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX-pileup.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX-pileup-without-duplicates.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX-intervals.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX-intervals-without-duplicates.txt && \
					 unlink $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bam && \
					 unlink $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX.bai && \
					 cd ..")

waltz/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX-pileup.txt.gz : bam/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"set -o pipefail && \
					 mkdir -p waltz && \
					 cd waltz && \
					 ln -sf ../bam/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.bam $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.bam && \
					 ln -sf ../bam/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.bai $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.bai && \
					 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.bam $$(REF_FASTA) $$(WALTZ_BED_FILE) && \
					 gzip $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX-pileup.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX-pileup-without-duplicates.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX-intervals.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX-intervals-without-duplicates.txt && \
					 unlink $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.bam && \
					 unlink $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX.bai && \
					 cd ..")

waltz/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX-pileup.txt.gz : bam/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.bam
	$$(call RUN,-c -n 4 -s 4G -m 6G,"set -o pipefail && \
					 mkdir -p waltz && \
					 cd waltz && \
					 ln -sf ../bam/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.bam $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.bam && \
					 ln -sf ../bam/$1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.bai $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.bai && \
					 $$(call WALTZ_CMD,2G,8G) org.mskcc.juber.waltz.Waltz PileupMetrics $$(WALTZ_MIN_MAPQ) $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.bam $$(REF_FASTA) $$(WALTZ_BED_FILE) && \
					 gzip $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX-pileup.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX-pileup-without-duplicates.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX-intervals.txt && \
					 gzip $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX-intervals-without-duplicates.txt && \
					 unlink $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.bam && \
					 unlink $1_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX.bai && \
					 cd ..")
					 
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call waltz-genotype,$(sample))))

summary/intervals_summary.txt : $(foreach sample,$(SAMPLES),waltz/$(sample)_cl_aln_srt_MD_IR_FX_BR-pileup.txt.gz) $(foreach sample,$(SAMPLES),waltz/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX-pileup.txt.gz) $(foreach sample,$(SAMPLES),waltz/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX-pileup.txt.gz) $(foreach sample,$(SAMPLES),waltz/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX-pileup.txt.gz)
	$(call RUN, -c -n 1 -s 64G -m 128G,"set -o pipefail && \
					    $(RSCRIPT) $(SCRIPTS_DIR)/qc/pileup_metrics.R --option 1 --sample_names '$(SAMPLES)'")
summary/pileup_summary_standard.txt : $(foreach sample,$(SAMPLES),waltz/$(sample)_cl_aln_srt_MD_IR_FX_BR-pileup.txt.gz) $(foreach sample,$(SAMPLES),waltz/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX-pileup.txt.gz) $(foreach sample,$(SAMPLES),waltz/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_SIMPLEX-pileup.txt.gz) $(foreach sample,$(SAMPLES),waltz/$(sample)_cl_aln_srt_MD_IR_FX_BR__grp_DC_MA_RG_IR_FX_DUPLEX-pileup.txt.gz)
	$(call RUN, -c -n 1 -s 64G -m 128G,"set -o pipefail && \
					    $(RSCRIPT) $(SCRIPTS_DIR)/qc/pileup_metrics.R --option 2 --sample_names '$(SAMPLES)'")

		
..DUMMY := $(shell mkdir -p version; \
	     $(JAVA8) -version &> version/pileup_metrics.txt)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: pileup_metrics
