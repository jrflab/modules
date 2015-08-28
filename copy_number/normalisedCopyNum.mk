# MSK Impact procedure
include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR = log/norm_copynum.$(NOW)
NORMALISE_COPYNUM = $(RSCRIPT) modules/copy_number/normaliseCopyNum.R
NORMALISE_MIN_COV ?= 50
NORMALISE_UNDO_SD ?= 2
NORMALISE_WINDOW_SIZE ?= 100
NORMALISE_OUTLIER_SD_SCALE ?= 2.5
NORMALISE_ALPHA ?= 0.05
NORMALISE_TRIM ?= 0.05
NORMALISE_COPYNUM_OPTS = --undoSD $(NORMALISE_UNDO_SD) --outlierSDscale $(NORMALISE_OUTLIER_SD_SCALE) --alpha $(NORMALISE_ALPHA) --minCov $(NORMALISE_MIN_COV) --trim $(NORMALISE_TRIM)

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : norm_copynum

norm_copynum : norm_copynum/seg.txt

norm_copynum/targets.bed : $(TARGETS_FILE)
	$(INIT) $(BEDTOOLS) makewindows -w $(NORMALISE_WINDOW_SIZE) -s $$(($(NORMALISE_WINDOW_SIZE) + 1)) -b $< > $@

%.nuc.bed : %.bed
	$(INIT) $(BEDTOOLS) nuc -fi $(REF_FASTA) -bed $< > $@

norm_copynum/doc/%.doc : bam/%.bam norm_copynum/targets.bed
	$(call LSCRIPT_MEM,9G,15G,"$(call GATK_MEM,7G) -T DepthOfCoverage -R $(REF_FASTA) -I $< -L $(<<) -o $@")

norm_copynum/seg.txt : norm_copynum/targets.nuc.bed $(foreach sample,$(SAMPLES),norm_copynum/doc/$(sample).doc)
	$(call LSCRIPT_MEM,8G,10G,"$(NORMALISE_COPYNUM) $(NORMALISE_COPYNUM_OPTS) --sampleSetsFile $(SAMPLE_SET_FILE) --nucFile $< --centromereFile $(CENTROMERE_TABLE2) --outDir $(@D) $(addsuffix .sample_interval_summary,$(filter %.doc,$^))")
