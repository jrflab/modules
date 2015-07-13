# MSK Impact procedure
include modules/Makefile.inc
include modules/variant_callers/gatk.inc

LOGDIR = log/norm_copynum.$(NOW)
NORMALISE_COPYNUM = $(RSCRIPT) modules/copy_number/normaliseCopyNum.R

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : norm_copynum

norm_copynum : norm_copynum/seg.txt

norm_copynum/targets.bed : $(TARGETS_FILE)
	$(INIT) $(BEDTOOLS) makewindows -w 100 -s 101 -b $< > $@

%.nuc.bed : %.bed
	$(INIT) $(BEDTOOLS) nuc -fi $(REF_FASTA) -bed $< > $@

norm_copynum/doc/%.doc : bam/%.bam norm_copynum/targets.bed
	$(call LSCRIPT_MEM,9G,15G,"$(call GATK_MEM,7G) -T DepthOfCoverage -R $(REF_FASTA) -I $< -L $(<<) -o $@")

norm_copynum/seg.txt : norm_copynum/targets.nuc.bed $(foreach sample,$(SAMPLES),norm_copynum/doc/$(sample).doc)
	$(call LSCRIPT_MEM,8G,10G,"$(NORMALISE_COPYNUM) --sampleSetsFile $(SAMPLE_SET_FILE) --nucFile $< --centromereFile $(CENTROMERE_TABLE2) --outDir $(@D) $(addsuffix .sample_interval_summary,$(filter %.doc,$^))")
