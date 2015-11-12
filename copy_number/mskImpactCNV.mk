# MSK Impact procedure
include modules/Makefile.inc

LOGDIR = log/msk_impact.$(NOW)


.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : msk_impact

MODIFIED_TARGETS_FILE = $(TARGETS_FILE:.bed=.w100.bed)
$(MODIFIED_TARGETS_FILE) : $(TARGETS_FILE)
	$(INIT) $(BEDTOOLS) makewindows -w 100 -s 101 -b $< > $@


NUC_FILE = $(MODIFIED_TARGETS_FILE:.bed=.nuc.bed)
$(NUC_FILE) : $(MODIFIED_TARGETS_FILE)
	$(INIT) $(BEDTOOLS) nuc -fi $(REF_FASTA) -bed $< > $@

gatk_cov/%.doc : bam/%.bam $(NUC_FILE)
	$(call LSCRIPT_MEM,8G,10G,$(call GATK_MEM,7G) -T DepthOfCoverage -R $(REF_FASTA) -I $< -L $(<<) -O $@)

