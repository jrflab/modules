include innovation-lab/Makefile.inc
include innovation-lab/config/gatk.inc
include innovation-lab/config/align.inc

LOGDIR ?= log/process_bam.$(NOW)

ifndef PROCESS_BAM_MK

MERGE_SPLIT_BAMS ?= false
BAM_CHR1_BASE_RECAL ?= false
BAM_BASE_RECAL_OPTS = -knownSites $(DBSNP) $(if $(findstring true,$(BAM_CHR1_BASE_RECAL)),-L $(word 1,$(CHROMOSOMES)))

ifneq ($(KNOWN_INDELS),)
BAM_REALN_OPTS = --knownAlleles $(KNOWN_INDELS)
BAM_REALN_TARGET_OPTS = --known $(KNOWN_INDELS)
endif

BAM_FILTER_FLAGS ?= 768

.DELETE_ON_ERROR:
.SECONDARY: 

BAM_REPROCESS ?= false

BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
ifeq ($(BAM_REPROCESS),true)
processed_bams : $(BAMS) $(addsuffix .bai,$(BAMS))
bam/%.bam : unprocessed_bam/%.$(BAM_SUFFIX)
	$(INIT) ln -f $< $@
else
ifeq ($(MERGE_SPLIT_BAMS),true)
merged_bams : $(BAMS) $(addsuffix .bai,$(BAMS))
bam/%.bam : unprocessed_bam/%$(if $(findstring true,$(BAM_FIX_RG)),.rg).bam
	$(INIT) ln -f $< $@
endif
endif

ifeq ($(MERGE_SPLIT_BAMS),true)
define bam-header
unprocessed_bam/$1.header.sam : $$(foreach split,$2,unprocessed_bam/$$(split).bam)
	$$(INIT) $$(SAMTOOLS) view -H $$< | grep -v '^@RG' > $$@.tmp; \
	for bam in $$(^M); do $$(SAMTOOLS) view -H $$$$bam | grep '^@RG' >> $$@.tmp; done; \
	uniq $$@.tmp > $$@ && $(RM) $$@.tmp
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call bam-header,$(sample),$(split.$(sample)))))

define merged-bam
unprocessed_bam/$1.bam : unprocessed_bam/$1.header.sam $$(foreach split,$2,unprocessed_bam/$$(split).bam)
	$$(call RUN,-s 12G -m 15G,"$$(SAMTOOLS) merge -f -h $$< $$@ $$(filter %.bam,$$^)")
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call merged-bam,$(sample),$(split.$(sample)))))
endif


%.bam.bai : %.bam
	$(call RUN,-c -s 4G -m 8G,"$(SAMTOOLS) index $<")

%.bai : %.bam.bai
	$(INIT) cp $< $@

%.filtered.bam : %.bam
	$(call RUN,-s 6G -m 7G,"$(SAMTOOLS) view -bF $(BAM_FILTER_FLAGS) $< > $@ && \
				$(RM) $<")

%.fixmate.bam : %.bam
	$(call RUN,-s 9G -m 14G,"$(FIX_MATE) \
				 INPUT=$< \
				 OUTPUT=$@ && \
				 $(RM) $<")

%.recal_report.grp : %.bam %.bai
	$(call RUN,-s 15G -m 15G,"$(call GATK_CMD,7G) \
				  -T BaseRecalibrator \
				  -R $(REF_FASTA) \
				  $(BAM_BASE_RECAL_OPTS) \
				  -I $< \
				  -o $@")

%.sorted.bam : %.bam
	$(call RUN,-s 30G -m 30G,"$(SORT_SAM) I=$< O=$@ SO=coordinate VERBOSITY=ERROR && \
				  $(RM) $<")

%.markdup.bam : %.bam
	$(call RUN,-s 24G -m 36G,"$(MKDIR) metrics && \
				  $(MARK_DUP) \
				  INPUT=$< \
				  OUTPUT=$@ \
				  METRICS_FILE=metrics/$(call strip-suffix,$(@F)).dup_metrics.txt && \
				  $(RM) $<")

%.rmdup.bam : %.bam
	$(call RUN,-s 4G -m 7G,"$(SAMTOOLS) rmdup $< $@ && \
				$(RM) $<")

%.clean.bam : %.bam
	$(call RUN,-s 6G -m 12G,"$(CLEAN_BAM) \
				 INPUT=$< \
				 OUTPUT=$@")

%.rg.bam : %.bam
	$(call RUN,-s 12G -m 16G,"$(ADD_RG) \
				  INPUT=$< \
				  OUTPUT=$@ \
				  RGLB=$(call strip-suffix,$(@F)) \
				  RGPL=$(SEQ_PLATFORM) RGPU=00000000 \
				  RGSM=$(call strip-suffix,$(@F)) \
				  RGID=$(call strip-suffix,$(@F)) \
				  VERBOSITY=ERROR && \
				  $(RM) $<")

ifeq ($(SPLIT_CHR),true)
define chr-target-realn
%.$1.chr_split.intervals : %.bam %.bam.bai
	$$(call RUN,-n 4 -s 4G -m 4G -v $(GATK_ENV),"$$(call GATK_CMD,5G) \
						     -T RealignerTargetCreator \
						     -I $$(<) \
						     -L $1 \
						     -nt 4 -R $$(REF_FASTA)  -o $$@ $$(BAM_REALN_TARGET_OPTS)")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-target-realn,$(chr))))

define chr-realn
%.$(1).chr_realn.bam : %.bam %.$(1).chr_split.intervals %.bam.bai
	$$(call RUN,-s 16G -m 16G  -v $(GATK_ENV),"if [[ -s $$(word 2,$$^) ]]; then $$(call GATK_CMD,4G) -T IndelRealigner \
						   -I $$(<) -R $$(REF_FASTA) -L $1 -targetIntervals $$(word 2,$$^) \
						   -o $$(@) $$(BAM_REALN_OPTS); \
						   else $$(call GATK_CMD,8G) -T PrintReads -R $$(REF_FASTA) -I $$< -L $1 -o $$@ ; fi")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-realn,$(chr))))

%.realn.bam : $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_realn.bam) $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_realn.bai)
	$(call RUN,-n 2 -s 10G -m 11G,"$(MERGE_SAM) $(foreach i,$(filter %.bam,$^), I=$(i)) SORT_ORDER=coordinate O=$@ USE_THREADING=true && \
				       $(RM) $^ $(@:.realn.bam=.bam)")

%.recal.bam : $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_recal.bam) $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_recal.bai)
	$(call RUN,-n 2 -s 10G -m 11G,"$(MERGE_SAM) $(foreach i,$(filter %.bam,$^), I=$(i)) SORT_ORDER=coordinate O=$@ USE_THREADING=true && \
				       $(RM) $^ $(@:.recal.bam=.bam)")

define chr-recal
%.$1.chr_recal.bam : %.bam %.recal_report.grp
	$$(call RUN,-s 11G -m 15G  -v $(GATK_ENV),"$$(call GATK_CMD,6G) -T PrintReads -L $1 -R $$(REF_FASTA) -I $$< -BQSR $$(<<) -o $$@")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-recal,$(chr))))

else

%.recal.bam : %.bam %.recal_report.grp
	$(call RUN,-s 14G -m 15G  -v $(GATK_ENV),"$(call GATK_CMD,7G) -T PrintReads -R $(REF_FASTA) -I $< -BQSR $(word 2,$^) -o $@ && \
						  $(RM) $<")

%.realn.bam : %.bam %.intervals %.bam.bai
	if [[ -s $(word 2,$^) ]]; then $(call RUN,-s 16G -m 16G -v $(GATK_ENV),"$(call GATK_CMD,7G) -T IndelRealigner \
										-I $< -R $(REF_FASTA) -targetIntervals $(<<) \
										-o $@ $(BAM_REALN_OPTS) && $(RM) $<") ; \
										else mv $< $@ ; fi

%.intervals : %.bam %.bam.bai
	$(call RUN,-n 4 -s 3G -m 3.5G -v $(GATK_ENV),"$(call GATK_CMD,6G) -T RealignerTargetCreator \
						      -I $< \
						      -nt 4 -R $(REF_FASTA) -o $@ $(BAM_REALN_TARGET_OPTS)")
endif

endif
PROCESS_BAM_MK = true
