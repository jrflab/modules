# various bam processing steps
##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc


LOGDIR ?= log/process_bam.$(NOW)

DUP_TYPE ?= markdup
NO_RECAL ?= false
NO_REALN ?= false
SPLIT_CHR ?= true
SPLIT_FASTQ ?= false
SPLIT_SORT ?= false
MERGE_SPLIT_BAMS ?= false # merge processed split bams
NUM_SORT_SPLITS ?= 50
SORT_SPLIT_SEQ = $(shell seq 0 $$(($(NUM_SORT_SPLITS) - 1)))

ifneq ($(KNOWN_INDELS),)
REALN_OPTS = --knownAlleles $(KNOWN_INDELS)
REALN_TARGET_OPTS = --known $(KNOWN_INDELS)
endif

# not primary alignment
# read fails platform/vendor quality checks
BAM_FILTER_FLAGS ?= 768


.DELETE_ON_ERROR:
.SECONDARY: 

REPROCESS ?= false
ifeq ($(REPROCESS),true)
BAM_SUFFIX := .sorted.filtered
ifeq ($(NO_REALN),false)
BAM_SUFFIX := $(BAM_SUFFIX).realn
endif

ifeq ($(DUP_TYPE),rmdup)
BAM_SUFFIX := $(BAM_SUFFIX).rmdup
else ifeq ($(DUP_TYPE),markdup) 
BAM_SUFFIX := $(BAM_SUFFIX).markdup
endif

ifeq ($(NO_RECAL),false)
BAM_SUFFIX := $(BAM_SUFFIX).recal
endif
endif

FIX_RG ?= false
ifeq ($(FIX_RG),true)
BAM_SUFFIX := $(BAM_SUFFIX).rg
endif

ifdef BAM_SUFFIX
BAM_SUFFIX := $(BAM_SUFFIX).bam
BAMS = $(foreach sample,$(SAMPLES),bam/$(sample).bam)
processed_bams : $(addsuffix .md5,$(BAMS)) $(addsuffix .bai,$(BAMS))

bam/%.bam.md5 : unprocessed_bam/%$(BAM_SUFFIX).md5
	$(INIT) cp $< $@ && ln -f $(<:.md5=) $(@:.md5=)
endif


ifeq ($(MERGE_SPLIT_BAMS),true)
define bam-header
unprocessed_bam/$1.header.sam : $$(foreach split,$2,unprocessed_bam/$$(split).bam.md5)
	$$(INIT) $$(SAMTOOLS) view -H $$(<M) | grep -v '^@RG' > $$@.tmp; \
	for bam in $$(^M); do $$(SAMTOOLS) view -H $$$$bam | grep '^@RG' >> $$@.tmp; done; \
	uniq $$@.tmp > $$@ && $(RM) $$@.tmp
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call bam-header,$(sample),$(split_lookup.$(sample)))))

define merged-bam
unprocessed_bam/$1.bam.md5 : unprocessed_bam/$1.header.sam $$(foreach split,$2,unprocessed_bam/$$(split).bam.md5)
	$$(call LSCRIPT_MEM,12G,15G,"$$(SAMTOOLS) merge -f -h $$(<M) $$(@M) $$(filter %.bam,$$(^M)) && $$(MD5)")
endef
$(foreach sample,$(SPLIT_SAMPLES),$(eval $(call merged-bam,$(sample),$(split_lookup.$(sample)))))
endif


# indices
# if bam file is a symlink, need to create a symlink to index
%.bam.bai : %.bam.md5
	$(call LSCRIPT_CHECK_MEM,4G,8G,"$(CHECK_MD5) $(SAMTOOLS) index $(<M)")

%.bai : %.bam.md5
	$(call LSCRIPT_CHECK_MEM,4G,8G,"$(CHECK_MD5) $(SAMTOOLS) index $(<M) $@")

%.bam.md5 : %.bam
	$(call LSCRIPT_CHECK,"$(MD5)")

# sam to bam
#%.bam : %.sam
#	$(call INIT_MEM,2G,3G) $(SAMTOOLS) view -bS $< > $@

# filter
%.filtered.bam.md5 : %.bam.md5
	$(call LSCRIPT_MEM,6G,7G,"$(CHECK_MD5) $(SAMTOOLS) view -bF $(BAM_FILTER_FLAGS) $(<M) > $(@M) && $(MD5) && $(RM) $(<M) $<")

%.fixmate.bam : %.bam
	$(call LSCRIPT_MEM,9G,14G,"$(call FIX_MATE_MEM,8G) I=$< O=$@ && $(RM) $<")

# recalibrate base quality
%.recal_report.grp : %.bam.md5 %.bai
	$(call LSCRIPT_MEM,11G,15G,"$(CHECK_MD5) $(call GATK_MEM,10G) -T BaseRecalibrator -R $(REF_FASTA) -knownSites $(DBSNP) -I $(<:.md5=) -o $@")

# recalibration
%.recal.bam.md5 : %.bam.md5 %.recal_report.grp
	$(call LSCRIPT_MEM,11G,15G,"$(CHECK_MD5) $(call GATK_MEM,10G) -T PrintReads -R $(REF_FASTA) -I $(<M) -BQSR $(word 2,$^) -o $(@:.md5=) && $(MD5) && $(RM) $< $(<:.md5=)")

#sort only if necessary
#%.sorted.bam.md5 : %.bam.md5
#	 if ! $(SAMTOOLS) view -H $(<M) | grep -q 'SO:coordinate' -; then $(call LSCRIPT_MEM,20G,25G,"$(CHECK_MD5) $(call SORT_SAM_MEM,19G) I=$(<M) O=$(@M) SO=coordinate"); else cp $< $@ && ln -v $(<M) $(@M); fi && $(RM) $<

#define split-sort
#%.sorted_split_$1.bam.md5 : %.bam.md5
#$$(call LSCRIPT_MEM,17G,19G,"$$(CHECK_MD5) cat <($$(SAMTOOLS) view -H $$(<M)) <($$(SAMTOOLS) view $$(<M) | awk 'NR % $$(NUM_SORT_SPLITS) == $1') | $$(call SORT_SAM_MEM,15G,3000000) I=/dev/stdin O=$$(@M) SO=coordinate && $$(MD5)")
#endef

ifeq ($(SPLIT_SORT),true)
define split-sort
%.sorted_split_$1.bam.md5 : %.bam.md5
	$$(call LSCRIPT_MEM,7G,10G,"$$(CHECK_MD5) cat <($$(SAMTOOLS) view -H $$(<M)) <($$(SAMTOOLS) view $$(<M) | awk 'NR % $$(NUM_SORT_SPLITS) == $1') | $$(SAMTOOLS) view -Sub - | $$(SAMTOOLS) sort -m 4G - $$(@:.bam.md5=) && $$(MD5)")
endef
$(foreach i,$(SORT_SPLIT_SEQ),$(eval $(call split-sort,$i)))

%.sorted.bam.md5 : $(foreach i,$(SORT_SPLIT_SEQ),%.sorted_split_$i.bam.md5)
	$(call LSCRIPT_MEM,8G,10G,"$(CHECK_MD5) $(SAMTOOLS) merge -h <($(SAMTOOLS) view -H $(<M)) $(@M) $(^M) && $(MD5) && $(RM) $^ $(^M) $(@:.sorted.bam=.bam) $(@:.sorted.bam.md5=.bam.md5)")
else
%.sorted.bam.md5 : %.bam.md5
	$(call LSCRIPT_MEM,20G,25G,"$(CHECK_MD5) $(call SORT_SAM_MEM,19G,4500000) I=$(<:.md5=) O=$(@:.md5=) SO=coordinate && $(MD5) && $(RM) $< $(<:.md5=)")
endif


# reorder only if necessary
#%.reord.bam : %.bam
#	if ! $(SAMTOOLS) view -H $< | grep '@SQ' | cut -f 2 | sed 's/^SN://' | diff -q - $(REF_CHR_ORD) > /dev/null; then $(call LAUNCH_MEM,6G,7G) $(call REORDER_SAM_MEM,6G) I=$< O=$@ R=$(REF_FASTA) &> $(LOG); else ln -v $< $@; fi && $(RM) $<


# mark duplicates
%.markdup.bam.md5 : %.bam.md5
	$(call LSCRIPT_MEM,14G,18G,"$(CHECK_MD5) $(MKDIR) metrics; $(call MARK_DUP_MEM,14G) I=$(<:.md5=) O=$(@:.md5=) METRICS_FILE=metrics/$(call strip-suffix,$(@F:.md5=)).dup_metrics.txt && $(MD5) && $(RM) $< $(<:.md5=)")

%.rmdup.bam.md5 : %.bam.md5
	$(call LSCRIPT_MEM,4G,7G,"$(CHECK_MD5) $(SAMTOOLS) rmdup $(<:.md5=) $(@:.md5=) && $(MD5) && $(RM) $< $(<:.md5=)")

# clean sam files
%.clean.bam.md5 : %.bam.md5
	$(call LSCRIPT_MEM,6G,12G,"$(CHECK_MD5) $(call CLEANBAM_MEM,6G) I=$(<:.md5=) O=$(@:.md5=) && $(MD5)")

# add rg
%.rg.bam.md5 : %.bam.md5
	$(call LSCRIPT_MEM,12G,16G,"$(call ADD_RG_MEM,10G) I=$(<:.md5=) O=$(@:.md5=) RGLB=$(call strip-suffix,$(@F)) RGPL=illumina RGPU=00000000 RGSM=$(call strip-suffix,$(@F)) RGID=$(call strip-suffix,$(@F)) && $(RM) $< $(<:.md5=) && $(MD5)")


# if SPLIT_CHR is set to true, we will split realn processing by chromosome
ifeq ($(SPLIT_CHR),true)
# indel realignment intervals (i.e. where to do MSA)
# split by samples and chromosomes
# %=sample
# $(eval $(call chr-target-aln,chromosome))
define chr-target-realn
%.$1.chr_split.intervals : %.bam.md5 %.bam.bai
	$$(call LSCRIPT_PARALLEL_MEM,4,3G,4G,"$$(CHECK_MD5) \
		$$(call GATK_MEM,11G) -T RealignerTargetCreator \
		-I $$(<:.md5=) \
		-L $1 \
		-nt 4 -R $$(REF_FASTA)  -o $$@ $$(REALN_TARGET_OPTS)")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-target-realn,$(chr))))

# indel realignment per chromosome
# only realign if intervals is non-empty
# %=sample
# $(eval $(call chr-aln,chromosome))
define chr-realn
%.$(1).chr_realn.bam.md5 : %.bam.md5 %.$(1).chr_split.intervals %.bam.bai
	$$(call LSCRIPT_MEM,9G,12G,"$$(CHECK_MD5) if [[ -s $$(word 2,$$^) ]]; then $$(call GATK_MEM,8G) -T IndelRealigner \
	-I $$(<:.md5=) -R $$(REF_FASTA) -L $1 -targetIntervals $$(word 2,$$^) \
	-o $$(@:.md5=) $$(REALN_OPTS); \
	else $$(call GATK_MEM,8G) -T PrintReads -R $$(REF_FASTA) -I $$(<:.md5=) -L $1 -o $$(@:.md5=) ; fi && $$(MD5)")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-realn,$(chr))))

# merge sample chromosome bams
%.realn.bam.md5 : $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_realn.bam.md5) $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_realn.bai)
	$(call LSCRIPT_PARALLEL_MEM,2,10G,11G,"$(CHECK_MD5) $(MERGE_SAMS) $(foreach i,$(filter %.bam.md5,$^), I=$(i:.md5=)) SORT_ORDER=coordinate O=$(@:.md5=) USE_THREADING=true && $(MD5) && $(RM) $^ $(^:.md5=) $(@:.realn.bam.md5=.bam) $(@:.realn.bam.md5=.bam.md5)")

else # no splitting by chr
%.realn.bam.md5 : %.bam.md5 %.intervals %.bam.bai
	if [[ -s $(word 2,$^) ]]; then $(call LSCRIPT_MEM,9G,12G,"$(CHECK_MD5) $(call GATK_MEM,8G) -T IndelRealigner \
	-I $(<:.md5=) -R $(REF_FASTA) -targetIntervals $(word 2,$^) \
	-o $(@:.md5=) $(REALN_OPTS) && $(MD5) && $(RM) $< $(<:.md5=)") ; \
	else mv $(<:.md5=) $(@:.md5=) && $(MD5) ; fi

%.intervals : %.bam.md5 %.bam.bai
	$(call LSCRIPT_PARALLEL_MEM,4,2.5G,3G,"$(CHECK_MD5) $(call GATK_MEM,8G) -T RealignerTargetCreator \
	-I $(<:.md5=) \
	-nt 4 -R $(REF_FASTA)  -o $@ $(REALN_TARGET_OPTS)")
endif

