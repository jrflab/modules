# various bam processing steps
##### MAKE INCLUDES #####
include ~/share/modules/gatk.inc

SAMPLE_FILE ?= samples.txt
SAMPLES ?= $(shell cat $(SAMPLE_FILE))

SPLIT_CHR ?= true

# not primary alignment
# read fails platform/vendor quality checks
BAM_FILTER_FLAGS ?= 768

.DELETE_ON_ERROR:
.SECONDARY: 

# indices
# if bam file is a symlink, need to create a symlink to index
%.bam.bai : %.bam.md5
	$(call LSCRIPT_MEM,4G,8G,"$(CHECK_MD5) $(SAMTOOLS) index $(<:.md5=)")

%.bai : %.bam.md5
	$(call LSCRIPT_MEM,4G,8G,"$(CHECK_MD5) $(SAMTOOLS) index $(<:.md5=) $@")

# sam to bam
#%.bam : %.sam
#	$(call INIT_MEM,2G,3G) $(SAMTOOLS) view -bS $< > $@

# filter
%.filtered.bam.md5 : %.bam.md5
	$(call LSCRIPT_MEM,6G,7G,"$(CHECK_MD5) $(SAMTOOLS) view -bF $(BAM_FILTER_FLAGS) $(<:.md5=) > $(@:.md5=) 2> $(LOG) && $(MD5) && $(RM) $(<:.md5=) $<")

%.fixmate.bam : %.bam
	$(call LSCRIPT_MEM,9G,14G,"$(call FIX_MATE_MEM,8G) I=$< O=$@ &> $(LOG) && $(RM) $<")

# recalibrate base quality
%.recal_report.grp : %.bam.md5 %.bai
	$(call LSCRIPT_MEM,11G,15G,"$(CHECK_MD5) $(call GATK_MEM,10G) -T BaseRecalibrator -R $(REF_FASTA) -knownSites $(DBSNP) -I $(<:.md5=) -o $@ &> $(LOG)")

# recalibration
%.recal.bam.md5 : %.bam.md5 %.recal_report.grp
	$(call LSCRIPT_MEM,11G,15G,"$(CHECK_MD5) $(call GATK_MEM,10G) -T PrintReads -R $(REF_FASTA) -I $(<:.md5=) -BQSR $(word 2,$^) -o $(@:.md5=) &> $(LOG) && $(MD5) && $(RM) $< $(<:.md5=)")

# sort only if necessary
#%.sorted.bam : %.bam
#	$(call INIT_MEM,20G,25G) if ! $(SAMTOOLS) view -H $< | grep -q 'SO:coordinate' -; then $(call SORT_SAM_MEM,19G) I=$< O=$@ SO=coordinate &> $(LOG); else ln -v $< $@; fi && $(RM) $<

%.sorted.bam.md5 : %.bam.md5
	$(call LSCRIPT_MEM,9G,14G,"$(CHECK_MD5) $(call SORT_SAM_MEM,8G) I=$(<:.md5=) O=$(@:.md5=) SO=coordinate &> $(LOG) && $(MD5) && $(RM) $< $(<:.md5=)")

# reorder only if necessary
#%.reord.bam : %.bam
#	if ! $(SAMTOOLS) view -H $< | grep '@SQ' | cut -f 2 | sed 's/^SN://' | diff -q - $(REF_CHR_ORD) > /dev/null; then $(call LAUNCH_MEM,6G,7G) $(call REORDER_SAM_MEM,6G) I=$< O=$@ R=$(REF_FASTA) &> $(LOG); else ln -v $< $@; fi && $(RM) $<


# mark duplicates
%.markdup.bam.md5 : %.bam.md5
	$(call LSCRIPT_MEM,14G,18G,"$(CHECK_MD5) $(call MARK_DUP_MEM,14G) I=$(<:.md5=) O=$(@:.md5=) METRICS_FILE=metrics/$(call strip-suffix,$(@F:.md5=)).dup_metrics.txt &> $(LOG) && $(MD5) && $(RM) $< $(<:.md5=)")

%.rmdup.bam.md5 : %.bam.md5
	$(call LSCRIPT_MEM,4G,7G,"$(CHECK_MD5) $(SAMTOOLS) rmdup $(<:.md5=) $(@:.md5=) &> $(LOG) && $(MD5) && $(RM) $< $(<:.md5=)")

# clean sam files
%.clean.bam.md5 : %.bam.md5
	$(call LSCRIPT_MEM,6G,12G,"$(CHECK_MD5) $(call CLEANBAM_MEM,6G) I=$(<:.md5=) O=$(@:.md5=) &> $(LOG) && $(MD5)")

# add rg
%.rg.bam.md5 : %.bam.md5
	$(call LSCRIPT_MEM,10G,12G,"$(call ADD_RG_MEM,10G) I=$(<:.md5=) O=$(@:.md5=) RGLB=$(call strip-suffix,$(@F)) RGPL=illumina RGPU=00000000 RGSM=$(call strip-suffix,$(@F)) RGID=$(call strip-suffix,$(@F)) &> $(LOG) && $(RM) $< $(<:.md5=) && $(MD5)")


# if SPLIT_CHR is set to true, we will split realn processing by chromosome
ifeq ($(SPLIT_CHR),true)
# indel realignment intervals (i.e. where to do MSA)
# split by samples and chromosomes
# %=sample
# $(eval $(call chr-target-aln,chromosome))
define chr-target-realn
%.$1.chr_split.intervals : %.bam.md5 %.bam.bai
	$$(call LSCRIPT_PARALLEL_MEM,4,2.5G,3G,"$$(CHECK_MD5) $$(call GATK_MEM,8G) -T RealignerTargetCreator \
	-I $$(<:.md5=) \
	-L $1 \
	-nt 4 -R $$(REF_FASTA)  -o $$@  --known $$(KNOWN_INDELS) \
	&> $$(LOG)")
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
	-o $$(@:.md5=) --knownAlleles $$(KNOWN_INDELS) &> $$(LOG); \
	else $$(call GATK_MEM,8G) -T PrintReads -R $$(REF_FASTA) -I $$(<:.md5=) -L $1 -o $$(@:.md5=) &> $$(LOG) ; fi && $$(MD5)")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-realn,$(chr))))

# merge sample chromosome bams
%.realn.bam.md5 : $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_realn.bam.md5) $(foreach chr,$(CHROMOSOMES),%.$(chr).chr_realn.bai)
	$(call LSCRIPT_PARALLEL_MEM,2,10G,11G,"$(CHECK_MD5) $(MERGE_SAMS) $(foreach i,$(filter %.bam.md5,$^), I=$(i:.md5=)) SORT_ORDER=coordinate O=$(@:.md5=) USE_THREADING=true &> $(LOG) && $(MD5) && $(RM) $^ $(^:.md5=)")

else # no splitting by chr
%.realn.bam.md5 : %.bam.md5 %.intervals %.bam.bai
	if [[ -s $(word 2,$^) ]]; then $(call LSCRIPT_MEM,9G,12G,"$(CHECK_MD5) $(call GATK_MEM,8G) -T IndelRealigner \
	-I $(<:.md5=) -R $(REF_FASTA) -targetIntervals $(word 2,$^) \
	-o $(@:.md5=) --knownAlleles $(KNOWN_INDELS) &> $(LOG) && $(MD5) && $(RM) $< $(<:.md5=)") ; \
	else mv $(<:.md5=) $(@:.md5=) && $(MD5) &> $(LOG) ; fi

%.intervals : %.bam.md5 %.bam.bai
	$(call LSCRIPT_PARALLEL_MEM,4,2.5G,3G,"$(CHECK_MD5) $(call GATK_MEM,8G) -T RealignerTargetCreator \
	-I $(<:.md5=) \
	-nt 4 -R $(REF_FASTA)  -o $@  --known $(KNOWN_INDELS) \
	&> $(LOG)")
endif

