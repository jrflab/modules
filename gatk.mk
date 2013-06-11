# common GATK steps for variant calling
# Author: Raymond Lim <raylim@mm.st> & Fong Chun Chan <fongchunchan@gmail.com>
#
# GATK needs sorted, properly ordered bams
# INPUTS: sample bams
# OUTPUTS: vcf file for SNPs and indels
# OPTIONS: HARD_FILTER_SNPS = true/false (default: true)
# 		   POOL_SNP_RECAL = true/false (default: false)
# 		   SPLIT_CHR = true/false (default: true)
#

include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

HARD_FILTER_SNPS ?= true
POOL_SNP_RECAL ?= false
SPLIT_CHR ?= true

SNP_EFF_FLAGS ?= -ud 0 -no-intron -no-intergenic

###### RECIPES #######

# if SPLIT_CHR is set to true, we will split gatk processing by chromosome
ifeq ($(SPLIT_CHR),true)
# indel realignment intervals (i.e. where to do MSA)
# split by samples and chromosomes
# %=sample
# $(eval $(call chr-target-aln,chromosome))
define chr-target-realn
gatk/realn/%.$1.intervals : bam/%.bam bam/%.bam.bai
	$$(call INIT_PARALLEL_MEM,6,2G,2.5G) $$(call GATK_MEM,8G) -T RealignerTargetCreator \
	-I $$< \
	-L $1 \
	-nt 6 -R $$(REF_FASTA)  -o $$@  --known $$(KNOWN_INDELS) \
	&> $$(LOG)
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-target-realn,$(chr))))

# indel realignment per chromosome
# only realign if intervals is non-empty
# %=sample
# $(eval $(call chr-aln,chromosome))
#define chr-realn
#gatk/chr_bam/%.$(1).realn.bam : bam/%.bam gatk/realn/%.$(1).intervals bam/%.bam.bai
#	$$(call INIT_MEM,19G,25G) if [[ -s $$(word 2,$$^) ]]; then $$(GATK) -T IndelRealigner \
#	-I $$< -R $$(REF_FASTA) -L $1 -targetIntervals $$(word 2,$$^) \
#	-o $$@ --knownAlleles $$(KNOWN_INDELS) &> $$(LOG); \
#	else $$(GATK) -T PrintReads -R $$(REF_FASTA) -I $$< -L $1 -o $$@ &> $$(LOG) ; #fi
#endef
#$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-realn,$(chr))))
define chr-realn
gatk/chr_bam/%.$(1).realn.bam : bam/%.bam gatk/realn/%.$(1).intervals bam/%.bam.bai
	$$(call INIT_MEM,9G,12G) $$(call GATK_MEM,8G) -T IndelRealigner \
	 -I $$< -R $$(REF_FASTA) -L $1 -targetIntervals $$(word 2,$$^) \
	-o $$@ --knownAlleles $$(KNOWN_INDELS) &> $$(LOG)
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-realn,$(chr))))


# merge sample chromosome bams
define merge-chr
gatk/bam/$1.%.bam : $$(foreach chr,$$(CHROMOSOMES),gatk/chr_bam/$1.$$(chr).%.bam) $$(foreach chr,$$(CHROMOSOMES),gatk/chr_bam/$1.$$(chr).%.bai)
	$$(call INIT_PARALLEL_MEM,2,10G,11G) $$(MERGE_SAMS) $$(foreach i,$$(filter %.bam,$$^), I=$$i) SORT_ORDER=coordinate O=$$@ USE_THREADING=true &> $$(LOGDIR)/$$@.log && $$(RM) $$^ 
endef
$(foreach sample,$(SAMPLES),$(eval $(call merge-chr,$(sample))))

# $(eval $(call chr-variants,chromosome))
define chr-variants
gatk/chr_vcf/%.$1.variants.vcf : gatk/chr_bam/%.$1.realn.bam gatk/chr_bam/%.$1.realn.bai
	$$(call INIT_PARALLEL_MEM,2,5G,6G) $$(call GATK_MEM,8G) -T UnifiedGenotyper \
	-L $1 -nt 2 -I $$< --dbsnp $$(DBSNP) -o $$@  -glm BOTH -rf BadCigar \
	-stand_call_conf $$(VARIANT_CONF_THRESHOLD)  -R $$(REF_FASTA) &> $$(LOGDIR)/$$@.log
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call chr-variants,$(chr))))

# merge variants 
gatk/vcf/%.variants.vcf : $(foreach chr,$(CHROMOSOMES),gatk/chr_vcf/%.$(chr).variants.vcf)
	$(call INIT_MEM,2G,3G) $(call GATK_MEM,2G) -T CombineVariants --assumeIdenticalSamples $(foreach i,$^, --variant $i) -R $(REF_FASTA) -o $@ &> $(LOG)

#$(call INIT_MEM,2G,3G) grep '^#' $< > $@; cat $^ | grep -v '^#' | vcfsorter.pl $(REF_DICT) - >> $@ 2> $(LOGDIR)/$@.log

else # no splitting by chr
gatk/bam/%.realn.bam : bam/%.bam gatk/realn/%.intervals bam/%.bam.bai
	$(call INIT_MEM,9G,12G) $(call GATK_MEM,8G) -T IndelRealigner \
	 -I $< -R $(REF_FASTA) -targetIntervals $(word 2,$^) \
	 -o $@ --knownAlleles $(KNOWN_INDELS) &> $(LOG)
#gatk/bam/%.realn.bam : bam/%.bam gatk/realn/%.intervals bam/%.bam.bai
#	$(call INIT_MEM,19G,25G) if [[ -s $(word 2,$^) ]]; then $(GATK) -T IndelRealigner \
#	-I $< -R $(REF_FASTA) -targetIntervals $(word 2,$^) \
#	-o $@ --knownAlleles $(KNOWN_INDELS) &> $(LOG); \
#	else cp $< $@ &> $(LOG) ; fi

gatk/realn/%.intervals : bam/%.bam bam/%.bam.bai
	$(call INIT_PARALLEL_MEM,6,2G,3G) $(call GATK_MEM,8G) -T RealignerTargetCreator \
	-I $< \
	-nt 6 -R $(REF_FASTA)  -o $@  --known $(KNOWN_INDELS) \
	&> $(LOG)

gatk/vcf/%.variants.vcf : gatk/bam/%.realn.bam gatk/bam/%.realn.bai
	$(call INIT_PARALLEL_MEM,2,4G,6G) $(call GATK_MEM,8G) -T UnifiedGenotyper \
	-nt 2 -I $< --dbsnp $(DBSNP) -o $@  -glm BOTH -rf BadCigar \
	-stand_call_conf $(VARIANT_CONF_THRESHOLD)  -R $(REF_FASTA) &> $(LOG)
endif

# select snps % = sample
gatk/vcf/%.snps.vcf : gatk/vcf/%.variants.vcf gatk/vcf/%.variants.vcf.idx
	$(call INIT_MEM,8G,12G) $(call GATK_MEM,8G) -T SelectVariants  -R $(REF_FASTA)  --variant $<  -o $@ \
	 -selectType SNP &> $(LOGDIR)/$@.log

# select indels % = indels
gatk/vcf/%.indels.vcf : gatk/vcf/%.variants.vcf gatk/vcf/%.variants.vcf.idx
	 $(call INIT_MEM,8G,12G) $(call GATK_MEM,8G) -T SelectVariants -R $(REF_FASTA) --variant $<  -o $@ \
	-selectType INDEL &> $(LOG)

%.bai : %.bam
	$(call INIT_MEM,4G,8G) $(SAMTOOLS) index $< $@

$(REF_FASTA).fai : $(REF_FASTA)
	$(call INIT_MEM,4G,8G) $(SAMTOOLS) faidx $< || ($(RM) $@ && false)

$(REF_FASTA:.fasta=.dict) : $(REF_FASTA)
	$(call INIT_MEM,5G,6G) $(CREATE_SEQ_DICT) REFERENCE=$< OUTPUT=$@

# recalibrate base quality
%.recal_report.grp : %.bam %.bai
	$(call INIT_MEM,8G,20G) $(call GATK_MEM,8G) -T BaseRecalibrator -R $(REF_FASTA) -knownSites $(DBSNP) -I $< -o $@ &> $(LOGDIR)/$@.log

# recalibration
%.recal.bam : %.bam %.recal_report.grp %.bai
	$(call INIT_MEM,8G,12G) $(call GATK_MEM,8G) -T PrintReads -R $(REF_FASTA) -I $< -BQSR $(word 2,$^) -o $@ &> $(LOGDIR)/$@.log && $(RM) $<

# reduce bam size
%.reduced.bam : %.bam %.bai
	$(call INIT_MEM,8G,12G) $(call GATK_MEM,8G) -T ReduceReads -R $(REF_FASTA) -I $< -o $@ &> $(LOGDIR)/$@.log # && $(RM) $<


# variant recalibrator for snps
# not recommended for indels (not enough data)
#ifeq ($(POOL_SNP_RECAL),true)
#gatk/variant_recal/%snps.recal gatk/variant_recal/%snps.tranches gatk/variant_recal/%snps.plots.R : $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).%snps.annotated.vcf)
#SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,2G,4G) -pe $(PARALLEL_ENV) 6" $(MKDIR) $(@D); \
#$(GATK) -T VariantRecalibrator -R $(REF_FASTA) $(foreach i,$^, -input $i)  -nt 6 \
#-resource:hapmap,known=false,training=true,truth=true,prior=15.0 $(HAPMAP) \
#-resource:omni,known=false,training=true,truth=false,prior=12.0 $(OMNI) \
#-resource:dbsnp,known=true,training=false,truth=false,prior=8.0 $(DBSNP)  $(foreach i,$(VARIANT_RECAL_ANNOTATIONS), -an $i) \
#-recalFile gatk/variant_recal/$*snps.recal \
#-tranchesFile gatk/variant_recal/$*snps.tranches \
#-rscriptFile gatk/variant_recal/$*snps.plots.R &> $(LOGDIR)/$@.log || ($(RM) $@ && false)
#else 
#endif

VARIANT_RECAL = $(call INIT_PARALLEL_MEM,6,2.5G,3G) $(call GATK_MEM,8G) -T VariantRecalibrator -R $(REF_FASTA) -nt 6 \
				 -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $(HAPMAP) \
				 -resource:omni,known=false,training=true,truth=false,prior=12.0 $(OMNI) \
				 -resource:dbsnp,known=true,training=false,truth=false,prior=8.0 $(DBSNP) \
				 $(foreach i,$(VARIANT_RECAL_ANNOTATIONS), -an $i) \
				 $(foreach i,$(filter %.vcf,$^), -input $i) \
				 -recalFile $@ -tranchesFile $(basename $@).tranches -rscriptFile $(basename $@).snps.plots.R &> $(LOGDIR)/$@.log

# apply variant recal function
# arguments: vcf, recal file
APPLY_VARIANT_RECAL = $(call INIT_MEM,9G,12G) \
	$(call GATK_MEM,8G) -T ApplyRecalibration  -R $(REF_FASTA) -input $< \
	--ts_filter_level $(VARIANT_RECAL_TRUTH_SENSITIVITY_LEVEL)  -recalFile $(word 2,$^) \
	-tranchesFile $(basename $(word 2,$^)).tranches -o $@ &> $(LOGDIR)/$@.log

# apply variant recal %=sample
ifeq ($(HARD_FILTER_SNPS),true)
gatk/vcf/%.snps.filtered.vcf : gatk/vcf/%.snps.vcf gatk/vcf/%.snps.vcf.idx
	$(call INIT_MEM,9G,12G) \
	$(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) $(SNP_FILTERS) -o $@ \
	--variant $< &> $(LOGDIR)/$@.log
else 
ifeq ($(POOL_SNP_RECAL),true)
gatk/vcf/snps.recal.vcf : $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).snps.annotated.vcf) $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).snps.annotated.vcf.idx)
	$(VARIANT_RECAL)
gatk/vcf/%.snps.filtered.vcf : gatk/vcf/%.snps.vcf gatk/vcf/snps.recal.vcf gatk/vcf/snps.recal.vcf.idx gatk/vcf/%.snps.vcf.idx 
	$(APPLY_VARIANT_RECAL)
#gatk/variant_recal/no_realn.snps.recal : $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).snps.annotated.vcf)
#$(VARIANT_RECAL)
#gatk/variant_recal/%.snps.annotated.filtered.vcf : gatk/vcf/%.snps.annotated.vcf gatk/variant_recal/no_realn.snps.recal
#$(APPLY_VARIANT_RECAL)
else 
gatk/vcf/%.snps.recal.vcf : gatk/vcf/%.snps.vcf gatk/vcf/%.snps.vcf.idx
	$(VARIANT_RECAL)
gatk/vcf/%.snps.filtered.vcf : gatk/vcf/%.snps.vcf gatk/vcf/%.snps.recal.vcf gatk/vcf/%.snps.vcf.idx gatk/vcf/%.snps.recal.vcf.idx
	$(APPLY_VARIANT_RECAL)
endif
endif

# hard filter indels %=sample
gatk/vcf/%.indels.filtered.vcf : gatk/vcf/%.indels.vcf gatk/vcf/%.indels.vcf.idx
	$(call INIT_MEM,9G,12G) $(call GATK_MEM,8G) -T VariantFiltration -R $(REF_FASTA) $(INDEL_FILTERS) -o $@ \
	--variant $< &> $(LOGDIR)/$@.log || ($(RM) $@ && false)

# filter for only novel snps/indels
%.novel.txt : %.txt
	/bin/awk 'NR == 1 || $$4 == "."' $< > $@

vcf/%.gatk_snps.vcf : gatk/vcf/%.snps.filtered.vcf
	$(INIT) cp $< $@

vcf/%.gatk_indels.vcf : gatk/vcf/%.indels.filtered.vcf
	$(INIT) cp $< $@


include ~/share/modules/processBam.mk
include ~/share/modules/vcftools.mk
