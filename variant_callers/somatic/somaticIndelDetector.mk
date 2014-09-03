# GATK Somatic indel detector for tumour-normal matched samples
# vim: set ft=make :

include ~/share/modules/Makefile.inc
include ~/share/modules/variant_callers/gatk.inc

REF ?= hg19
LOGDIR = log/gatk_som_indel.$(NOW)
SAMPLE_PAIR_FILE ?= sample_pairs.txt
SAMPLE_FILE ?= samples.txt
TUMOR_SAMPLES ?= $(shell cut -f 1 $(SAMPLE_PAIR_FILE))
NORMAL_SAMPLES ?= $(shell cut -f 2 $(SAMPLE_PAIR_FILE))
SAMPLES ?= $(shell cat $(SAMPLE_FILE))
SPLIT_CHR ?= true

GATK_OLD_JAR = /opt/common/gatk/GenomeAnalysisTK-2.3-9-ge5ebf34/GenomeAnalysisTK.jar
SOM_INDEL = $(call JAVA) -Xmx8G -jar $(GATK_OLD_JAR) -T SomaticIndelDetector

DEPTH_FILTER = 10

VCF_SAMPLES = 0 1
VCF_GEN_IDS = GT AD DP MM MQS NQSBQ NQSMM REnd RStart SC
#VCF_ANNOTATIONS = RMSMappingQuality ReadPositionRankSumTest QualByDepth MappingQualityZero MappingQualityRankSumTest HaplotypeScore FisherStrand DepthOfCoverage BaseQualityRankSumTest HomopolymerRun

$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval normal_lookup.$(word $i,$(TUMOR_SAMPLES)) := $(word $i,$(NORMAL_SAMPLES))))
$(foreach i,$(shell seq 1 $(words $(TUMOR_SAMPLES))),$(eval tumor_lookup.$(word $i,$(NORMAL_SAMPLES)) := $(word $i,$(TUMOR_SAMPLES))))

VPATH = bam

LOGDIR = log/gatk.$(NOW)

INDEL_WINDOW_SIZE = 200

FILTER_SUFFIX := dp_ft.dbsnp.nsfp.hrun.eff

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all som_indel_vcfs som_indel_tables

all : som_indel_vcfs som_indel_tables
som_indel_vcfs : $(foreach tumor,$(TUMOR_SAMPLES),vcf/$(tumor)_$(normal_lookup.$(tumor)).gatk_som_indels.$(FILTER_SUFFIX).vcf)
som_indel_tables : $(foreach tumor,$(TUMOR_SAMPLES),tables/$(tumor)_$(normal_lookup.$(tumor)).gatk_som_indels.$(FILTER_SUFFIX).som_indel_ft.txt)
	
mutect_som_indel_tables : $(foreach tumor,$(TUMOR_SAMPLES),tables/$(tumor)_$(normal_lookup.$(tumor)).mutect_som_indels.txt)

ifeq ($(SPLIT_CHR),true)
#$(call som-indel,tumor,normal,chr)
define som-indel-tumor-normal-chr
gatk/chr_vcf/$1_$2.$3.som_indels.vcf : bam/$1.bam bam/$2.bam bam/$1.bai bam/$2.bai
	$$(call LSCRIPT_MEM,9G,14G,"$$(MKDIR) gatk/metrics; $$(SOM_INDEL) -R $$(REF_FASTA) -I:tumor $$(word 1,$$^) -I:normal $$(word 2,$$^) -o $$@ --metrics_file gatk/metrics/$1_$2.som_indels.grp -L $3 --window_size $$(INDEL_WINDOW_SIZE) &> $$(LOG)")
endef
$(foreach chr,$(CHROMOSOMES),$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call som-indel-tumor-normal-chr,$(tumor),$(normal_lookup.$(tumor)),$(chr)))))

define merge-som-indel-tumor-normal
gatk/vcf/$1_$2.som_indels.vcf : $$(foreach chr,$$(CHROMOSOMES),gatk/chr_vcf/$1_$2.$$(chr).som_indels.vcf)
	$$(INIT) grep '^#' $$< > $$@; cat $$^ | grep -v '^#' | vcfsorter.pl $$(REF_DICT) - >> $$@ 2> $$(LOG)
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call merge-som-indel-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))

else # no splitting by chromosome

#$(call som-indel-tumor-normal,tumor,normal)
define som-indel-tumor-normal
gatk/chr_vcf/$1_$2.som_indels.vcf : bam/$1.bam bam/$2.bam bam/$1.bai bam/$2.bai
	$$(call LSCRIPT_MEM,9G,14G,"$$(MKDIR) gatk/metrics; $$(SOM_INDEL) -R $$(REF_FASTA) -I:tumor $$(word 1,$$^) -I:normal $$(word 2,$$^) -o $$@ --metrics_file gatk/metrics/$1_$2.som_indels.grp --window_size $$(INDEL_WINDOW_SIZE) &> $$(LOG)")
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call som-indel-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))
endif

#$(call pedigree-som-indel-tumor-normal,tumor,normal)
define pedigree-som-indel-tumor-normal
vcf/$1_$2.gatk_som_indels.vcf : gatk/vcf/$1_$2.som_indels.vcf
	$$(INIT) grep '^##' $$< > $$@; echo "##PEDIGREE=<Derived=$1,Original=$2>" >> $$@; grep '^#[^#]' $$< >> $$@; cat $$^ | grep -v '^#' >> $$@ 2> $$(LOG)
endef
$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call pedigree-som-indel-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))

tables/%.som_indel_ft.txt : tables/%.txt
	$(INIT) head -1 $< > $@; sed '1d' $< | awk -F$$'\t' '$$8 != $$9 && $$11 >= $(DEPTH_FILTER) && $$12 >= $(DEPTH_FILTER)' >> $@

tables/%.mutect_som_indels.txt : tables/%.mutect.dp_ft.annotated.nsfp.pass.novel.txt tables/%.gatk_som_indels.annotated.nsfp.som_indel_ft.txt
	$(call LSCRIPT,"$(RSCRIPT) $(RBIND) $^ > $@")


include ~/share/modules/variant_callers/gatk.mk
