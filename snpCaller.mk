# Run samtools mpileup on all snp positions and bams
##### DEFAULTS ######
REF ?= hg19
LOGDIR = log/snpCaller.$(NOW)

##### MAKE INCLUDES #####
include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc
VPATH ?= bam

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all

ifeq ($(EXOME),true)
DBSNP_SUBSET ?= $(HOME)/share/reference/dbsnp_137_exome.bed
else
DBSNP_SUBSET = $(HOME)/share/reference/dbsnp_tseq_intersect.bed
endif

CLUSTER_VCF = $(RSCRIPT) $(HOME)/share/scripts/clusterVcf.R

all : snp_vcf/snps_filtered.sdp_ft.clust.png

#snp_vcf/snps.vcf : $(foreach sample,$(SAMPLES),bam/$(sample).bam)
#$(call LSCRIPT_MEM,4G,8G,"$(SAMTOOLS) mpileup -f $(REF_FASTA) -g -l <(sed '/^#/d' $(DBSNP) | cut -f 1,2) $^ | $(BCFTOOLS) view -g - > $@")

snp_vcf/snps.vcf : $(foreach sample,$(SAMPLES),snp_vcf/$(sample).snps.vcf)
	$(call LSCRIPT_MEM,16G,20G,"$(call GATK_MEM,14G) -T CombineVariants $(foreach vcf,$^,--variant $(vcf) ) -o $@ --genotypemergeoption UNSORTED -R $(REF_FASTA)")

snp_vcf/snps_filtered.vcf : snp_vcf/snps.vcf
	$(INIT) grep '^#' $< > $@ && grep -e '0/1' -e '1/1' $< >> $@

snp_vcf/%.snps.vcf : bam/%.bam 
	$(call LSCRIPT_PARALLEL_MEM,4,2.5G,3G,"$(call GATK_MEM,8G) -T UnifiedGenotyper -nt 4 -R $(REF_FASTA) --dbsnp $(DBSNP) $(foreach bam,$(filter %.bam,$^),-I $(bam) ) -L $(DBSNP_SUBSET) -o $@ --output_mode EMIT_ALL_SITES")

snp_vcf/%.clust.png : snp_vcf/%.vcf
	$(INIT) $(CLUSTER_VCF) --outPrefix snp_vcf/$* $<

include ~/share/modules/vcftools.mk
