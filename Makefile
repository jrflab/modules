# Top-level Makefile
#
#
# Author: Raymond Lim <raylim@mm.st>
#

ifneq ("$(wildcard config.inc)", "")
include config.inc
endif
ifneq ("$(wildcard project_config.inc)", "")
include project_config.inc
endif
include modules/config.inc

export

NUM_ATTEMPTS ?= 20
NOW := $(shell date +"%F")
MAKELOG = log/$(@).$(NOW).log

USE_CLUSTER ?= true
ifeq ($(USE_CLUSTER),true)
MAKE = modules/scripts/qmake.pl -n $@.$(NOW) $(if $(SLACK_CHANNEL),-c $(SLACK_CHANNEL)) -r $(NUM_ATTEMPTS) -m -s -- make
endif
NUM_JOBS ?= 100

define RUN_MAKE_J
$(MAKE) -e -f $1 -j $2 $(TARGET) && \
	mkdir -p completed_tasks && \
	touch completed_tasks/$@
endef

RUN_MAKE = $(call RUN_MAKE_J,$1,$(NUM_JOBS))


TARGETS += merge_fastq
merge_fastq : 
	$(call RUN_MAKE,modules/fastq_tools/mergeFastq.mk)

TARGETS += gatk
gatk : 
	$(call RUN_MAKE,modules/variant_callers/gatkVariantCaller.mk)

TARGETS += bwa
bwa : NUM_ATTEMPTS = 50
bwa :
	$(call RUN_MAKE,modules/aligners/bwaAligner.mk)

TARGETS += bwamem
bwamem :
	$(call RUN_MAKE,modules/aligners/bwamemAligner.mk)

TARGETS += bowtie
bowtie : NUM_ATTEMPTS = 50
bowtie :
	$(call RUN_MAKE,modules/aligners/bowtieAligner.mk)

TARGETS += tmap
tmap : NUM_ATTEMPTS = 50
tmap :
	$(call RUN_MAKE,modules/aligners/tmapAligner.mk)

TARGETS += tophat_fusion
tophat_fusion : 
	$(call RUN_MAKE,modules/sv_callers/tophatFusion.mk)

TARGETS += tophat
tophat : 
	$(call RUN_MAKE,modules/aligners/tophatAligner.mk)

TARGETS += hisat
hisat : 
	$(call RUN_MAKE,modules/aligners/hisatAligner.mk)

TARGETS += cufflinks
cufflinks : 
	$(call RUN_MAKE,modules/rnaseq/cufflinks.mk)

TARGETS += process_bam
process_bam : 
	$(call RUN_MAKE,modules/bam_tools/processBam.mk)

TARGETS += bam_metrics
bam_metrics :
	$(call RUN_MAKE,modules/qc/bamMetrics.mk)

TARGETS += bam_interval_metrics
bam_interval_metrics :
	$(call RUN_MAKE,modules/qc/bamIntervalMetrics.mk)

TARGETS += rnaseq_metrics
rnaseq_metrics :
	$(call RUN_MAKE,modules/qc/rnaseqMetrics.mk)

TARGETS += fastqc
fastqc :
	$(call RUN_MAKE,modules/qc/fastqc.mk)

# not tested on the cluster
# requires x11 for graphics
TARGETS += interval_qc
interval_qc :
	$(call RUN_MAKE,modules/qc/intervalBamQC.mk)

TARGETS += jsm
jsm :
	$(call RUN_MAKE,modules/variant_callers/somatic/jsm.mk)

TARGETS += mutect
mutect :
	$(call RUN_MAKE,modules/variant_callers/somatic/mutectVariantCaller.mk)

TARGETS += varscan_cnv
varscan_cnv :
	$(call RUN_MAKE,modules/copy_number/varscanCNV.mk)

TARGETS += varscan_fpfilter
varscan_fpfilter :
	$(call RUN_MAKE,modules/variant_callers/varscanFpfilter.mk)

TARGETS += varscanTN
varscanTN :
	$(call RUN_MAKE,modules/variant_callers/somatic/varscanTN.mk)

TARGETS += varscan
varscan :
	$(call RUN_MAKE,modules/variant_callers/varscan.mk)


# single sample mutation seq
TARGETS += museqTN
museqTN :
	$(call RUN_MAKE,modules/variant_callers/somatic/museqTN.mk)

TARGETS += merge_vcfTN
merge_vcfTN :
	$(call RUN_MAKE,modules/vcf_tools/vcfMergeTN.mk)

TARGETS += somatic_sniper
somatic_sniper :
	$(call RUN_MAKE,modules/variant_callers/somatic/somaticSniper.mk)


TARGETS += som_indels
som_indels :
	$(call RUN_MAKE,modules/variant_callers/somatic/somaticIndelDetector.mk)

TARGETS += snvmix
snvmix :
	$(call RUN_MAKE,modules/variant_callers/snvmix.mk)

TARGETS += pyrohmm
pyrohmm :
	$(call RUN_MAKE,modules/variant_callers/pyroHMMVar.mk)


TARGETS += compare_vcf
compare_vcf :
	$(call RUN_MAKE,modules/vcf_tools/vcfCompare.mk)

TARGETS += merge_vcf_platform
merge_vcf_platform :
	$(call RUN_MAKE,modules/vcf_tools/vcfMergePlatform.mk)

TARGETS += compare_vcfTN
compare_vcfTN :
	$(call RUN_MAKE,modules/vcf_tools/vcfCompareTN.mk)

TARGETS += read_depth
read_depth :
	$(call RUN_MAKE,modules/qc/readDepth.mk)

TARGETS += qualimap
qualimap :
	$(call RUN_MAKE,modules/qc/qualimap.mk)

TARGETS += hmmcopy
hmmcopy :
	$(call RUN_MAKE,modules/copy_number/hmmCopy.mk)

TARGETS += nfuse_wgss_wtss
nfuse_wgss_wtss :
	$(call RUN_MAKE,modules/copy_number/nfuseWGSSWTSS.mk)

TARGETS += sum_reads
sum_reads :
	$(call RUN_MAKE,modules/rnaseq/sumRNASeqReads.mk)


TARGETS += exomecnv
exomecnv : 
	$(call RUN_MAKE,modules/copy_number/exomeCNV.mk)

TARGETS += exomecnvloh
exomecnvloh : 
	$(call RUN_MAKE,modules/copy_number/exomeCNVLOH.mk)

TARGETS += gistic
gistic :
	$(call RUN_MAKE,modules/copy_number/gistic.mk)

TARGETS += freec
freec : 
	$(call RUN_MAKE,modules/copy_number/controlFreeC.mk)

TARGETS += freecTN
freecTN : 
	$(call RUN_MAKE,modules/copy_number/controlFreeCTN.mk)

TARGETS += freec_lohTN
freec_lohTN : 
	$(call RUN_MAKE,modules/copy_number/controlFreeCLOHTN.mk)

NUM_DEFUSE_JOBS ?= 5
TARGETS += defuse
defuse :
	$(call RUN_MAKE_J,modules/sv_callers/defuse.mk,$(NUM_DEFUSE_JOBS))

NUM_CHIMSCAN_JOBS ?= 5
TARGETS += chimscan
chimscan :
	$(call RUN_MAKE_J,modules/sv_callers/chimerascan.mk,$(NUM_CHIMSCAN_JOBS))

TARGETS += oncofuse
oncofuse :
	$(call RUN_MAKE,modules/sv_callers/oncofuse.mk)

TARGETS += lumpy
lumpy :
	$(call RUN_MAKE,modules/sv_callers/lumpy.mk)

TARGETS += hydra
hydra :
	$(call RUN_MAKE,modules/sv_callers/hydra.mk)

TARGETS += pindel
pindel :
	$(call RUN_MAKE,modules/variant_callers/pindel.mk)

TARGETS += scalpel
scalpel :
	$(call RUN_MAKE,modules/variant_callers/somatic/scalpel.mk)

TARGETS += snp6
snp6 :
	$(call RUN_MAKE,modules/snp6/snp6.mk)

TARGETS += snpcaller
snpcaller :
	$(call RUN_MAKE,modules/variant_callers/snpCaller.mk)

TARGETS += soapfuse
soapfuse :
	$(call RUN_MAKE,modules/sv_callers/soapFuse.mk)

TARGETS += mapsplice
mapsplice :
	$(call RUN_MAKE,modules/sv_callers/mapsplice.mk)

TARGETS += fusioncatcher
fusioncatcher :
	$(call RUN_MAKE,modules/sv_callers/fusioncatcher.mk)

TARGETS += strelka
strelka :
	$(call RUN_MAKE,modules/variant_callers/somatic/strelka.mk)

TARGETS += crest
crest :
	$(call RUN_MAKE,modules/sv_callers/crest.mk)

TARGETS += pyloh
pyloh :
	$(call RUN_MAKE,modules/ploidy/pyloh.mk)

TARGETS += clonehd
clonehd :
	$(call RUN_MAKE,modules/clonality/clonehd.mk)

TARGETS += emu
emu :
	$(call RUN_MAKE,modules/mut_sigs/emu.mk)


TARGETS += extract_fastq
extract_fastq :
	$(call RUN_MAKE,modules/fastq_tools/extractFastq.mk)

TARGETS += titan
titan :
	$(call RUN_MAKE,modules/copy_number/titan.mk)

TARGETS += ann_titan
ann_titan :
	$(call RUN_MAKE,modules/copy_number/annotateTitan.mk)

TARGETS += samtools_het
samtools_het :
	$(call RUN_MAKE,modules/variant_callers/samtoolsHet.mk)

TARGETS += absolute_seq
absolute_seq :
	$(call RUN_MAKE,modules/clonality/absoluteSeq.mk)

TARGETS += merge_strelka_varscan
merge_strelka_varscan :
	$(call RUN_MAKE,modules/variant_callers/somatic/mergeStrelkaVarscanIndels.mk)

TARGETS += rseqc
rseqc :
	$(call RUN_MAKE,modules/qc/rseqc.mk)

TARGETS += integrate_rnaseq
integrate_rnaseq :
	$(call RUN_MAKE,modules/sv_callers/integrateRnaseq.mk)

TARGETS += integrate
integrate :
	$(call RUN_MAKE,modules/sv_callers/integrate.mk)


TARGETS += merge_split_fastq
merge_split_fastq :
	$(call RUN_MAKE,modules/fastq_tools/mergeSplitFastq.mk)

TARGETS += contest
contest :
	$(call RUN_MAKE,modules/contamination/contest.mk)

TARGETS += virus_detection_bowtie2
virus_detection_bowtie2 :
	$(call RUN_MAKE,modules/virus/virus_detection_bowtie2.mk)

TARGETS += fix_rg
fix_rg :
	$(call RUN_MAKE,modules/bam_tools/fixRG.mk)

TARGETS += gatk_validation
gatk_validation :
	$(call RUN_MAKE,modules/variant_callers/somatic/gatkValidation.mk)

TARGETS += samtools_validation
samtools_validation :
	$(call RUN_MAKE,modules/variant_callers/somatic/samtoolsValidation.mk)

TARGETS += norm_copynum
norm_copynum :
	$(call RUN_MAKE,modules/copy_number/normalisedCopyNum.mk)

TARGETS += mutation_summary
mutation_summary :
	$(call RUN_MAKE,modules/summary/mutationSummary.mk)

TARGETS += recurrent_mutations
recurrent_mutations :
	$(call RUN_MAKE,modules/recurrent_mutations/report.mk)

TARGETS += facets
facets :
	$(call RUN_MAKE,modules/copy_number/facets.mk)

TARGETS += brass
brass :
	$(call RUN_MAKE,modules/sv_callers/brass.mk)

TARGETS += mutsig_report
mutsig_report :
	$(call RUN_MAKE,modules/mut_sigs/mutSigReport.mk)

# standalone bam file merger
TARGETS += merge_bam
merge_bam :
	$(call RUN_MAKE,modules/bam_tools/mergeBam.mk)

# annotate external vcfs
TARGETS += ann_ext_vcf
ann_ext_vcf: 
	$(call RUN_MAKE,modules/vcf_tools/annotateExtVcf.mk)

# call hotspots
TARGETS += hotspot
hotspot: 
	$(call RUN_MAKE,modules/variant_callers/hotspot.mk)

TARGETS += strelka_varscan_indels
strelka_varscan_indels:
	$(call RUN_MAKE,modules/variant_callers/somatic/strelkaVarscanIndels.mk)

.PHONY : $(TARGETS)
