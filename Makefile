ifneq ("$(wildcard config.inc)", "")
	include config.inc
endif
ifneq ("$(wildcard project_config.inc)", "")
	include project_config.inc
endif
include modules/config.inc

export

NUM_ATTEMPTS ?= 10
NOW := $(shell date +"%F")
MAKELOG = log/$(@).$(NOW).log

USE_CLUSTER ?= true
QMAKE = modules/scripts/qmake.pl -n $@.$(NOW) $(if $(SLACK_CHANNEL),-c $(SLACK_CHANNEL)) -r $(NUM_ATTEMPTS) -m -s -- make
NUM_JOBS ?= 100

define RUN_QMAKE
$(QMAKE) -e -f $1 -j $2 $(TARGET) && \
	mkdir -p completed_tasks && \
	touch completed_tasks/$@
endef

RUN_MAKE = $(if $(findstring false,$(USE_CLUSTER))$(findstring n,$(MAKEFLAGS)),+$(MAKE) -f $1,$(call RUN_QMAKE,$1,$(NUM_JOBS)))

#==================================================
# workflows
#==================================================

TARGETS += somatic_indels
somatic_indels :
	$(call RUN_MAKE,modules/variant_callers/somatic/somaticIndels.mk)
	
TARGETS += somatic_variants
somatic_variants :
	$(call RUN_MAKE,modules/variant_callers/somatic/somaticVariants.mk)
	

#==================================================
# aligners
#==================================================

TARGETS += bwamem
bwamem :
	$(call RUN_MAKE,modules/aligners/bwamemAligner.mk)

TARGETS += bwa
bwa : NUM_ATTEMPTS = 50
bwa :
	$(call RUN_MAKE,modules/aligners/bwaAligner.mk)

TARGETS += bowtie
bowtie : NUM_ATTEMPTS = 50
bowtie :
	$(call RUN_MAKE,modules/aligners/bowtieAligner.mk)

TARGETS += tmap
tmap : NUM_ATTEMPTS = 50
tmap :
	$(call RUN_MAKE,modules/aligners/tmapAligner.mk)

TARGETS += hisat
hisat : 
	$(call RUN_MAKE,modules/aligners/hisatAligner.mk)

TARGETS += tophat
tophat : 
	$(call RUN_MAKE,modules/aligners/tophatAligner.mk)

TARGETS += star
star :
	$(call RUN_MAKE,modules/aligners/starAligner.mk)

TARGETS += star_fusion_aligner
star_fusion_aligner :
	$(call RUN_MAKE,modules/aligners/starFusionAligner.mk)
	
TARGETS += blast_reads
blast_reads :
	$(call RUN_MAKE,modules/fastq_tools/blastReads.mk)


#==================================================
# variant callers
#==================================================

TARGETS += mutect
mutect :
	$(call RUN_MAKE,modules/variant_callers/somatic/mutect.mk)
	
TARGETS += mutect2
mutect2 :
	$(call RUN_MAKE,modules/variant_callers/somatic/mutect2.mk)
	
TARGETS += somatic_sniper
somatic_sniper :
	$(call RUN_MAKE,modules/variant_callers/somatic/somaticSniper.mk)

TARGETS += snvmix
snvmix :
	$(call RUN_MAKE,modules/variant_callers/snvmix.mk)
	
TARGETS += tvcTN
tvcTN :
	$(call RUN_MAKE,modules/variant_callers/somatic/tvcTN.mk)

TARGETS += tvc
tvc :
	$(call RUN_MAKE,modules/variant_callers/tvc.mk)

TARGETS += varscanTN
varscanTN :
	$(call RUN_MAKE,modules/variant_callers/somatic/varscanTN.mk)
	
TARGETS += varscan
varscan :
	$(call RUN_MAKE,modules/variant_callers/varscan.mk)
	
TARGETS += strelka
strelka :
	$(call RUN_MAKE,modules/variant_callers/somatic/strelka.mk)

TARGETS += scalpel
scalpel :
	$(call RUN_MAKE,modules/variant_callers/somatic/scalpel.mk)
    
TARGETS += lancet
lancet :
	$(call RUN_MAKE,modules/variant_callers/somatic/lancet.mk)

TARGETS += pindelTN
pindelTN :
	$(call RUN_MAKE,modules/variant_callers/somatic/pindelTN.mk)
	
TARGETS += pindel
pindel :
	$(call RUN_MAKE,modules/variant_callers/pindel.mk)
	
TARGETS += gatk
gatk : 
	$(call RUN_MAKE,modules/variant_callers/gatk.mk)
	
TARGETS += haplotype_caller
haplotype_caller : 
	$(call RUN_MAKE,modules/variant_callers/haplotypeCaller.mk)
	
TARGETS += samtools_het
samtools_het :
	$(call RUN_MAKE,modules/variant_callers/samtoolsHet.mk)

TARGETS += platypus
platypus :
	$(call RUN_MAKE,modules/variant_callers/somatic/platypus.mk)
	
TARGETS += msisensor
msisensor :
	$(call RUN_MAKE,modules/variant_callers/somatic/msisensor.mk)	

TARGETS += hla_polysolver
hla_polysolver :
	$(call RUN_MAKE,modules/variant_callers/somatic/polysolver.mk)
	
TARGETS += pyrohmm
pyrohmm :
	$(call RUN_MAKE,modules/variant_callers/pyroHMMVar.mk)

TARGETS += museqTN
museqTN :
	$(call RUN_MAKE,modules/variant_callers/somatic/museqTN.mk)
	
TARGETS += hotspot
hotspot : 
	$(call RUN_MAKE,modules/variant_callers/hotspot.mk)
	
TARGETS += jsm
jsm :
	$(call RUN_MAKE,modules/variant_callers/somatic/jsm.mk)

TARGETS += sufam
sufam:
	$(call RUN_MAKE,modules/variant_callers/sufamsampleset.mk)

TARGETS += sufam_gt
sufam_gt :
	$(call RUN_MAKE,modules/variant_callers/sufam_gt.mk)

TARGETS += get_basecount
get_basecount :
	$(call RUN_MAKE,modules/variant_callers/get_basecounts.mk)
	
TARGETS += strelka_varscan_indels
strelka_varscan_indels :
	$(call RUN_MAKE,modules/variant_callers/somatic/strelkaVarscanIndels.mk)


#==================================================
# copy number
#==================================================

TARGETS += facets
facets :
	$(call RUN_MAKE,modules/copy_number/facets.mk)

TARGETS += facets_suite
facets_suite :
	$(call RUN_MAKE,modules/copy_number/facets_suite.mk)

TARGETS += ascat
ascat :
	$(call RUN_MAKE,modules/copy_number/ascat.mk)

TARGETS += norm_copynum
norm_copynum :
	$(call RUN_MAKE,modules/copy_number/normalisedCopyNum.mk)

TARGETS += titan
titan :
	$(call RUN_MAKE,modules/copy_number/titan.mk)

TARGETS += varscan_cnv
varscan_cnv :
	$(call RUN_MAKE,modules/copy_number/varscanCNV.mk)

TARGETS += hmmcopy
hmmcopy :
	$(call RUN_MAKE,modules/copy_number/hmmCopy.mk)

TARGETS += freec
freec : 
	$(call RUN_MAKE,modules/copy_number/controlFreeC.mk)

TARGETS += freecTN
freecTN : 
	$(call RUN_MAKE,modules/copy_number/controlFreeCTN.mk)

TARGETS += freec_lohTN
freec_lohTN : 
	$(call RUN_MAKE,modules/copy_number/controlFreeCLOHTN.mk)

TARGETS += exomecnv
exomecnv : 
	$(call RUN_MAKE,modules/copy_number/exomeCNV.mk)

TARGETS += exomecnvloh
exomecnvloh : 
	$(call RUN_MAKE,modules/copy_number/exomeCNVLOH.mk)

TARGETS += gistic
gistic :
	$(call RUN_MAKE,modules/copy_number/gistic.mk)
	
TARGETS += snp6
snp6 :
	$(call RUN_MAKE,modules/snp6/snp6.mk)
	
TARGETS += cnv_kit
cnv_kit :
	$(call RUN_MAKE,modules/copy_number/cnvkit.mk)


#==================================================
# RNAseq structural variant callers
#==================================================

TARGETS += star_fusion
star_fusion :
	$(call RUN_MAKE,modules/sv_callers/starFusion.mk)

TARGETS += tophat_fusion
tophat_fusion : 
	$(call RUN_MAKE,modules/sv_callers/tophatFusion.mk)

TARGETS += manta_rnaseq
manta_rnaseq :
	$(call RUN_MAKE,modules/sv_callers/mantaRnaseq.mk)
	
TARGETS += integrate_rnaseq
integrate_rnaseq :
	$(call RUN_MAKE,modules/sv_callers/integrateRnaseq.mk)

TARGETS += soapfuse
soapfuse :
	$(call RUN_MAKE,modules/sv_callers/soapFuse.mk)

TARGETS += mapsplice
mapsplice :
	$(call RUN_MAKE,modules/sv_callers/mapsplice.mk)

TARGETS += fusioncatcher
fusioncatcher :
	$(call RUN_MAKE,modules/sv_callers/fusioncatcher.mk)
	
TARGETS += oncofuse
oncofuse :
	$(call RUN_MAKE,modules/sv_callers/oncofuse.mk)


#==================================================
# DNA structural variant callers
#==================================================	

TARGETS += manta_tumor_normal
manta_tumor_normal :
	$(call RUN_MAKE,modules/sv_callers/manta_tumor_normal.mk)
	
TARGETS += svaba_tumor_normal
svaba_tumor_normal :
	$(call RUN_MAKE,modules/sv_callers/svaba_tumor_normal.mk)
	
TARGETS += gridss_tumor_normal
gridss_tumor_normal :
	$(call RUN_MAKE,modules/sv_callers/gridss_tumor_normal.mk)
	
TARGETS += manta
manta :
	$(call RUN_MAKE,modules/sv_callers/manta.mk)

TARGETS += brass
brass :
	$(call RUN_MAKE,modules/sv_callers/brass.mk)

TARGETS += integrate
integrate :
	$(call RUN_MAKE,modules/sv_callers/integrate.mk)

TARGETS += defuse
defuse :
	$(call RUN_MAKE,modules/sv_callers/defuse.mk)

NUM_CHIMSCAN_JOBS ?= 5
TARGETS += chimscan
chimscan :
	$(call RUN_MAKE_J,modules/sv_callers/chimerascan.mk,$(NUM_CHIMSCAN_JOBS))

TARGETS += lumpy
lumpy :
	$(call RUN_MAKE,modules/sv_callers/lumpy.mk)

TARGETS += hydra
hydra :
	$(call RUN_MAKE,modules/sv_callers/hydra.mk)

TARGETS += nfuse_wgss_wtss
nfuse_wgss_wtss :
	$(call RUN_MAKE,modules/sv_callers/nfuseWGSSWTSS.mk)

TARGETS += crest
crest :
	$(call RUN_MAKE,modules/sv_callers/crest.mk)

TARGETS += delly
delly :
	$(call RUN_MAKE,modules/sv_callers/delly.mk)
	

#==================================================
# BAM tools
#==================================================

TARGETS += fix_bam
fix_bam :
	$(call RUN_MAKE,modules/bam_tools/fix_bam.mk)

TARGETS += fix_rg
fix_rg :
	$(call RUN_MAKE,modules/bam_tools/fix_rg.mk)

TARGETS += fix_mate
fix_mate :
	$(call RUN_MAKE,modules/bam_tools/fix_mate.mk)

TARGETS += merge_bam
merge_bam :
	$(call RUN_MAKE,modules/bam_tools/merge_bam.mk)
	
TARGETS += process_bam
process_bam : 
	$(call RUN_MAKE,modules/bam_tools/processBam.mk)
	
TARGETS += getbam_irb_mirror
getbam_irb_mirror : 
	$(call RUN_MAKE,modules/bam_tools/getbam_irb_mirror.mk)
	
TARGETS += getbam_data_mirror
getbam_data_mirror : 
	$(call RUN_MAKE,modules/bam_tools/getbam_data_mirror.mk)
	
TARGETS += putbam_data_mirror
putbam_data_mirror : 
	$(call RUN_MAKE,modules/bam_tools/putbam_data_mirror.mk)

#==================================================
# VCF tools
#==================================================

TARGETS += merge_sv
merge_sv : 
	$(call RUN_MAKE,modules/vcf_tools/merge_sv.mk)
	
TARGETS += annotate_sv
annotate_sv : 
	$(call RUN_MAKE,modules/vcf_tools/annotate_sv.mk)
	

#==================================================
# FASTQ tools
#==================================================

TARGETS += merge_fastq
merge_fastq : 
	$(call RUN_MAKE,modules/fastq_tools/mergeFastq.mk)

TARGETS += merge_split_fastq
merge_split_fastq :
	$(call RUN_MAKE,modules/fastq_tools/mergeSplitFastq.mk)

TARGETS += extract_fastq
extract_fastq :
	$(call RUN_MAKE,modules/fastq_tools/extractFastq.mk)
	
TARGETS += extract_unmapped
extract_unmapped :
	$(call RUN_MAKE,modules/fastq_tools/extractReads.mk)
	
TARGETS += extract_unmapped_pairs
extract_unmapped_pairs :
	$(call RUN_MAKE,modules/fastq_tools/extractunmappedpairs.mk)

TARGETS += bam_to_fasta
bam_to_fasta :
	$(call RUN_MAKE,modules/fastq_tools/bamtoFasta.mk)


#==================================================
# QC
#==================================================

TARGETS += bam_metrics
bam_metrics :
	$(call RUN_MAKE,modules/qc/bam_metrics.mk)

TARGETS += bam_interval_metrics
bam_interval_metrics :
	$(call RUN_MAKE,modules/qc/bam_interval_metrics.mk)

TARGETS += wgs_metrics
wgs_metrics :
	$(call RUN_MAKE,modules/qc/wgs_metrics.mk)

TARGETS += rnaseq_metrics
rnaseq_metrics :
	$(call RUN_MAKE,modules/qc/rnaseqMetrics.mk)

TARGETS += fastqc
fastqc :
	$(call RUN_MAKE,modules/qc/fastqc.mk)

TARGETS += interval_qc
interval_qc :
	$(call RUN_MAKE,modules/qc/intervalBamQC.mk)

TARGETS += rseqc
rseqc :
	$(call RUN_MAKE,modules/qc/rseqc.mk)

TARGETS += qualimap
qualimap :
	$(call RUN_MAKE,modules/qc/qualimap.mk)

TARGETS += bam_stats
bam_stats :
	$(call RUN_MAKE,modules/qc/bamStats.mk)


#==================================================
# RNA sequencing
#==================================================

TARGETS += sum_reads
sum_reads :
	$(call RUN_MAKE,modules/rnaseq/sumreads.mk)

TARGETS += kallisto
kallisto :
	$(call RUN_MAKE,modules/rnaseq/kallisto.mk)
	
TARGETS += immune_deconv
immune_deconv :
	$(call RUN_MAKE,modules/rnaseq/immunedeconv.mk)
	

#==================================================
# Ploidy / Clonality
#==================================================

TARGETS += pyloh
pyloh :
	$(call RUN_MAKE,modules/ploidy/pyloh.mk)

TARGETS += clonehd
clonehd :
	$(call RUN_MAKE,modules/clonality/clonehd.mk)

TARGETS += absolute_seq
absolute_seq :
	$(call RUN_MAKE,modules/clonality/absoluteSeq.mk)
	
TARGETS += pyclone_13
pyclone_13 :
	$(call RUN_MAKE,modules/clonality/pyclone_13.mk)
	
TARGETS += pyclone_vi
pyclone_vi :
	$(call RUN_MAKE,modules/clonality/pyclone_vi.mk)

#==================================================
# mutational signatures
#==================================================

TARGETS += deconstruct_sigs
deconstruct_sigs :
	$(call RUN_MAKE,modules/signatures/deconstruct_sigs.mk)

TARGETS += sv_signature
sv_signature :
	$(call RUN_MAKE,modules/signatures/sv_signature.mk)
	
TARGETS += star_fish
star_fish :
	$(call RUN_MAKE,modules/signatures/star_fish.mk)
	
TARGETS += hr_detect
hr_detect :
	$(call RUN_MAKE,modules/signatures/hr_detect.mk)


#==================================================
# miscellaneous
#==================================================

TARGETS += cluster_samples
cluster_samples :
	$(call RUN_MAKE,modules/contamination/clusterSamples.mk)

TARGETS += contest
contest :
	$(call RUN_MAKE,modules/contamination/contest.mk)

TARGETS += virus_detection_bowtie2
virus_detection_bowtie2 :
	$(call RUN_MAKE,modules/virus/virus_detection_bowtie2.mk)
	
TARGETS += viral_detection
viral_detection :
	$(call RUN_MAKE,modules/test/workflows/viral_detection.mk)
	
TARGETS += krona_classify
krona_classify :
	$(call RUN_MAKE,modules/virus/krona_classify.mk)
	
TARGETS += medicc2
medicc2 :
	$(call RUN_MAKE,modules/copy_number/medicc2.mk)


#==================================================
# reports
#==================================================

TARGETS += recurrent_mutations
recurrent_mutations :
	$(call RUN_MAKE,modules/recurrent_mutations/report.mk)
	
TARGETS += genome_summary
genome_summary :
	$(call RUN_MAKE,modules/summary/genomesummary.mk)

TARGETS += mutation_summary
mutation_summary :
	$(call RUN_MAKE,modules/summary/mutationsummary.mk)
	
TARGETS += delmh_summary
delmh_summary :
	$(call RUN_MAKE,modules/summary/delmh_summary.mk)


#==================================================
# annotations
#==================================================

TARGETS += ann_ext_vcf
ann_ext_vcf : 
	$(call RUN_MAKE,modules/vcf_tools/annotateExtVcf.mk)

TARGETS += ann_somatic_vcf
ann_somatic_vcf : 
	$(call RUN_MAKE,modules/vcf_tools/annotateSomaticVcf.mk)

TARGETS += ann_vcf
ann_vcf : 
	$(call RUN_MAKE,modules/vcf_tools/annotateVcf.mk)
	
TARGETS += cravat_annotate
cravat_annotate :
	$(call RUN_MAKE,modules/vcf_tools/cravat_annotation.mk)
	
TARGETS += cravat_summary
cravat_summary :
	$(call RUN_MAKE,modules/summary/cravat_summary.mk)
	
TARGETS += ann_summary_vcf
ann_summary_vcf : 
	$(call RUN_MAKE,modules/vcf_tools/annotateSummaryVcf.mk)


#==================================================
# beta testing
#==================================================

TARGETS += hotspot_summary
hotspot_summary :
	$(MAKE) -f modules/variant_callers/genotypehotspots.mk -j $(NUM_JOBS)
	$(call RUN_MAKE,modules/summary/hotspotsummary.mk)
	
	
.PHONY : $(TARGETS)
