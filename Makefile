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
NUM_JOBS ?= 500

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
somatic_indels:
	$(call RUN_MAKE,modules/variant_callers/somatic/somaticIndels.mk)
	
TARGETS += somatic_variants
somatic_variants:
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
star:
	$(call RUN_MAKE,modules/aligners/starAligner.mk)

TARGETS += star_fusion_aligner
star_fusion_aligner:
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
tvcTN:
	$(call RUN_MAKE,modules/variant_callers/somatic/tvcTN.mk)

TARGETS += tvc
tvc:
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
platypus:
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
hotspot: 
	$(call RUN_MAKE,modules/variant_callers/hotspot.mk)
	
TARGETS += genotype_hotspot
genotype_hotspot:
	$(call RUN_MAKE,modules/variant_callers/genotypehotspots.mk)
	
TARGETS += genotype_pdx
genotype_pdx:
	$(call RUN_MAKE,modules/variant_callers/genotypepdx.mk)
	
TARGETS += jsm
jsm :
	$(call RUN_MAKE,modules/variant_callers/somatic/jsm.mk)

TARGETS += sufam
sufam:
	$(call RUN_MAKE,modules/variant_callers/sufamsampleset.mk)
	
TARGETS += get_basecount
get_basecount:
	$(call RUN_MAKE,modules/variant_callers/getBaseCount.mk)
	
TARGETS += strelka_varscan_indels
strelka_varscan_indels:
	$(call RUN_MAKE,modules/variant_callers/somatic/strelkaVarscanIndels.mk)


#==================================================
# copy number
#==================================================

TARGETS += facets
facets :
	$(call RUN_MAKE,modules/copy_number/facets.mk)
	
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
# structural variant callers
#==================================================

TARGETS += star_fusion
star_fusion:
	$(call RUN_MAKE,modules/sv_callers/starFusion.mk)

TARGETS += tophat_fusion
tophat_fusion : 
	$(call RUN_MAKE,modules/sv_callers/tophatFusion.mk)

TARGETS += manta_rnaseq
manta_rnaseq :
	$(call RUN_MAKE,modules/sv_callers/mantaRnaseq.mk)

TARGETS += manta
manta :
	$(call RUN_MAKE,modules/sv_callers/manta.mk)

TARGETS += mantaTN
mantaTN :
	$(call RUN_MAKE,modules/sv_callers/mantaTN.mk)
	
TARGETS += svabaTN
svabaTN :
	$(call RUN_MAKE,modules/sv_callers/svabaTN.mk)	

TARGETS += brass
brass :
	$(call RUN_MAKE,modules/sv_callers/brass.mk)

TARGETS += integrate_rnaseq
integrate_rnaseq :
	$(call RUN_MAKE,modules/sv_callers/integrateRnaseq.mk)

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

TARGETS += oncofuse
oncofuse :
	$(call RUN_MAKE,modules/sv_callers/oncofuse.mk)

TARGETS += lumpy
lumpy :
	$(call RUN_MAKE,modules/sv_callers/lumpy.mk)

TARGETS += hydra
hydra :
	$(call RUN_MAKE,modules/sv_callers/hydra.mk)

TARGETS += nfuse_wgss_wtss
nfuse_wgss_wtss :
	$(call RUN_MAKE,modules/sv_callers/nfuseWGSSWTSS.mk)

TARGETS += soapfuse
soapfuse :
	$(call RUN_MAKE,modules/sv_callers/soapFuse.mk)

TARGETS += mapsplice
mapsplice :
	$(call RUN_MAKE,modules/sv_callers/mapsplice.mk)

TARGETS += fusioncatcher
fusioncatcher :
	$(call RUN_MAKE,modules/sv_callers/fusioncatcher.mk)

TARGETS += crest
crest :
	$(call RUN_MAKE,modules/sv_callers/crest.mk)

TARGETS += delly
delly :
	$(call RUN_MAKE,modules/sv_callers/delly.mk)


#==================================================
# pre-processing
#==================================================

TARGETS += merge_fastq
merge_fastq : 
	$(call RUN_MAKE,modules/fastq_tools/mergeFastq.mk)

TARGETS += fix_bam
fix_bam :
	$(call RUN_MAKE,modules/bam_tools/fixBam.mk)

TARGETS += fix_rg
fix_rg :
	$(call RUN_MAKE,modules/bam_tools/fixRG.mk)

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
	
TARGETS += process_bam
process_bam : 
	$(call RUN_MAKE,modules/bam_tools/processBam.mk)

TARGETS += merge_bam
merge_bam :
	$(call RUN_MAKE,modules/bam_tools/mergeBam.mk)
	
TARGETS += split_rg
split_rg :
	$(call RUN_MAKE,modules/bam_tools/splitRG.mk)



#==================================================
# quality control
#==================================================

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
# rna sequencing
#==================================================

TARGETS += cufflinks
cufflinks : 
	$(call RUN_MAKE,modules/rnaseq/cufflinks.mk)

TARGETS += sum_reads
sum_reads :
	$(call RUN_MAKE,modules/rnaseq/sumRNASeqReads.mk)

TARGETS += exon_counts
exon_counts :
	$(call RUN_MAKE,modules/rnaseq/dexseq.mk)
	

#==================================================
# chip sequencing
#==================================================
	
TARGETS += macs2TN
macs2TN:
	$(call RUN_MAKE,modules/variant_callers/somatic/macs2TN.mk)


#==================================================
# ploidy
#==================================================

TARGETS += pyloh
pyloh :
	$(call RUN_MAKE,modules/ploidy/pyloh.mk)


#==================================================
# clonality
#==================================================

TARGETS += clonehd
clonehd :
	$(call RUN_MAKE,modules/clonality/clonehd.mk)

TARGETS += absolute_seq
absolute_seq :
	$(call RUN_MAKE,modules/clonality/absoluteSeq.mk)
	
TARGETS += ms_pyclone
ms_pyclone :
	$(call RUN_MAKE,modules/test/workflows/mspyclone.mk)
	
TARGETS += ss_pyclone
ss_pyclone :
	$(call RUN_MAKE,modules/test/workflows/pyclone.mk)
	

#==================================================
# mutational signatures
#==================================================

TARGETS += emu
emu :
	$(call RUN_MAKE,modules/signatures/emu.mk)
	
TARGETS += mut_sig
mut_sig :
	$(call RUN_MAKE,modules/signatures/mut_sig.mk)
	
TARGETS += deconstruct_sigs
deconstruct_sigs :
	$(call RUN_MAKE,modules/signatures/deconstruct_sigs.mk)


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
viral_detection:
	$(call RUN_MAKE,modules/test/workflows/viral_detection.mk)
	
TARGETS += krona_classify
krona_classify :
	$(call RUN_MAKE,modules/virus/krona_classify.mk)
	

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
	
TARGETS += cravat_summary
cravat_summary :
	$(call RUN_MAKE,modules/summary/cravat_summary.mk)
	
TARGETS += delmh_summary
delmh_summary :
	$(call RUN_MAKE,modules/summary/delmh_summary.mk)


#==================================================
# annotations
#==================================================

TARGETS += ann_ext_vcf
ann_ext_vcf: 
	$(call RUN_MAKE,modules/vcf_tools/annotateExtVcf.mk)

TARGETS += ann_somatic_vcf
ann_somatic_vcf: 
	$(call RUN_MAKE,modules/vcf_tools/annotateSomaticVcf.mk)

TARGETS += ann_vcf
ann_vcf: 
	$(call RUN_MAKE,modules/vcf_tools/annotateVcf.mk)
	
TARGETS += cravat_annotation
cravat_annotation :
	$(call RUN_MAKE,modules/test/workflows/cravat_annotation.mk)
	
TARGETS += cravat_annotate
cravat_annotate :
	$(call RUN_MAKE,modules/vcf_tools/cravat_annotation.mk)
	
TARGETS += ann_summary_vcf
ann_summary_vcf: 
	$(call RUN_MAKE,modules/vcf_tools/annotateSummaryVcf.mk)


#==================================================
# beta testing
#==================================================

TARGETS += hotspot_summary
hotspot_summary:
	$(MAKE) -f modules/variant_callers/genotypehotspots.mk -j $(NUM_JOBS)
	$(call RUN_MAKE,modules/summary/hotspotsummary.mk)
	

.PHONY : $(TARGETS)
