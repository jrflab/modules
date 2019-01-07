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
QMAKE = modules/scripts/qmake.pl -n $@.$(NOW) $(if $(SLACK_CHANNEL),-c $(SLACK_CHANNEL)) -r $(NUM_ATTEMPTS) -m -s -- make
NUM_JOBS ?= 50

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

TARGETS += copynumber_summary
copynumber_summary:
	$(MAKE) -f modules/copy_number/genomealtered.mk -j $(NUM_JOBS)
	$(MAKE) -f modules/copy_number/lstscore.mk -j $(NUM_JOBS)
	$(MAKE) -f modules/copy_number/ntaiscore.mk -j $(NUM_JOBS)
	$(MAKE) -f modules/copy_number/myriadhrdscore.mk -j $(NUM_JOBS)
	$(call RUN_MAKE,modules/summary/genomesummary.mk)
	
TARGETS += hotspot_summary
hotspot_summary:
	$(MAKE) -f modules/variant_callers/genotypehotspots.mk -j $(NUM_JOBS)
	$(call RUN_MAKE,modules/summary/hotspotsummary.mk)
	
TARGETS += viral_detection
viral_detection:
	$(MAKE) -f modules/fastq_tools/extractReads.mk -j $(NUM_JOBS)
	$(MAKE) -f modules/fastq_tools/bamtoFasta.mk -j $(NUM_JOBS)
	$(MAKE) -f modules/fastq_tools/blastReads.mk -j $(NUM_JOBS)
	$(call RUN_MAKE,modules/virus/kronaClassify.mk)
	
TARGETS += multisample_pyclone
multisample_pyclone:
	$(MAKE) -f modules/copy_number/ascat.mk -j $(NUM_JOBS)
	$(MAKE) -f modules/variant_callers/sufammultisample.mk -j $(NUM_JOBS)
	$(MAKE) -f modules/clonality/setuppyclone.mk -j $(NUM_JOBS)
	$(MAKE) -f modules/clonality/runpyclone.mk -j $(NUM_JOBS)
	$(call RUN_MAKE,modules/clonality/plotpyclone.mk)
	
TARGETS += run_cnvkit
run_cnvkit :
	$(MAKE) -f modules/copy_number/cnvkitcoverage.mk -j $(NUM_JOBS)
	$(MAKE) -f modules/copy_number/cnvkitreference.mk -j $(NUM_JOBS)
	$(MAKE) -f modules/copy_number/cnvkitfix.mk -j $(NUM_JOBS)
	$(call RUN_MAKE,modules/copy_number/cnvkitplot.mk)
	
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
	$(call RUN_MAKE,modules/variant_callers/sufamSampleSet.mk)
	
TARGETS += sufam_multisample
sufam_multisample:
	$(call RUN_MAKE,modules/variant_callers/sufammultisample.mk)


#==================================================
# copy number
#==================================================

TARGETS += facets
facets :
	$(call RUN_MAKE,modules/copy_number/facets.mk)
	
TARGETS += facets_plot
facets_plot :
	$(call RUN_MAKE,modules/copy_number/facetsplot.mk)
	
TARGETS += ascat
ascat :
	$(call RUN_MAKE,modules/copy_number/ascat.mk)

TARGETS += norm_copynum
norm_copynum :
	$(call RUN_MAKE,modules/copy_number/normalisedCopyNum.mk)

TARGETS += titan
titan :
	$(call RUN_MAKE,modules/copy_number/titan.mk)

TARGETS += strelka_varscan_indels
strelka_varscan_indels:
	$(call RUN_MAKE,modules/variant_callers/somatic/strelkaVarscanIndels.mk)

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
	
TARGETS += genome_altered
genome_altered :
	$(call RUN_MAKE,modules/copy_number/genomealtered.mk)
	
TARGETS += lst_score
lst_score :
	$(call RUN_MAKE,modules/copy_number/lstscore.mk)
	
TARGETS += ntai_score
ntai_score :
	$(call RUN_MAKE,modules/copy_number/ntaiscore.mk)
	
TARGETS += myriad_score
myriad_score :
	$(call RUN_MAKE,modules/copy_number/myriadhrdscore.mk)
	
TARGETS += snp6
snp6 :
	$(call RUN_MAKE,modules/snp6/snp6.mk)

TARGETS += cnvkit_coverage
cnvkit_coverage :
	$(call RUN_MAKE,modules/copy_number/cnvkitcoverage.mk)
	
TARGETS += cnvkit_reference
cnvkit_reference :
	$(call RUN_MAKE,modules/copy_number/cnvkitreference.mk)
	
TARGETS += cnvkit_fix
cnvkit_fix :
	$(call RUN_MAKE,modules/copy_number/cnvkitfix.mk)

TARGETS += cnvkit_plot
cnvkit_plot :
	$(call RUN_MAKE,modules/copy_number/cnvkitplot.mk)
	
TARGETS += cnvkit_heatmap
cnvkit_heatmap :
	$(call RUN_MAKE,modules/copy_number/cnvkitheatmap.mk)
	
TARGETS += cnvkit_pca
cnvkit_pca :
	$(call RUN_MAKE,modules/copy_number/cnvkitprcomp.mk)
	
TARGETS += cnvkit_qc
cnvkit_qc :
	$(call RUN_MAKE,modules/copy_number/cnvkitqc.mk)


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

TARGETS += check_bam
check_bam :
	$(call RUN_MAKE,modules/bam_tools/checkBam.mk)
	

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
	
TARGETS += setup_pyclone
setup_pyclone :
	$(call RUN_MAKE,modules/clonality/setuppyclone.mk)

TARGETS += run_pyclone
run_pyclone :
	$(call RUN_MAKE,modules/clonality/runpyclone.mk)
	
TARGETS += plot_pyclone
plot_pyclone :
	$(call RUN_MAKE,modules/clonality/plotpyclone.mk)

	
#==================================================
# mutational signatures
#==================================================

TARGETS += emu
emu :
	$(call RUN_MAKE,modules/mut_sigs/emu.mk)


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
	
TARGETS += krona_classify
krona_classify :
	$(call RUN_MAKE,modules/virus/krona_classify.mk)


#==================================================
# reports
#==================================================

TARGETS += recurrent_mutations
recurrent_mutations :
	$(call RUN_MAKE,modules/recurrent_mutations/report.mk)

TARGETS += mutsig_report
mutsig_report :
	$(call RUN_MAKE,modules/mut_sigs/mutSigReport.mk)
	
TARGETS += genome_summary
genome_summary :
	$(call RUN_MAKE,modules/summary/genomesummary.mk)

TARGETS += mutation_summary
mutation_summary :
	$(call RUN_MAKE,modules/summary/mutationSummary.mk)


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

	
#==================================================
# beta testing
#==================================================

TARGETS += sufam_multisample_test
sufam_multisample_test:
	$(call RUN_MAKE,modules/test/variant_callers/sufammultisample.mk)
	
TARGETS += qdnaseq_extract_test
qdnaseq_extract_test:
	$(call RUN_MAKE,modules/test/copy_number/qdnaseqextract.mk)
	
TARGETS += qdnaseq_copynumber_test
qdnaseq_copynumber_test:
	$(call RUN_MAKE,modules/test/copy_number/qdnaseqcopynumber.mk)
	
TARGETS += cnvkit_coverage_test
cnvkit_coverage_test :
	$(call RUN_MAKE,modules/test/copy_number/cnvkitcoverage.mk)
	
TARGETS += cnvkit_reference_test
cnvkit_reference_test :
	$(call RUN_MAKE,modules/test/copy_number/cnvkitreference.mk)
	
TARGETS += cnvkit_fix_test
cnvkit_fix_test :
	$(call RUN_MAKE,modules/test/copy_number/cnvkitfix.mk)

TARGETS += cnvkit_plot_test
cnvkit_plot_test :
	$(call RUN_MAKE,modules/test/copy_number/cnvkitplot.mk)
	
TARGETS += cravat_annotation
cravat_annotation :
	$(call RUN_MAKE,modules/test/annotations/cravat_annotation.mk)

TARGETS += cravat_summary
cravat_summary :
	$(call RUN_MAKE,modules/test/annotations/cravat_summary.mk)
	
#==================================================
# alpha testing
#==================================================
	
TARGETS += run_qdnaseq
run_qdnaseq :
	$(MAKE) -f modules/test/copy_number/qdnaseqextract.mk -j $(NUM_JOBS)
	$(call RUN_MAKE,modules/test/copy_number/qdnaseqcopynumber.mk)



.PHONY : $(TARGETS)
