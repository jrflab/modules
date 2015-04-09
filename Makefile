# Top-level Makefile
#
#
# Author: Raymond Lim <raylim@mm.st>
#

export

NUM_ATTEMPTS ?= 20
NOW := $(shell date +"%F")
MAKELOG = log/$(@).$(NOW).log

QMAKE = scripts/qmake.pl -n $@.$(NOW) -r $(NUM_ATTEMPTS) -m -- $(QMAKE_BINARY)
MAKE = scripts/qmake.pl -n $@.$(NOW) -r $(NUM_ATTEMPTS) -m -- make
FLAGS ?= -j 100

#SAMPLE_DIRS = $(HOME)/share/references/sample_dirs.txt
#FIND_LANES = ssh xhost08 sh $(HOME)/share/scripts/findLanes.sh $(SAMPLE_DIRS)

PHRED64 ?= false
HARD_FILTER_SNPS ?= true


TARGETS = get_bams 
get_bams :
	$(MAKE) -f modules/getBams.mk -j 50 >> $@.log

# create sample.lanes.txt [sample-id lane-id lane-bam-path]
%.lane.txt : %.txt
	$(FIND_LANES) $(SAMPLE_DIRS) < $< > $@

# retrieve lane bams
TARGETS += get_lanes
get_lanes :
	$(MAKE) -f modules/getLaneBams.mk -j 50 >> $@.log 2>&1


TARGETS += merge_fastq
merge_fastq : 
	$(MAKE) -e -f modules/fastq_tools/mergeFastq.mk $(FLAGS) $(TARGET)

TARGETS += merge_split_fastq
merge_split_fastq : 
	$(MAKE) -e -f modules/fastq_tools/fastq.mk MERGE_SPLIT_FASTQ=true $(FLAGS) $(TARGET)

TARGETS += gatk
gatk : 
	$(MAKE) $(MAKEFLAGS) -e -f modules/variant_callers/gatkVariantCaller.mk $(FLAGS) $(TARGET)

TARGETS += bwa
bwa : NUM_ATTEMPTS = 50
bwa :
	$(MAKE) $(MAKEFLAGS) -e -f modules/aligners/bwaAligner.mk $(FLAGS) $(TARGET)

TARGETS += bwamem
bwamem :
	$(MAKE) $(MAKEFLAGS) -e -f modules/aligners/bwamemAligner.mk $(FLAGS) $(TARGET)


TARGETS += bowtie
bowtie : NUM_ATTEMPTS = 50
bowtie :
	$(MAKE) $(MAKEFLAGS) -e -f modules/aligners/bowtieAligner.mk $(FLAGS) $(TARGET)

TARGETS += tmap
tmap : NUM_ATTEMPTS = 50
tmap :
	$(MAKE) $(MAKEFLAGS) -e -f modules/aligners/tmapAligner.mk $(FLAGS) $(TARGET)

TARGETS += tophat_fusion
tophat_fusion : 
	$(MAKE) $(MAKEFLAGS) -e -f modules/sv_callers/tophatFusion.mk $(FLAGS) $(TARGET)

TARGETS += tophat
tophat : 
	$(MAKE) $(MAKEFLAGS) -e -f modules/aligners/tophatAligner.mk $(FLAGS) $(TARGET)

TARGETS += cufflinks
cufflinks : 
	$(MAKE) $(MAKEFLAGS) -e -f modules/rnaseq/cufflinks.mk $(FLAGS) $(TARGET)

TARGETS += process_bam
process_bam : 
	$(MAKE) $(MAKEFLAGS) -e -f modules/bam_tools/processBam.mk $(FLAGS) $(TARGET)

TARGETS += bam_metrics
bam_metrics :
	$(MAKE) $(MAKEFLAGS) -e -f modules/qc/bamMetrics.mk $(FLAGS) $(TARGET)

TARGETS += bam_interval_metrics
bam_interval_metrics :
	$(MAKE) $(MAKEFLAGS) -e -f modules/qc/bamIntervalMetrics.mk $(FLAGS) $(TARGET)

TARGETS += rnaseq_metrics
rnaseq_metrics :
	$(MAKE) $(MAKEFLAGS) -e -f modules/qc/rnaseqMetrics.mk $(FLAGS) $(TARGET)

TARGETS += fastqc
fastqc :
	$(MAKE) $(MAKEFLAGS) -e -f modules/qc/fastqc.mk $(FLAGS) $(TARGET)

# not tested on the cluster
# requires x11 for graphics
TARGETS += interval_qc
interval_qc :
	$(MAKE) -e -f modules/qc/intervalBamQC.mk -j 50 >> $@.log

TARGETS += miso
miso :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f modules/isoforms/miso.mk $(FLAGS) metrics && \
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f modules/isoforms/miso.mk $(FLAGS) summaries && \
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f modules/isoforms/miso.mk $(FLAGS) tar

TARGETS += jsm
jsm :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f modules/variant_callers/somatic/jsm.mk $(FLAGS) $(TARGET)

TARGETS += mutect
mutect :
	$(MAKE) $(MAKEFLAGS) -e -f modules/variant_callers/somatic/mutectVariantCaller.mk $(FLAGS) $(TARGET)

TARGETS += varscan_cnv
varscan_cnv :
	$(MAKE) $(MAKEFLAGS) -e -f modules/copy_number/varscanCNV.mk $(FLAGS) $(TARGET)

TARGETS += varscan_fpfilter
varscan_fpfilter :
	$(MAKE) $(MAKEFLAGS) -e -f modules/variant_callers/varscanFpfilter.mk $(FLAGS) $(TARGET)

TARGETS += varscanTN
varscanTN :
	$(MAKE) $(MAKEFLAGS) -f modules/variant_callers/somatic/varscanTN.mk $(FLAGS) $(TARGET)

TARGETS += varscan
varscan :
	$(MAKE) $(MAKEFLAGS) -f modules/variant_callers/varscan.mk $(FLAGS)


# single sample mutation seq
TARGETS += museqTN
museqTN :
	$(MAKE) $(MAKEFLAGS) -e -f modules/variant_callers/somatic/museqTN.mk $(FLAGS) $(TARGET)

TARGETS += merge_vcfTN
merge_vcfTN :
	$(MAKE) $(MAKEFLAGS) -e -f modules/vcf_tools/vcfMergeTN.mk $(FLAGS) $(TARGET)

TARGETS += somatic_sniper
somatic_sniper :
	$(MAKE) $(MAKEFLAGS) -e -f modules/variant_callers/somatic/somaticSniper.mk $(FLAGS) $(TARGET)


TARGETS += som_indels
som_indels :
	$(MAKE) $(MAKEFLAGS) -e -f modules/variant_callers/somatic/somaticIndelDetector.mk $(FLAGS) $(TARGET)

TARGETS += snvmix
snvmix :
	$(MAKE) $(MAKEFLAGS) -e -f modules/variant_callers/snvmix.mk $(FLAGS) $(TARGET)

TARGETS += pyrohmm
pyrohmm :
	$(MAKE) $(MAKEFLAGS) -e -f modules/variant_callers/pyroHMMVar.mk $(FLAGS) $(TARGET)


TARGETS += compare_vcf
compare_vcf :
	$(MAKE) $(MAKEFLAGS) -e -f modules/vcf_tools/vcfCompare.mk $(FLAGS) $(TARGET)

TARGETS += merge_vcf_platform
merge_vcf_platform :
	$(MAKE) $(MAKEFLAGS) -e -f modules/vcf_tools/vcfMergePlatform.mk $(FLAGS) $(TARGET)

TARGETS += compare_vcfTN
compare_vcfTN :
	$(MAKE) $(MAKEFLAGS) -e -f modules/vcf_tools/vcfCompareTN.mk $(FLAGS) $(TARGET)

TARGETS += read_depth
read_depth :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f modules/qc/readDepth.mk $(FLAGS) $(TARGET)

TARGETS += qualimap
qualimap :
	$(MAKE) $(MAKEFLAGS) -e -f modules/qc/qualimap.mk $(FLAGS) $(TARGET)

TARGETS += hmmcopy
hmmcopy :
	$(MAKE) $(MAKEFLAGS) -e -f modules/copy_number/hmmCopy.mk $(FLAGS) $(TARGET)

TARGETS += nfuse_wgss_wtss
nfuse_wgss_wtss :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f modules/copy_number/nfuseWGSSWTSS.mk $(FLAGS) $(TARGET)

TARGETS += sum_reads
sum_reads :
	$(MAKE) $(MAKEFLAGS) -e -f modules/rnaseq/sumRNASeqReads.mk $(FLAGS) $(TARGET)

TARGETS += dindel
dindel :
	$(MAKE) $(MAKEFLAGS) -e -f modules/variant_callers/dindel.mk $(FLAGS) windows && \
	$(MAKE) $(MAKEFLAGS) -e -f modules/variant_callers/dindel.mk $(FLAGS) vcf 2>&1

TARGETS += exomecnv
exomecnv : 
	$(MAKE) $(MAKEFLAGS) -e -f modules/copy_number/exomeCNV.mk $(FLAGS) $(TARGET) 

TARGETS += exomecnvloh
exomecnvloh : 
	$(MAKE) $(MAKEFLAGS) -e -f modules/copy_number/exomeCNVLOH.mk $(FLAGS) $(TARGET) 

TARGETS += gistic
gistic :
	$(MAKE) $(MAKEFLAGS) -e -f modules/copy_number/gistic.mk $(FLAGS) $(TARGET) 

TARGETS += freec
freec : 
	$(MAKE) $(MAKEFLAGS) -e -f modules/copy_number/controlFreeC.mk $(FLAGS) $(TARGET) 

TARGETS += freecTN
freecTN : 
	$(MAKE) $(MAKEFLAGS) -e -f modules/copy_number/controlFreeCTN.mk $(FLAGS) $(TARGET) 

TARGETS += freec_lohTN
freec_lohTN : 
	$(MAKE) $(MAKEFLAGS) -e -f modules/copy_number/controlFreeCLOHTN.mk $(FLAGS) $(TARGET) 

NUM_DEFUSE_JOBS ?= 5
TARGETS += defuse
defuse :
	$(MAKE) -e -f modules/sv_callers/defuse.mk -j$(NUM_DEFUSE_JOBS) -k $(TARGET)

TARGETS += oncofuse
oncofuse :
	$(MAKE) -e -f modules/sv_callers/oncofuse.mk $(FLAGS) $(TARGET)

TARGETS += lumpy
lumpy :
	$(MAKE) -e -f modules/sv_callers/lumpy.mk -j5 -k $(TARGET)

TARGETS += hydra
hydra :
	$(MAKE) -e -f modules/sv_callers/hydra.mk $(FLAGS) $(TARGET)

NUM_CHIMSCAN_JOBS ?= 5
TARGETS += chimscan
chimscan :
	$(MAKE) -e -f modules/sv_callers/chimerascan.mk -j$(NUM_CHIMSCAN_JOBS) $(TARGET)

TARGETS += pindel
pindel :
	$(MAKE) -e -f modules/variant_callers/pindel.mk $(FLAGS) $(TARGET)

TARGETS += scalpel
scalpel :
	$(MAKE) -e -f modules/variant_callers/somatic/scalpel.mk $(FLAGS) $(TARGET)

TARGETS += snp6
snp6 :
	$(MAKE) -e -f modules/snp6/snp6.mk $(FLAGS) $(TARGET)

TARGETS += snpcaller
snpcaller :
	$(MAKE) -e -f modules/variant_callers/snpCaller.mk $(FLAGS) $(TARGET)

TARGETS += soapfuse
soapfuse :
	$(MAKE) -e -f modules/sv_callers/soapFuse.mk $(FLAGS) $(TARGET)

TARGETS += mapsplice
mapsplice :
	$(MAKE) -e -f modules/sv_callers/mapsplice.mk $(FLAGS) $(TARGET)

TARGETS += fusioncatcher
fusioncatcher :
	$(MAKE) -e -f modules/sv_callers/fusioncatcher.mk $(FLAGS) $(TARGET)

TARGETS += strelka
strelka :
	$(MAKE) -e -f modules/variant_callers/somatic/strelka.mk $(FLAGS) $(TARGET)

TARGETS += crest
crest :
	$(MAKE) -e -f modules/variant_callers/somatic/crest.mk $(FLAGS) $(TARGET)

TARGETS += pyloh
pyloh :
	$(MAKE) -e -f modules/ploidy/pyloh.mk $(FLAGS) $(TARGET)

TARGETS += clonehd
clonehd :
	$(MAKE) -e -f modules/clonality/clonehd.mk $(FLAGS) $(TARGET)

TARGETS += emu
emu :
	$(MAKE) -e -f modules/mut_sigs/emu.mk $(FLAGS) $(TARGET)


TARGETS += extract_fastq
extract_fastq :
	$(MAKE) -e -f modules/fastq_tools/extractFastq.mk $(FLAGS) $(TARGET)

TARGETS += titan
titan :
	$(MAKE) -e -f modules/copy_number/titan.mk $(FLAGS) $(TARGET)

TARGETS += samtools_het
samtools_het :
	$(MAKE) -e -f modules/variant_callers/samtoolsHet.mk $(FLAGS) $(TARGET)

TARGETS += absolute_seq
absolute_seq :
	$(MAKE) -e -f modules/clonality/absoluteSeq.mk $(FLAGS) $(TARGET)

TARGETS += merge_strelka_varscan
merge_strelka_varscan :
	$(MAKE) -e -f modules/variant_callers/somatic/mergeStrelkaVarscanIndels.mk $(FLAGS) $(TARGET)

TARGETS += clean
clean :
	$(RM) tmp; \
	$(RM) gsc_bam; \
	$(RM) fastq; \
	$(RM) bwa/sai; \
	$(RM) gatk/chr_bam gatk/chr_vcf

.PHONY : $(TARGETS)

