# Top-level Makefile
#
#
# Author: Raymond Lim <raylim@mm.st>
#

export

NUM_ATTEMPTS ?= 3
NOW := $(shell date +"%F")
MAKELOG = log/$(@).$(NOW).log

#QMAKE_BINARY = /common/sge/bin/lx24-amd64/qmake
QMAKE_BINARY = /common/sge/bin/lx-amd64/qmake
QMAKE = ~/share/scripts/qmake.pl -n $@.$(NOW) -r $(NUM_ATTEMPTS) -m -- $(QMAKE_BINARY)
MAKE = ~/share/scripts/qmake.pl -n $@.$(NOW) -r $(NUM_ATTEMPTS) -m -- make
QMAKEFLAGS = -cwd -v -inherit -q jrf.q 
FLAGS = -j 75

#SAMPLE_DIRS = $(HOME)/share/references/sample_dirs.txt
#FIND_LANES = ssh xhost08 sh $(HOME)/share/scripts/findLanes.sh $(SAMPLE_DIRS)

PHRED64 ?= false
HARD_FILTER_SNPS ?= true


TARGETS = get_bams 
get_bams :
	$(MAKE) -f ~/share/modules/getBams.mk -j 50 >> $@.log

# create sample.lanes.txt [sample-id lane-id lane-bam-path]
%.lane.txt : %.txt
	$(FIND_LANES) $(SAMPLE_DIRS) < $< > $@

# retrieve lane bams
TARGETS += get_lanes
get_lanes :
	$(MAKE) -f ~/share/modules/getLaneBams.mk -j 50 >> $@.log 2>&1

# lane level bwa
TARGETS += lane_bwa
lane_bwa : 
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/bwaAligner.mk SAMPLE_FILE=$(LANE_SAMPLE_FILE) $(FLAGS) $(TARGET)

# merge lanes
TARGETS += merge
merge : 
	$(MAKE) -e -f ~/share/modules/merge.mk $(FLAGS) $(TARGET)
	
TARGETS += merge_fastq
merge_fastq : 
	$(MAKE) -e -f ~/share/modules/mergeFastq.mk $(FLAGS) $(TARGET)

TARGETS += gatk
gatk : 
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/gatkVariantCaller.mk $(FLAGS) $(TARGET)

# OPTIONS: HARD_FILTER_SNPS = true/false
# 		   POOL_SNP_RECAL = true/false
TARGETS += gsnap
gsnap : 
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e  -k -f ~/share/modules/gsnapAligner.mk PHRED64=$(PHRED64) $(FLAGS) $(TARGET)

TARGETS += gsnap_iadb
gsnap_iadb : 
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e  -k -f ~/share/modules/gsnapIADB.mk PHRED64=$(PHRED64) $(FLAGS) $(TARGET)

TARGETS += bwa
bwa : 
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/bwaAlignerMD5.mk $(FLAGS) $(TARGET)

#TARGETS += bwa_md5
#bwa_md5 : 
#$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/bwaAlignerMD5.mk $(FLAGS) $(TARGET)

TARGETS += bowtie
bowtie : 
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/bowtieAlignerMD5.mk $(FLAGS) $(TARGET)

TARGETS += tmap
tmap :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/tmapAligner.mk $(FLAGS) $(TARGET)

TARGETS += tophat_fusion
tophat_fusion : 
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/tophatFusion.mk $(FLAGS) $(TARGET)

TARGETS += novoalign_iadb
novoalign_iadb : 
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/novoalignIADB.mk $(FLAGS) $(TARGET)

TARGETS += process_bam
process_bam : 
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/processBamMD5.mk $(FLAGS) $(TARGET)

TARGETS += bam_metrics
bam_metrics :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/bamMetrics.mk $(FLAGS) $(TARGET)

TARGETS += bam_interval_metrics
bam_interval_metrics :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/bamIntervalMetrics.mk $(FLAGS) $(TARGET)

TARGETS += rnaseq_metrics
rnaseq_metrics :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/rnaseqMetrics.mk $(FLAGS) $(TARGET)

TARGETS += fastqc
fastqc :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/fastqc.mk $(FLAGS) $(TARGET)

# not tested on the cluster
# requires x11 for graphics
TARGETS += interval_qc
interval_qc :
	$(MAKE) -e -f ~/share/modules/intervalBamQC.mk -j 50 >> $@.log

TARGETS += miso
miso :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/miso.mk $(FLAGS) metrics && \
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/miso.mk $(FLAGS) summaries && \
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/miso.mk $(FLAGS) tar

TARGETS += jsm
jsm :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/jsm.mk $(FLAGS) $(TARGET)

TARGETS += mutect
mutect :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/mutect.mk $(FLAGS) $(TARGET)

TARGETS += varscan_cnv
varscan_cnv :
	$(MAKE) $(MAKEFLAGS) -f ~/share/modules/varscanTN.mk $(FLAGS) cnv

TARGETS += varscanTN
varscanTN :
	$(MAKE) $(MAKEFLAGS) -f ~/share/modules/varscanTN.mk $(FLAGS) vcfs tables

# single sample mutation seq
TARGETS += museqTN
museqTN :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/museqTN.mk $(FLAGS) $(TARGET)

TARGETS += merge_vcfTN
merge_vcfTN :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/vcfMergeTN.mk $(FLAGS) $(TARGET)

TARGETS += somatic_sniper
somatic_sniper :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/somaticSniper.mk $(FLAGS) $(TARGET)


TARGETS += som_indels
som_indels :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/somaticIndelDetector.mk $(FLAGS) $(TARGET)

TARGETS += snvmix
snvmix :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/snvmix.mk $(FLAGS) $(TARGET)

TARGETS += pyrohmm
pyrohmm :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/pyroHMMVar.mk $(FLAGS) $(TARGET)


TARGETS += compare_vcf
compare_vcf :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/vcfCompare.mk $(FLAGS) $(TARGET)

TARGETS += merge_vcf_platform
merge_vcf_platform :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/vcfMergePlatform.mk $(FLAGS) $(TARGET)

TARGETS += compare_vcfTN
compare_vcfTN :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/vcfCompareTN.mk $(FLAGS) $(TARGET)

TARGETS += read_depth
read_depth :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/readDepth.mk $(FLAGS) $(TARGET)

TARGETS += qualimap
qualimap :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/qualimap.mk $(FLAGS) $(TARGET)

TARGETS += hmm_copy
hmm_copy :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/hmmCopy.mk $(FLAGS) $(TARGET)

TARGETS += nfuse_wgss_wtss
nfuse_wgss_wtss :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/nfuseWGSSWTSS.mk $(FLAGS) $(TARGET)

TARGETS += sum_reads
sum_reads :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/sumRNASeqReads.mk $(FLAGS) $(TARGET)

TARGETS += dindel
dindel :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/dindel.mk $(FLAGS) windows && \
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/dindel.mk $(FLAGS) vcf 2>&1

TARGETS += exomecnv
exomecnv : 
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/exomeCNV.mk $(FLAGS) $(TARGET) 

TARGETS += exomecnvloh
exomecnvloh : 
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/exomeCNVLOH.mk $(FLAGS) $(TARGET) 


TARGETS += freec
freec : 
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/controlFreeC.mk $(FLAGS) $(TARGET) 

TARGETS += freecTN
freecTN : 
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/controlFreeCTN.mk $(FLAGS) $(TARGET) 

TARGETS += freec_lohTN
freec_lohTN : 
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/controlFreeCLOHTN.mk $(FLAGS) $(TARGET) 

TARGETS += defuse
defuse :
	$(MAKE) -e -f ~/share/modules/defuse.mk -j15 -k $(TARGET)

TARGETS += chimscan
chimscan :
	$(MAKE) -e -f ~/share/modules/chimerascan.mk $(FLAGS) $(TARGET)

TARGETS += pindel
pindel :
	$(MAKE) -e -f ~/share/modules/pindel.mk $(FLAGS) $(TARGET)

TARGETS += scalpel
scalpel :
	$(MAKE) -e -f ~/share/modules/scalpel.mk $(FLAGS) $(TARGET)

TARGETS += snp6
snp6 :
	$(MAKE) -e -f ~/share/modules/snp6.mk $(FLAGS) $(TARGET)

TARGETS += cleanlinks
cleanlinks :
	symlinks -dr .

TARGETS += clean
clean :
	$(RM) tmp; \
	$(RM) gsc_bam; \
	$(RM) fastq; \
	$(RM) bwa/sai; \
	$(RM) gatk/chr_bam gatk/chr_vcf

.PHONY : $(TARGETS)

