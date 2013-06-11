# Top-level Makefile
#
#
# Author: Raymond Lim <raylim@mm.st>
#

export

#REF = hg19
#include ~/share/modules/Makefile.inc

SAMPLE_FILE = samples.txt
#LANE_SAMPLE_FILE = lane_samples.txt
SAMPLE_PAIR_FILE = sample_pairs.txt # pairs tumor-normal

NUM_ATTEMPTS = 1
NOW := $(shell date +"%F")
QMAKE_BINARY = /common/sge/bin/lx-amd64/qmake
QMAKE = ~/share/scripts/qmake.pl -n $@.$(NOW) -r $(NUM_ATTEMPTS) -m -- $(QMAKE_BINARY)
QMAKEFLAGS = -cwd -V PATH -inherit -q jrf.q
FLAGS = -j 300

#SAMPLE_DIRS = $(HOME)/share/references/sample_dirs.txt
#FIND_LANES = ssh xhost08 sh $(HOME)/share/scripts/findLanes.sh $(SAMPLE_DIRS)

PHRED64 ?= false
HARD_FILTER_SNPS ?= true

.PHONY : gatk gsnap gsnap_iadb bwa jsm hmm_copy dindel merge_lanes process_bam get_bams miso qualimap bam_metrics museq bam_interval_metrics amplicon_qc defuse nfuse_wgss_wtss mutect fastqc somatic_sniper varscan


get_bams :
	$(MAKE) -f ~/share/modules/getBams.mk -j 50 >> $@.log

# create sample.lanes.txt [sample-id lane-id lane-bam-path]
%.lane.txt : %.txt
	$(FIND_LANES) $(SAMPLE_DIRS) < $< > $@

# retrieve lane bams
get_lanes :
	$(MAKE) -f ~/share/modules/getLaneBams.mk -j 50 >> $@.log 2>&1

# lane level bwa
lane_bwa : 
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/bwaAligner.mk SAMPLE_FILE=$(LANE_SAMPLE_FILE) $(FLAGS) $(TARGET)

# merge lanes
merge_lanes : 
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/mergeLanes.mk $(FLAGS) $(TARGET)

gatk : 
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e  -f ~/share/modules/gatkVariantCaller.mk $(FLAGS) $(TARGET)

# OPTIONS: HARD_FILTER_SNPS = true/false
# 		   POOL_SNP_RECAL = true/false
gsnap : 
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e  -k -f ~/share/modules/gsnapAligner.mk PHRED64=$(PHRED64) $(FLAGS) $(TARGET)

gsnap_iadb : 
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e  -k -f ~/share/modules/gsnapIADB.mk PHRED64=$(PHRED64) $(FLAGS) $(TARGET)


bwa : 
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/bwaAligner.mk $(FLAGS) $(TARGET)

novoalign_iadb : 
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/novoalignIADB.mk $(FLAGS) $(TARGET)

#process_bam : 
#$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/processBam.mk $(FLAGS) $(TARGET) >> $@.log

bam_metrics :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/bamMetrics.mk $(FLAGS) $(TARGET)

bam_interval_metrics :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/bamIntervalMetrics.mk $(FLAGS) $(TARGET)

fastqc :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/fastqc.mk $(FLAGS) $(TARGET)


# not tested on the cluster
# requires x11 for graphics
amplicon_qc :
	$(MAKE) -f ~/share/modules/ampliconBamQC.mk -j 50 >> $@.log

miso :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/miso.mk $(FLAGS) metrics && \
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/miso.mk $(FLAGS) summaries && \
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/miso.mk $(FLAGS) tar

jsm :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/jsm.mk $(FLAGS) $(TARGET)

mutect :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/mutect.mk $(FLAGS) $(TARGET)

varscan :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/varscanTN.mk $(FLAGS) $(TARGET)

# single sample mutation seq
museq :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/museq.mk $(FLAGS) $(TARGET)

somatic_sniper :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/somaticSniper.mk $(FLAGS) $(TARGET)


som_indels :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/somaticIndelDetector.mk $(FLAGS) $(TARGET)

read_depth :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/readDepth.mk $(FLAGS) $(TARGET)

qualimap :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/qualimap.mk $(FLAGS) $(TARGET)

hmm_copy :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/hmmCopy.mk $(FLAGS) $(TARGET)

nfuse_wgss_wtss :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/nfuseWGSSWTSS.mk $(FLAGS) $(TARGET)

sum_reads :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/sumRNASeqReads.mk $(FLAGS) $(TARGET) SAMPLE_FILE=${SAMPLE_FILE} >>

dindel :
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/dindel.mk $(FLAGS) windows && \
	$(QMAKE) $(QMAKEFLAGS) -N qmake.$@ -- -e -f ~/share/modules/dindel.mk $(FLAGS) vcf 2>&1

defuse :
	$(MAKE) -f ~/share/modules/defuse.mk -j 300 >> $@.log 2>&1 

cleanlinks :
	symlinks -dr .

clean :
	$(RM) tmp; \
	$(RM) gsc_bam; \
	$(RM) fastq; \
	$(RM) bwa/sai; \
	$(RM) gatk/chr_bam gatk/chr_vcf
