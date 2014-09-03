# run pyrohmmvar: realignment-based variant calling method for 454 and ion torrent

include ~/share/modules/Makefile.inc
include ~/share/modules/variant_callers/gatk.inc

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY: all

LOGDIR = log/qsnp.$(NOW)
QSNP = $(JAVA) $(JARDIR)/qsnp-1.0.jar

define QSNP_CONFIG
[general]\n\
chrFiles=$(FREEC_REF)\n\
chrLenFile=$(CHR_LEN)\n\
maxThreads=$(FREEC_THREADS)\n\
samtools=$(SAMTOOLS)\n\
outputDir=$3\n\
noisyData=$(NOISY_DATA)\n\
ploidy=2\n\
window=$(FREEC_WINDOW_SIZE)\n\
gemMappabilityFile=$(GEM_MAP_FILE)\n\
printNA=$(PRINT_NA)\n\
[sample]\n\
mateFile=$1\n\
inputFormat=BAM\n\
mateOrientation=FR\n\
[control]\n\
mateFile=$2\n\
inputFormat=BAM\n\
mateOrientation=FR\n\
[BAF]\n\
shiftInQuality=33\n\
SNPfile=$(SNP_TXT)\n\
$(FREEC_TARGET_CONFIG)
endef

