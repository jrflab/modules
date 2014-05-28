# Use ExomeCNV to detect copy number variants and LOH
# vim: set ft=make :

include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc


LOGDIR = log/control_freec.$(NOW)
FREEC = $(HOME)/share/usr/bin/freec
FREEC_THREADS = 4
# mem is per thread
FREEC_MEM = 4G
FREEC_HMEM = 6G

MAKE_GRAPH = $(HOME)/share/scripts/makeGraph.R

FREEC_WINDOW_SIZE = 10000


#GC_CONTENT_NORM = 0

VPATH ?= bam

ifeq ($(EXOME),true)
FREEC_TARGET_CONFIG =[target]\n\
captureRegions=$(EXOME_BED)
NOISY_DATA = false
PRINT_NA = false
else
ifdef TARGETS_FILE
FREEC_TARGET_CONFIG =[target]\n\
captureRegions=$(TARGETS_FILE)
NOISY_DATA = false
PRINT_NA = false
else
FREEC_TARGET_CONFIG = 
NOISY_DATA = false
PRINT_NA = true
endif
endif

#usage $(call FREEC_CONFIG,tumor-bam,normal-bam,output-dir)
define FREEC_CONFIG
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

PLOT_FREEC_COPY_NUM = $(RSCRIPT) $(HOME)/share/scripts/plotFreeCCopyNum.R
ANNOTATE_FREEC = $(RSCRIPT) $(HOME)/share/scripts/annotateFreeC.R
CBIND_CNV = $(RSCRIPT) $(HOME)/share/scripts/cbindCNVs.R

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all cnv config plots png tables

all : cnv config plots png tables
cnv : $(foreach pair,$(SAMPLE_PAIRS),freec/$(pair)/$(tumor.$(pair)).bam_ratio.txt.png)
#cnv : $(foreach i,$(SETS_SEQ),$(foreach tumor,$(call get_tumors,$(set.$i)),freec/$(tumor).bam_ratio.txt.png))
config : $(foreach pair,$(SAMPLE_PAIRS),freec/config/$(pair).config.txt)
png : freec/cnvs.png
tables : freec/recurrent_cnv.txt freec/annotated_cnv.txt

#$(call config-tumor-normal,tumor,normal)
define freec-config-tumor-normal
freec/config/$1_$2.config.txt : $1.bam $2.bam
	$$(INIT) echo -e "$$(call FREEC_CONFIG,$$<,$$(word 2,$$^),freec/$1_$2)" | sed 's/ //' >  $$@
endef
#$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call freec-config-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call freec-config-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))

define freec-tumor-normal
freec/$1_$2/$1.bam_ratio.txt : freec/config/$1_$2.config.txt
	$$(call LSCRIPT_PARALLEL_MEM,$$(FREEC_THREADS),$$(FREEC_MEM),$$(FREEC_HMEM),"$$(FREEC) -conf $$<")
	
freec/$1_$2/$1.bam_CNVs : freec/$1_$2/$1.bam_ratio.txt
endef
#$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call freec-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call freec-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))


freec/%.bam_ratio.txt.png : freec/%.bam_ratio.txt
	$(call LSCRIPT_MEM,2G,4G,"cat $(MAKE_GRAPH) | $(R) --slave --args 2 $<")

freec/annotated_cnv.txt : $(foreach pair,$(SAMPLE_PAIRS),freec/$(pair)/$(tumor.$(pair)).bam_CNVs)
	$(call LSCRIPT_MEM,2G,4G,"$(ANNOTATE_FREEC) --outDir $(@D) --txdb $(ENSEMBL_TXDB) --knownVariants $(KNOWN_CNVS) $^")

freec/cnvs.png : $(foreach pair,$(SAMPLE_PAIRS),freec/$(pair)/$(tumor.$(pair)).bam_ratio.txt)
	$(INIT) $(PLOT_FREEC_COPY_NUM) --outPrefix $(@:.png=) --centromereTable $(CENTROMERE_TABLE) $^

freec/recurrent_cnv.txt : $(foreach pair,$(SAMPLE_PAIRS),freec/$(pair)/$(tumor.$(pair)).bam_ratio.txt)
	$(INIT) $(CBIND_CNV) --ensemblTxdb $(ENSEMBL_TXDB) --outDir $(@D) $(^:ratio.txt=CNVs)
