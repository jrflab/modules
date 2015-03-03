# Use ExomeCNV to detect copy number variants and LOH
# vim: set ft=make :

include modules/Makefile.inc
include modules/variant_callers/gatk.inc


LOGDIR = log/control_freec_loh.$(NOW)
FREEC = $(HOME)/share/usr/bin/freec
FREEC_THREADS = 4
# mem is per thread
FREEC_MEM = 4G
FREEC_HMEM = 6G

MAKE_GRAPH = scripts/makeGraph.R

FREEC_WINDOW_SIZE = 10000


#GC_CONTENT_NORM = 0

ifeq ($(EXOME),true)
FREEC_TARGET_CONFIG =[target]\n\
captureRegions=$(EXOME_BED)
NOISY_DATA = true
PRINT_NA = false
else
ifdef TARGETS_FILE
FREEC_TARGET_CONFIG =[target]\n\
captureRegions=$(TARGETS_FILE)
NOISY_DATA = true
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
coefficientOfVariation=0.05\n\
window=$(FREEC_WINDOW_SIZE)\n\
gemMappabilityFile=$(GEM_MAP_FILE)\n\
printNA=$(PRINT_NA)\n\
[sample]\n\
mateFile=$1\n\
inputFormat=pileup\n\
mateOrientation=FR\n\
[control]\n\
mateFile=$2\n\
inputFormat=pileup\n\
mateOrientation=FR\n\
[BAF]\n\
shiftInQuality=33\n\
SNPfile=$(SNP_TXT)\n\
$(FREEC_TARGET_CONFIG)
endef

ifdef TARGETS_FILE
PILEUP_OPTS += -l $(TARGETS_FILE)
endif

PLOT_FREEC_COPY_NUM = $(RSCRIPT) scripts/plotFreeCCopyNum.R
CBIND_CNV = $(RSCRIPT) scripts/cbindCNVs.R

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all cnv config plots png tables

all : cnv config plots png tables
cnv : $(foreach i,$(SETS_SEQ),$(foreach tumor,$(call get_tumors,$(set.$i)),freec/$(tumor).pileup_ratio.txt.png))
config : $(foreach pair,$(SAMPLE_PAIRS),freec/$(pair).config.txt)
png : freec/cnvs.png
tables : freec/recurrent_cnv.txt

pileup/%.pileup.gz : bam/%.bam
	$(call LSCRIPT,"$(SAMTOOLS) mpileup $(PILEUP_OPTS) -f $(REF_FASTA) $< | sed '/^MT/d' | gzip -c > $@")

#$(call config-tumor-normal,tumor,normal)
define freec-config-tumor-normal
freec/$1_$2.config.txt : pileup/$1.pileup.gz pileup/$2.pileup.gz
	$$(INIT) echo -e "$$(call FREEC_CONFIG,$$<,$$(word 2,$$^),$$(@D))" | sed 's/ //' >  $$@
endef
#$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call freec-config-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call freec-config-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))

define freec-tumor-normal
freec/$1.pileup_ratio.txt : freec/$1_$2.config.txt
	$$(call LSCRIPT_PARALLEL_MEM,$$(FREEC_THREADS),$$(FREEC_MEM),$$(FREEC_HMEM),"$$(FREEC) -conf $$< &> $$(LOG)")
endef
#$(foreach tumor,$(TUMOR_SAMPLES),$(eval $(call freec-tumor-normal,$(tumor),$(normal_lookup.$(tumor)))))
$(foreach i,$(SETS_SEQ),\
	$(foreach tumor,$(call get_tumors,$(set.$i)), \
		$(eval $(call freec-tumor-normal,$(tumor),$(call get_normal,$(set.$i))))))

freec/%.pileup_ratio.txt.png : freec/%.pileup_ratio.txt
	$(call LSCRIPT_MEM,2G,4G,"cat $(MAKE_GRAPH) | $(R) --slave --args 2 $< &> $(LOG)")

freec/cnvs.png : $(foreach i,$(SETS_SEQ),$(foreach tumor,$(call get_tumors,$(set.$i)),freec/$(tumor).pileup_ratio.txt))
	$(INIT) $(PLOT_FREEC_COPY_NUM) --outFile $@ --centromereTable $(CENTROMERE_TABLE) $^

freec/recurrent_cnv.txt : $(foreach tumor,$(TUMOR_SAMPLES),freec/$(tumor).pileup_ratio.txt)
	$(INIT) $(CBIND_CNV) --ensemblTxdb $(ENSEMBL_TXDB) --outDir $(@D) $(^:ratio.txt=CNVs)
