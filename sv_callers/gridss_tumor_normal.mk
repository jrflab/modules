include modules/Makefile.inc

LOGDIR = log/gridss_tumor_normal.$(NOW)

GRIDSS_CORES ?= 8
GRIDSS_MEM_CORE ?= 6G
GRIDSS_REF ?= $(HOME)/share/lib/ref_files/b37/human_g1k_v37.fasta
GRIDSS_BLACKLIST ?= $(HOME)/share/lib/resource_files/gridss/example/ENCFF001TDO.bed
GRIDSS ?= gridss
GRIDSS_FILTER ?= gridss_somatic_filter
GRIDSS_PON_DIR ?= $(HOME)/share/lib/resource_files/gridss/pon

gridss : $(foreach pair,$(SAMPLE_PAIRS),gridss/$(pair)/$(pair).gridss_sv.vcf)

define gridss-tumor-normal
gridss/$1_$2/$1_$2.gridss_sv.vcf : bam/$1.bam bam/$2.bam
	$$(call RUN,-c -n $(GRIDSS_CORES) -s 4G -m $(GRIDSS_MEM_CORE) -v $(GRIDSS_ENV) -w 72:00:00,"set -o pipefail && \
												    mkdir -p gridss/$1_$2 && \
												    cd gridss/$1_$2 && \
												    $$(GRIDSS) \
												    -t $$(GRIDSS_CORES) \
												    -r $$(GRIDSS_REF) \
												    -o $1_$2.gridss_sv.vcf \
												    -b $$(GRIDSS_BLACKLIST) \
												    ../../bam/$2.bam \
												    ../../bam/$1.bam")
												    
gridss/$1_$2/$1_$2.gridss_sv_ft.vcf : gridss/$1_$2/$1_$2.gridss_sv.vcf
	$$(call RUN,-c -n $(GRIDSS_CORES) -s 4G -m $(GRIDSS_MEM_CORE) -v $(GRIDSS_ENV),"set -o pipefail && \
											cd gridss/$1_$2 && \
											$$(GRIDSS_FILTER) \
											--pondir  $$(GRIDSS_PON_DIR) \
											--ref $$(GRIDSS_REF) \
											--input $1_$2.gridss.sv.vcf \
											--output $1_$2.gridss_sv_ft.vcf.gz \
											--fulloutput $1_$2.gridss_sv_high_and_low_confidence_somatic.vcf.gz \
											--scriptdir $$(GRIDSS_ENV)/bin \
											-n 1 \
											-t 2")
												    


#svaba/$1_$2.svaba.somatic.sv.vcf : svaba/$1_$2.svaba.somatic.indel.vcf

#svaba/$1_$2.svaba.unfiltered.somatic.sv.vcf : svaba/$1_$2.svaba.somatic.indel.vcf

#vcf/$1_$2.svaba_sv.vcf : svaba/$1_$2.svaba.somatic.sv.vcf
#	$$(INIT) cat $$< > $$@

endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call gridss-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))


..DUMMY := $(shell mkdir -p version; \
	     echo 'gridss' > version/gridss_tumor_normal.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: gridss
