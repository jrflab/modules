# define the alignment postprocess steps
ifndef ALIGN_INC
PHRED64 ?= false
DUP_TYPE ?= markdup
NO_FILTER ?= false
NO_RECAL ?= false
NO_REALN ?= false
SPLIT_CHR ?= true
SPLIT_FASTQ ?= false
NO_SORT ?= false
FIX_RG ?= false # fix the readgroups so that it corresponds with the bam filename
SEQ_PLATFORM ?= illumina

BAM_SUFFIX := $(subst $( ),.,$(strip \
        $(if $(findstring true,$(FIX_RG)),rg)\
        $(if $(findstring false,$(NO_SORT)),sorted)\
        $(if $(findstring false,$(NO_FILTER)),filtered)\
        $(if $(findstring false,$(NO_REALN)),realn)\
        $(if $(findstring rmdup,$(DUP_TYPE)),rmdup)\
        $(if $(findstring markdup,$(DUP_TYPE)),markdup)\
        $(if $(findstring false,$(NO_RECAL)),recal)\
        bam))
endif
ALIGN_INC = true
