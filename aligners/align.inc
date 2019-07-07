ifndef ALIGN_INC
BAM_PHRED64 ?= false
BAM_DUP_TYPE ?= markdup
BAM_NO_FILTER ?= false
BAM_NO_RECAL ?= false
BAM_NO_REALN ?= false
SPLIT_CHR ?= true
SPLIT_FASTQ ?= false
BAM_NO_SORT ?= false
BAM_FIX_RG ?= false
SEQ_PLATFORM ?= illumina

BAM_SUFFIX := $(subst $( ),.,$(strip \
        $(if $(findstring false,$(BAM_NO_SORT)),sorted)\
        $(if $(findstring false,$(BAM_NO_FILTER)),filtered)\
        $(if $(findstring true,$(PDX)),pdx_filtered)\
        $(if $(findstring true,$(BAM_FIX_RG)),rg)\
        $(if $(findstring false,$(BAM_NO_REALN)),realn)\
        $(if $(findstring rmdup,$(BAM_DUP_TYPE)),rmdup)\
        $(if $(findstring markdup,$(BAM_DUP_TYPE)),markdup)\
        $(if $(findstring false,$(BAM_NO_RECAL)),recal)\
        bam))
endif
ALIGN_INC = true
