include modules/Makefile.inc
include modules/config.inc

LOGDIR = log/viral_detection.$(NOW)
PHONY += unmapped_reads

VIRUS_WORKFLOW += extract_unmapped
VIRUS_WORKFLOW += bam_to_fasta
VIRUS_WORKFLOW += blast_reads
VIRUS_WORKFLOW += krona_classify

viral_detection_workflow : $(VIRUS_WORKFLOW)

include modules/fastq_tools/extractReads.mk
include modules/fastq_tools/bamtoFasta.mk
include modules/fastq_tools/blastReads.mk
include modules/virus/krona_classify.mk

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)


