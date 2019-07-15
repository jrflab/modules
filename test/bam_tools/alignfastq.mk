include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/align_fastq.$(NOW)
PHONY += marianas

BWA_ALN_OPTS ?= -M
BWAMEM_REF_FASTA ?= $(REF_FASTA)
BWAMEM_THREADS = 12
BWAMEM_MEM_PER_THREAD = 2G

align_fastq : $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample).intervals)

define fastq-to-bam
marianas/$1/$1.bwamem.bam : marianas/$1/$1_R1_umi-clipped.fastq.gz marianas/$1/$1_R2_umi-clipped.fastq.gz
	$$(call RUN,-c -n $(BWAMEM_THREADS) -s 1G -m $(BWAMEM_MEM_PER_THREAD),"set -o pipefail && \
																		   $(BWA) mem -t $(BWAMEM_THREADS) $(BWA_ALN_OPTS) -R \"@RG\tID:$1\tLB:$1\tPL:${SEQ_PLATFORM}\tSM:$1\" $(BWAMEM_REF_FASTA) $$(^) | $(SAMTOOLS) view -bhS - > $$(@)")

marianas/$1/$1.sorted.bam : marianas/$1/$1.bwamem.bam
	$$(call RUN,-c -n 12 -s 1G -m 2G,"set -o pipefail && \
									  samtools sort -@ 12 -m 12G $$(^) -o $$(@) -T $(TMPDIR) && \
									  samtools index $$(@) && \
									  cp marianas/$1/$1.bam.bai marianas/$1/$1.bai")

marianas/$1/$1.intervals : marianas/$1/$1.sorted.bam
	$$(call RUN,-c -n 8 -s 1G -m 2G,"set -o pipefail && \
									   /home/$(USER)/share/usr/jdk1.8.0_121/bin/java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx12G -jar /home/$(USER)/share/usr/lib/java/GenomeAnalysisTK-3.7.jar \
									   -S LENIENT -T RealignerTargetCreator -I $$(^) -nt 8 -R $(REF_FASTA) -o $$(@) --known /home/$(USER)/share/reference/GATK_bundle/2.3/Mills_and_1000G_gold_standard.indels.b37.vcf.gz")


endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call fastq-to-bam,$(sample))))


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
