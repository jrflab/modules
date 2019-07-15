include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/align_fastq.$(NOW)
PHONY += marianas

BWA_ALN_OPTS ?= -M
BWAMEM_REF_FASTA ?= $(REF_FASTA)
BWAMEM_THREADS = 12
BWAMEM_MEM_PER_THREAD = 2G

SAMTOOLS_THREADS = 12
SAMTOOLS_MEM = 12G
SAMTOOLS_MEM_THREAD = 2G

GATK_THREADS = 8
GATK_MEM_THREAD = 2G

align_fastq : $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample).bam)

define fastq-to-bam
marianas/$1/$1.bwamem.bam : marianas/$1/$1_R1_umi-clipped.fastq.gz marianas/$1/$1_R2_umi-clipped.fastq.gz
	$$(call RUN,-c -n $(BWAMEM_THREADS) -s 1G -m $(BWAMEM_MEM_PER_THREAD) -w 1440,"set -o pipefail && \
																		           $(BWA) mem -t $(BWAMEM_THREADS) $(BWA_ALN_OPTS) -R \"@RG\tID:$1\tLB:$1\tPL:${SEQ_PLATFORM}\tSM:$1\" $(BWAMEM_REF_FASTA) $$(^) | $(SAMTOOLS) view -bhS - > $$(@)")
																		           
marianas/$1/$1.sorted.bam : marianas/$1/$1.bwamem.bam
	$$(call RUN,-c -n $(SAMTOOLS_THREADS) -s 1G -m $(SAMTOOLS_MEM_THREAD) -w 1440,"set -o pipefail && \
									  									   		   samtools sort -@ $(SAMTOOLS_THREADS) -m $(SAMTOOLS_MEM) $$(^) -o $$(@) -T $(TMPDIR) && \
									  									   		   samtools index $$(@) && \
									  									   		   cp marianas/$1/$1.bam.bai marianas/$1/$1.bai")
									  									   		   
marianas/$1/$1.fixed.bam : marianas/$1/$1.sorted.bam
	$$(call RUN,-c -n 1 -s 12G -m 18G -w 1440,"set -o pipefail && \
									  		   /home/${USER}/share/usr/jdk1.8.0_74/bin/java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx16G -jar /home/${USER}/share/usr/picard/bin/picard.jar \
									  		   FixMateInformation \
									  		   I=$$(<) \
									  		   O=$$(@) \
									  		   SORT_ORDER=coordinate \
									  		   TMP_DIR=$(TMPDIR) \
									  		   COMPRESSION_LEVEL=0 \
									  		   CREATE_INDEX=true \
									  		   VALIDATION_STRINGENCY=LENIENT")

marianas/$1/$1.intervals : marianas/$1/$1.fixed.bam
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD) -w 1440,"set -o pipefail && \
									   							   		   /home/$(USER)/share/usr/jdk1.8.0_121/bin/java -Djava.io.tmpdir=$(TMPDIR) -Xms1G -Xmx12G -jar /home/$(USER)/share/usr/lib/java/GenomeAnalysisTK-3.7.jar \
									   							   		   -S LENIENT -T RealignerTargetCreator -I $$(^) -nt 8 -R $(REF_FASTA) -o $$(@) --known /home/$(USER)/share/reference/GATK_bundle/2.3/Mills_and_1000G_gold_standard.indels.b37.vcf.gz")

marianas/$1/$1.realn.bam : marianas/$1/$1.sorted.bam marianas/$1/$1.intervals
	$$(call RUN,-c -n $(GATK_THREADS) -s 1G -m $(GATK_MEM_THREAD) -w 1440,"set -o pipefail && \
									   							   		   /home/$(USER)/share/usr/jdk1.8.0_121/bin/java -Djava.io.tmpdir=$(TMPDIR) -Xms1G -Xmx12G -jar /home/$(USER)/share/usr/lib/java/GenomeAnalysisTK-3.7.jar \
									   							   		   -S LENIENT -T IndelRealigner -I $$(<) -R $(REF_FASTA) -targetIntervals $$(<<) -o $$(@) --knownAlleles /home/brownd7/share/reference/GATK_bundle/2.3/Mills_and_1000G_gold_standard.indels.b37.vcf.gz")
									   							   
marianas/$1/$1.dedup.bam : marianas/$1/$1.realn.bam
	$$(call RUN, -c -n 1 -s 12G -m 18G -w 1440,"java -Djava.io.tmpdir=$$(TMPDIR) -Xms2G -Xmx16G -jar $$(PICARD_JAR) MarkDuplicates \
												I=$$(<) \
												O=$$(@) \
												M=marianas/$1/$1.txt \
												TMP_DIR=$(TMPDIR)")

marianas/$1/$1.bam : marianas/$1/$1.dedup.bam
	$$(call RUN, -c -n 1 -s 12G -m 18G -w 1440,"java -Djava.io.tmpdir=$(TMPDIR) -Xms2G -Xmx16G -jar $$(PICARD_JAR) AddOrReplaceReadGroups \
												I=$$(<) \
												O=$$(@) \
												RGID=$1 \
												RGLB=$1 \
												RGPL=illumina \
												RGPU=NA \
												RGSM=$1 \
												TMP_DIR=$(TMPDIR) && \
												samtools index $$(@) && \
												cp marianas/$1/$1.bam.bai marianas/$1/$1.bai")


endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call fastq-to-bam,$(sample))))


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)


