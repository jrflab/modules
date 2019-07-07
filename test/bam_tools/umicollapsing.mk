include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/umicollapsing.$(NOW)
PHONY += fgbio

collapsed_umi : $(foreach sample,$(SAMPLES),fgbio/$(sample).ubam)

TMPDIR ?= /home/${USER}/share/data/${USER}/tmp
JAVA = /home/${USER}/share/usr/jdk1.8.0_74/bin/java
PICARD = /home/${USER}/share/usr/picard/bin/picard.jar
REF_FASTA ?= /home/${USER}/share/reference/GATK_bundle/2.3/human_g1k_v37.fasta
POOL_A_INTERVAL ?= /home/${USER}/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.list
POOL_B_INTERVAL ?= /home/${USER}/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.list

define fastq-to-ubam
fgbio/%.ubam : 
	$$(call RUN,-c -n 1 -s 8G -m 16G -v $(FGBIO_ENV),"set -o pipefail && \
													  fgbio --tmp-dir $(TMPDIR) -Xms1g -Xmx12g FastqToBam --input CD10_P_Post_IGO_09342_B_8_S128_R1_001.fastq.gz CD10_P_Post_IGO_09342_B_8_S128_R2_001.fastq.gz \
													  --read-structures 3M+T 3M+T \
													  --sample CD10-P-Post \
													  --output CD10-P-Post.ubam \
													  --library CD10-P-Post && \
													  $JAVA -Xmx8G -jar $PICARD SortSam \
													  I=CD10-P-Post.ubam \
													  O=CD10-P-Post.qn.sorted.ubam \
													  SORT_ORDER=queryname \
													  TMP_DIR=$TMPDIR")
endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call fix-bam,$(sample))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
