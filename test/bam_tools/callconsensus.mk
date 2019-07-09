include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/call_consensus.$(NOW)
PHONY += fgbio

call_consensus : $(foreach sample,$(SAMPLES),fgbio/$(sample)-pool-A.duplex_qc.pdf) $(foreach sample,$(SAMPLES),fgbio/$(sample)-pool-B.duplex_qc.pdf) $(foreach sample,$(SAMPLES),fgbio/$(sample).filtered.bam)

TMPDIR ?= /home/${USER}/share/data/${USER}/tmp
JAVA = /home/${USER}/share/usr/jdk1.8.0_74/bin/java
PICARD = /home/${USER}/share/usr/picard/bin/picard.jar
REF_FASTA ?= /home/${USER}/share/reference/GATK_bundle/2.3/human_g1k_v37.fasta
POOL_A_INTERVAL ?= /home/${USER}/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.list
POOL_B_INTERVAL ?= /home/${USER}/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.list

define call-consensus
fgbio/%.groupedbyumi.bam : fgbio/%.regrouped.bam
	$$(call RUN,-c -n 1 -s 12G -m 24G -v $(FGBIO_ENV),"set -o pipefail && \
													   fgbio --tmp-dir $(TMPDIR) -Xms1G -Xmx24G GroupReadsByUmi --strategy Paired \
													   --input fgbio/$$(*).regrouped.bam \
													   --output fgbio/$$(*).groupedbyumi.bam \
													   --family-size-histogram fgbio/$$(*)-pool-AB-family_sizes.txt")

fgbio/%-pool-A.duplex_qc.pdf : fgbio/%.groupedbyumi.bam
	$$(call RUN,-c -n 1 -s 12G -m 24G -v $(FGBIO_ENV),"set -o pipefail && \
													   fgbio --tmp-dir $(TMPDIR) -Xms1G -Xmx24G CollectDuplexSeqMetrics \
													   --input fgbio/$$(*).groupedbyumi.bam \
													   --output fgbio/$$(*)-pool-A \
													   --intervals $(POOL_A_INTERVAL) \
													   --duplex-umi-counts true")

fgbio/%-pool-B.duplex_qc.pdf : fgbio/%.groupedbyumi.bam
	$$(call RUN,-c -n 1 -s 12G -m 24G -v $(FGBIO_ENV),"set -o pipefail && \
													   fgbio --tmp-dir $(TMPDIR) -Xms1G -Xmx24G CollectDuplexSeqMetrics \
													   --input fgbio/$$(*).groupedbyumi.bam \
													   --output fgbio/$$(*)-pool-B \
													   --intervals $(POOL_B_INTERVAL) \
													   --duplex-umi-counts true")
									  
fgbio/%.consensus.bam : fgbio/%.groupedbyumi.bam
	$$(call RUN,-c -n 1 -s 12G -m 24G -v $(FGBIO_ENV),"set -o pipefail && \
													   fgbio --tmp-dir $(TMPDIR) -Xms1G -Xmx24G CallDuplexConsensusReads \
													   --input fgbio/$$(*).groupedbyumi.bam \
													   --output fgbio/$$(*).consensus.bam \
													   --sort-order Queryname \
													   --min-reads 1 1 1")
													   
fgbio/%.filtered.bam : fgbio/%.consensus.bam
	$$(call RUN,-c -n 1 -s 12G -m 24G -v $(FGBIO_ENV),"set -o pipefail && \
													   fgbio --tmp-dir $TMPDIR -Xms1G -Xmx24G FilterConsensusReads \
													   --input fgbio/$$(*).consensus.bam \
													   --output fgbio/$$(*).filtered.bam \
													   --ref $(REF_FASTA) \
													   --min-reads 1 \
													   --max-read-error-rate 0.025 \
													   --max-base-error-rate 0.1 \
													   --min-base-quality 10")


endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call call-consensus,$(sample))))


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
