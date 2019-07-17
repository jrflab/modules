include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/collapseumi.$(NOW)
PHONY += marianas

umi_collapse : $(foreach sample,$(SAMPLES),marianas/$(sample)/first-pass.mate-position-sorted.txt)

JAVA = /home/${USER}/share/usr/jdk1.8.0_74/bin/java
MARIANAS = /home/${USER}/share/usr/marianas-1.8.1/Marianas-1.8.1.jar
WALTZ = /home/${USER}/share/usr/marianas-1.8.1/Waltz-3.0.jar
BED_FILE=/home/${USER}/share/reference/target_panels/MSK-ACCESS-v1_0-probe-AB.waltz.bed

define genotype-and-collapse
marianas/$1/$1-pileup.txt : marianas/$1/$1.bam
	$$(call RUN,-c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  cd marianas/$1 && \
									  $(JAVA) -server -Xms2G -Xmx8G -cp $(WALTZ) org.mskcc.juber.waltz.Waltz PileupMetrics 20 $1.bam $(REF_FASTA) $(BED_FILE) && \
									  cd ../..")
									  
marianas/$1/first-pass.mate-position-sorted.txt : marianas/$1/$1-pileup.txt
	$$(call RUN,-c -n 1 -s 8G -m 12G,"set -o pipefail && \
									  cd marianas/$1 && \
									  $(JAVA) -server -Xms2G -Xmx8G -cp $(MARIANAS) org.mskcc.marianas.umi.duplex.DuplexUMIBamToCollapsedFastqFirstPass $1.bam $1-pileup.txt 1 20 0 1 90 $(REF_FASTA) && \
									  sort -S 8G -k \"6,6n\" -k \"8,8n\" first-pass.txt > first-pass.mate-position-sorted.txt && \
									  cd ../..")

endef
$(foreach sample,$(SAMPLES),\
	$(eval $(call genotype-and-collapse,$(sample))))
