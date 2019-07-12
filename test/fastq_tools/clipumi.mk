include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/clip_umi.$(NOW)
PHONY += marianas

clip_umi : $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)_R1_umi-clipped.fastq.gz) \
		   $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)_R2_umi-clipped.fastq.gz)

JAVA = /home/${USER}/share/usr/jdk1.8.0_74/bin/java
PICARD = /home/${USER}/share/usr/picard/bin/picard.jar
MARIANAS = /home/${USER}/share/usr/marianas/Marianas-1.8.1.jar

define copy-fastq
marianas/$1/$1_R1.fastq.gz marianas/$1/$1_R2.fastq.gz : $3
	$$(call RUN,-n 1 -s 2G -m 4G,"set -o pipefail && \
									  mkdir -p marianas/$1 && \
									  $(RSCRIPT) modules/test/fastq_tools/copyfastq.R --sample_name $1 --fastq_files $3")

endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),$(eval $(call clip-umi-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))
	
	
define clip-umi
marianas/%/%_R1_umi-clipped.fastq.gz marianas/%/%_R2_umi-clipped.fastq.gz : marianas/%/%_R1.fastq.gz marianas/%/%_R2.fastq.gz
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  $(JAVA) -Djava.io.tmpdir=$(TMPDIR) -server -Xms8G -Xmx8G -cp Marianas-1.8.1.jar \
									  org.mskcc.marianas.umi.duplex.fastqprocessing.ProcessLoopUMIFastq \
									  marianas/$$(*)/$$(*)_R1_umi-clipped.fastq.gz marianas/$$(*)/$$(*)_R2_umi-clipped.fastq.gz \
									  3")

endef
 $(foreach sample,$(SAMPLES),\
		$(eval $(call fix-bam,$(sample))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)