include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/clip_umi.$(NOW)
PHONY += marianas

clip_umi : $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)_R1.fastq.gz) \
		   $(foreach sample,$(SAMPLES),marianas/$(sample)/$(sample)_R2.fastq.gz)

JAVA = /home/${USER}/share/usr/jdk1.8.0_74/bin/java
PICARD = /home/${USER}/share/usr/picard/bin/picard.jar
MARIANAS = /home/${USER}/share/usr/marianas/Marianas-1.8.1.jar

define clip-umi-fastq
marianas/$1/$1_R1.fastq.gz marianas/$1/$1_R2.fastq.gz : $3
	$$(call RUN,-c -n 1 -s 8G -m 16G,"set -o pipefail && \
									  mkdir -p marianas/$1 && \
									  a=`echo $3 | tr ',' \"\n\"`")
endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),$(eval $(call clip-umi-fastq,$(split.$(ss)),$(ss),$(fq.$(ss))))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
