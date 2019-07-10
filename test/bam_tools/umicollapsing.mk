include modules/Makefile.inc
include modules/genome_inc/b37.inc

LOGDIR ?= log/umi_collapsing.$(NOW)
PHONY += fgbio

collapsed_umi : $(foreach sample,$(SAMPLES),fgbio/$(sample).qn.sorted.ubam)

JAVA = /home/${USER}/share/usr/jdk1.8.0_74/bin/java
PICARD = /home/${USER}/share/usr/picard/bin/picard.jar
POOL_A_INTERVAL ?= /home/${USER}/share/reference/target_panels/MSK-ACCESS-v1_0-probe-A.sorted.list
POOL_B_INTERVAL ?= /home/${USER}/share/reference/target_panels/MSK-ACCESS-v1_0-probe-B.sorted.list

define fastq-to-ubam
fgbio/$1.qn.sorted.ubam : $3
	$$(call RUN,-c -n 1 -s 8G -m 16G -v $(FGBIO_ENV),"set -o pipefail && \
													  fgbio --tmp-dir $(TMPDIR) -Xms1G -Xmx12G FastqToBam --input $$^ \
													  --read-structures 3M2S+T 3M2S+T \
													  --sample $1 \
													  --output fgbio/$1.ubam \
													  --library $1 && \
													  $(JAVA) -Xmx8G -jar $(PICARD) SortSam \
													  I=fgbio/$1.ubam \
													  O=fgbio/$1.qn.sorted.ubam \
													  SORT_ORDER=queryname \
													  TMP_DIR=$(TMPDIR)")
endef
$(foreach ss,$(SPLIT_SAMPLES),\
	$(if $(fq.$(ss)),$(eval $(call fastq-to-ubam,$(split.$(ss)),$(ss),$(fq.$(ss))))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
