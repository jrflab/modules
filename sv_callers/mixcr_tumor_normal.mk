include modules/Makefile.inc

LOGDIR = log/mixcr_tumor_normal.$(NOW)

mixcr : $(foreach sample,$(SAMPLES),mixcr/$(sample)/alignments.vdjca) \
		$(foreach sample,$(SAMPLES),mixcr/$(sample)/alignments_rescued_1.vdjca) \
		$(foreach sample,$(SAMPLES),mixcr/$(sample)/alignments_rescued_2.vdjca) \
		$(foreach sample,$(SAMPLES),mixcr/$(sample)/alignments_rescued_2_extended.vdjca) \
		$(foreach sample,$(SAMPLES),mixcr/$(sample)/clones.clns) \
		$(foreach sample,$(SAMPLES),mixcr/$(sample)/clones.tsv)

define mixcr-tumor-normal
mixcr/$1/alignments.vdjca : bam/$1.bam
	$$(call RUN,-n 8 -s 6G -m 8G -v $(MIXCR_ENV) -w 72:00:00,"set -o pipefail && \
															  mixcr align \
															  -Xmx48g \
															  --species hsa \
															  --preset rna-seq \
															  --dna \
															  -OallowPartialAlignments=true \
															  --threads 8 \
															  --verbose \
															  $$(<) \
															  $$(@)")

mixcr/$1/alignments_rescued_1.vdjca : mixcr/$1/alignments.vdjca
	$$(call RUN,-n 8 -s 4G -m 6G -v $(MIXCR_ENV) -w 24:00:00,"set -o pipefail && \
															  mixcr assemblePartial \
															  -OkValue=12 \
															  -OkOffset=-7 \
															  -OminimalAssembleOverlap=12 \
															  -OminimalNOverlap=5 \
															  $$(<) \
															  $$(@)")
								  
mixcr/$1/alignments_rescued_2.vdjca : mixcr/$1/alignments_rescued_1.vdjca
	$$(call RUN,-n 8 -s 4G -m 6G -v $(MIXCR_ENV) -w 24:00:00,"set -o pipefail && \
															  mixcr assemblePartial \
															  -OkValue=10 \
															  -OkOffset=-5 \
															  -OminimalAssembleOverlap=10 \
															  -OminimalNOverlap=3 \
															  $$(<) \
															  $$(@)")
								  
mixcr/$1/alignments_rescued_2_extended.vdjca : mixcr/$1/alignments_rescued_2.vdjca
	$$(call RUN,-n 8 -s 4G -m 6G -v $(MIXCR_ENV) -w 24:00:00,"set -o pipefail && \
															  mixcr extend \
															  $$(<) \
															  $$(@)")

mixcr/$1/clones.clns : mixcr/$1/alignments_rescued_2_extended.vdjca
	$$(call RUN,-n 8 -s 2G -m 4G -v $(MIXCR_ENV),"set -o pipefail && \
											      mixcr assemble \
											      -OminimalClonalSequenceLength=10 \
											      -ObadQualityThreshold=0 \
											      $$(<) \
											      $$(@)")
								  
mixcr/$1/clones.tsv : mixcr/$1/clones.clns
	$$(call RUN,-n 8 -s 2G -m 4G -v $(MIXCR_ENV),"set -o pipefail && \
											      mixcr exportClones \
											      $$(<) \
											      $$(@)")
								  
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call mixcr-tumor-normal,$(sample))))


..DUMMY := $(shell mkdir -p version; \
	     echo 'mixcr' > version/mixcr_tumor_normal.txt)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: mixcr
