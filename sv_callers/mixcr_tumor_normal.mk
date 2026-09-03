include modules/Makefile.inc

LOGDIR = log/mixcr_tumor_normal.$(NOW)

mixcr : $(foreach sample,$(SAMPLES),mixcr/$(sample)/$(sample).1.fastq.gz) \
		$(foreach sample,$(SAMPLES),mixcr/$(sample)/alignments.vdjca) \
		$(foreach sample,$(SAMPLES),mixcr/$(sample)/alignments_rescued_1.vdjca) \
		$(foreach sample,$(SAMPLES),mixcr/$(sample)/alignments_rescued_2.vdjca) \
		$(foreach sample,$(SAMPLES),mixcr/$(sample)/alignments_rescued_2_extended.vdjca) \
		$(foreach sample,$(SAMPLES),mixcr/$(sample)/clones.clns) \
		$(foreach sample,$(SAMPLES),mixcr/$(sample)/clones.tsv)

VDJ_LOCI_BED ?= $(HOME)/share/lib/bed_files/hg19_vdj_loci.bed

define extract-fastq
mixcr/$1/$1.txt : bam/$1.bam
	$$(call RUN,-n 1 -s 4G -m 6G,"set -o pipefail && \
								  mkdir -p mixcr/$1 && \
								  { $$(SAMTOOLS) view -L $(VDJ_LOCI_BED) $$(<) | cut -f1; \
								    $$(SAMTOOLS) view -f 4 $$(<) | cut -f1; } | sort -u > $$(@)")

mixcr/$1/$1.bam : bam/$1.bam mixcr/$1/$1.txt
	$$(call RUN,-n 4 -s 4G -m 9G,"set -o pipefail && \
								  mkdir -p mixcr/$1 && \
								  $$(SAMTOOLS) view -b -N $$(<<) $$(<) | \
								  $$(SAMTOOLS) sort -T mixcr/$1/$1 -O bam -n -@ 4 -m 6G -o $$(@) -")

mixcr/$1/$1.1.fastq : mixcr/$1/$1.bam
	$$(call RUN,-n 4 -s 4G -m 9G,"set -o pipefail && \
								  bedtools bamtofastq -i $$(<) -fq mixcr/$1/$1.1.fastq -fq2 mixcr/$1/$1.2.fastq")

mixcr/$1/$1.1.fastq.gz : mixcr/$1/$1.1.fastq
	$$(call RUN,-n 4 -s 4G -m 9G,"set -o pipefail && \
								      gzip mixcr/$1/$1.1.fastq && \
								      gzip mixcr/$1/$1.2.fastq")
					      
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call extract-fastq,$(sample))))

define mixcr-tumor-normal
mixcr/$1/alignments.vdjca : mixcr/$1/$1.1.fastq.gz
	$$(call RUN,-n 8 -s 4G -m 6G -v $(MIXCR_ENV) -w 24:00:00,"set -o pipefail && \
															  mixcr align \
															  --species hsa \
															  --preset rna-seq \
															  --dna \
															  -OallowPartialAlignments=true \
															  --threads 8 \
															  --verbose \
															  mixcr/$1/$1.1.fastq.gz mixcr/$1/$1.2.fastq.gz \
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
