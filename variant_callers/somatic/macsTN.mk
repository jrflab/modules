include modules/Makefile.inc

LOGDIR ?= log/macsTN.$(NOW)

macsTN : $(foreach pair,$(SAMPLE_PAIRS),macs2/$(pair).timestamp)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: macsTN

VPATH = macs2

define macs2-case-control
macs2/$1_$2.timestamp : bam/$1.bam bam/$2.bam
	$$(call RUN,-c -s 12G -m 16G,"macs2 callpeak -t bam/$1.bam -c bam/$2.bam -f BAM -g hs --keep-dup all --outdir macs2 -n $1_$2 -B --verbose 2 --nomodel -p 0.1 && touch macs2/$1_$2.timestamp")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call macs2-case-control,$(tumor.$(pair)),$(normal.$(pair)))))
