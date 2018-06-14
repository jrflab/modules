include modules/Makefile.inc

LOGDIR ?= log/macs2broadTN.$(NOW)
PHONY += macs2 macs2/broad

macs2broadTN : $(foreach pair,$(SAMPLE_PAIRS),macs2/broad/$(pair).timestamp)

define macs2-case-control
macs2/broad/$1_$2.timestamp : bam/$1.bam bam/$2.bam
	$$(call RUN,-c -s 8G -m 12G,"macs2 callpeak -t $$< -c $$(<<) -f BAM -g hs --keep-dup all --broad --outdir macs2/broad -n $1_$2 -B --verbose 2 --nomodel -p 0.1 && echo $$< $$(<<) > macs2/broad/$1_$2.timestamp")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call macs2-case-control,$(tumor.$(pair)),$(normal.$(pair)))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
