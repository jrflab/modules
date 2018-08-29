include modules/Makefile.inc

LOGDIR ?= log/macs2TN.$(NOW)
PHONY += macs2 macs2/broadpeaks macs2/narrowpeaks

macs2TN : $(foreach pair,$(SAMPLE_PAIRS),macs2/$(pair).timestamp)

define macs2-case-control
macs2/broadpeaks/$1_$2.task.complete : bam/$1.bam bam/$2.bam
	$$(call RUN,-c -s 8G -m 12G,"macs2 callpeak -t $$< -c $$(<<) -f BAM -g hs --keep-dup all --broad --outdir macs2/broadpeaks -n $1_$2 -B --verbose 2 --nomodel -p 0.1 && echo $$< $$(<<) > macs2/broadpeaks/$1_$2.task.complete")
	
macs2/narrowpeaks/$1_$2.task.complete : bam/$1.bam bam/$2.bam
	$$(call RUN,-c -s 8G -m 12G,"macs2 callpeak -t $$< -c $$(<<) -f BAM -g hs --keep-dup all --outdir macs2/narrowpeaks -n $1_$2 -B --verbose 2 --nomodel -p 0.1 && echo $$< $$(<<) > macs2/narrowpeaks/$1_$2.task.complete")

macs2/$1_$2.timestamp : macs2/broadpeaks/$1_$2.task.complete macs2/narrowpeaks/$1_$2.task.complete
	$$(call RUN,-c -s 1G -m 1G,"echo $$< $$(<<) > macs2/$1_$2.timestamp")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call macs2-case-control,$(tumor.$(pair)),$(normal.$(pair)))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
