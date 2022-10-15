include modules/Makefile.inc

LOGDIR = log/kallisto.$(NOW)

kallisto : $(foreach sample,$(SAMPLES),kallisto/$(sample)/$(sample)_R1.fastq.gz) \
	   $(foreach sample,$(SAMPLES),kallisto/$(sample)/$(sample)_R2.fastq.gz) \
	   $(foreach sample,$(SAMPLES),kallisto/$(sample)/abundance.tsv) \
	   kallisto/tpm_bygene.txt

SLEUTH_ANNOT ?= $(HOME)/share/lib/resource_files/Hugo_ENST_ensembl75_fixed.txt
KALLISTO_INDEX ?= $(HOME)/share/lib/ref_files/b37/ensembl_v75-0.43.0_kallisto_index

define merge-fastq
kallisto/$1/$1_R1.fastq.gz : $$(foreach split,$2,$$(word 1, $$(fq.$$(split))))
	$$(call RUN,-c -n 12 -s 0.5G -m 1G -w 24:00:00 -v $(PIGZ_ENV),"set -o pipefail && \
								       $$(PIGZ) -cd -p 12 $$(^) | $$(PIGZ) -c -p 12 > $$(@)")
	
kallisto/$1/$1_R2.fastq.gz : $$(foreach split,$2,$$(word 2, $$(fq.$$(split))))
	$$(call RUN,-c -n 12 -s 0.5G -m 1G -w 24:00:00 -v $(PIGZ_ENV),"set -o pipefail && \
								       $$(PIGZ) -cd -p 12 $$(^) | $$(PIGZ) -c -p 12 > $$(@)")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call merge-fastq,$(sample),$(split.$(sample)))))

define fastq-to-kallisto
kallisto/$1/abundance.tsv : kallisto/$1/$1_R1.fastq.gz kallisto/$1/$1_R2.fastq.gz
	$$(call RUN,-c -n 12 -s 2G -m 3G -v $(KALLISTO_ENV),"set -o pipefail && \
							     kallisto quant \
							     -i $$(KALLISTO_INDEX) \
							     -o kallisto/$1 \
							     --bias -b 100 -t 12\
							     --fusion $$(<) $$(<<)")

endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call fastq-to-kallisto,$(sample))))
		
kallisto/tpm_bygene.txt : $(foreach sample,$(SAMPLES),kallisto/$(sample)/abundance.tsv)
	$(call RUN, -c -n 24 -s 1G -m 2G -v $(KALLISTO_ENV),"set -o pipefail && \
							     $(RSCRIPT) $(SCRIPTS_DIR)/summarize_sleuth.R --annotation $(SLEUTH_ANNOT) --samples '$(SAMPLES)'")

..DUMMY := $(shell mkdir -p version; \
	     $(SAMTOOLS) --version > version/kallisto.txt; \
	     ~/share/usr/env/kallisto-0.46.2/bin/kallisto version >> version/kallisto.txt; \
	     ~/share/usr/env/kallisto-0.46.2/bin/R --version >> version/kallisto.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: kallisto
