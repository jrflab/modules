include modules/Makefile.inc

LOGDIR ?= log/cnvkit_qc.$(NOW)
PHONY += cnvkit cnvkit/qc

CNVKIT_CNR ?= $(wildcard $(foreach sample,$(SAMPLES),cnvkit/cnr/$(sample).cnr))

cnvkit : cnvkit/qc/ontarget.tsv cnvkit/qc/offtarget.tsv

cnvkit/qc/ontarget.tsv : $(wildcard cnvkit/cnr/$(SAMPLES).cnr)
	$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 8G -m 16G,"$(RSCRIPT) modules/copy_number/cnvkitqc.R --in_file '$(CNVKIT_CNR)' --out_file cnvkit/qc/ontarget.tsv --option 1")
	
cnvkit/qc/offtarget.tsv : $(wildcard cnvkit/cnr/$(SAMPLES).cnr)
	$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 8G -m 16G,"$(RSCRIPT) modules/copy_number/cnvkitqc.R --in_file '$(CNVKIT_CNR)' --out_file cnvkit/qc/offtarget.tsv --option 0")
				
.PHONY: $(PHONY)
