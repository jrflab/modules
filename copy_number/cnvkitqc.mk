include modules/Makefile.inc

LOGDIR ?= log/cnvkit_qc.$(NOW)
PHONY += cnvkit cnvkit/qc

CNVKIT_NORMAL ?= $(wildcard $(foreach sample,$(NORMAL_SAMPLES),cnvkit/cnr/$(sample).cnr))
CNVKIT_TUMOR ?= $(wildcard $(foreach sample,$(TUMOR_SAMPLES),cnvkit/cnr/$(sample).cnr))
CNVKIT_NORMAL_ON_TARGET ?= $(wildcard $(foreach sample,$(NORMAL_SAMPLES),cnvkit/cnn/normal/$(sample).targetcoverage.cnn))
CNVKIT_NORMAL_OFF_TARGET ?= $(wildcard $(foreach sample,$(NORMAL_SAMPLES),cnvkit/cnn/normal/$(sample).antitargetcoverage.cnn))
CNVKIT_TUMOR_ON_TARGET ?= $(wildcard $(foreach sample,$(TUMOR_SAMPLES),cnvkit/cnn/tumor/$(sample).targetcoverage.cnn))
CNVKIT_TUMOR_OFF_TARGET ?= $(wildcard $(foreach sample,$(TUMOR_SAMPLES),cnvkit/cnn/tumor/$(sample).antitargetcoverage.cnn))

cnvkit : cnvkit/qc/qc_ontarget.tsv cnvkit/qc/qc_offtarget.tsv cnvkit/qc/bin_qc_ontarget.tsv cnvkit/qc/bin_qc_offtarget.tsv

cnvkit/qc/qc_ontarget.tsv : $(wildcard cnvkit/cnr/$(SAMPLES).cnr)
	$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 8G -m 16G,"$(RSCRIPT) modules/copy_number/cnvkitqc.R --normal_files '$(CNVKIT_NORMAL)' --tumor_files '$(CNVKIT_TUMOR)' --out_file cnvkit/qc/qc_ontarget.tsv --option 1")
	
cnvkit/qc/qc_offtarget.tsv : $(wildcard cnvkit/cnr/$(SAMPLES).cnr)
	$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 8G -m 16G,"$(RSCRIPT) modules/copy_number/cnvkitqc.R --normal_files '$(CNVKIT_NORMAL)' --tumor_files '$(CNVKIT_TUMOR)' --out_file cnvkit/qc/qc_offtarget.tsv --option 0")
	
cnvkit/qc/bin_qc_ontarget.tsv : $(wildcard cnvkit/cnn/normal/$(NORMAL_SAMPLES).targetcoverage.cnn) $(wildcard cnvkit/cnn/tumor/$(TUMOR_SAMPLES).targetcoverage.cnn)
	$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 32G -m 48G,"$(RSCRIPT) modules/copy_number/cnvkitbinqc.R --normal_files '$(CNVKIT_NORMAL_ON_TARGET)' --tumor_files '$(CNVKIT_TUMOR_ON_TARGET)' --out_file cnvkit/qc/bin_qc_ontarget.tsv")
	
cnvkit/qc/bin_qc_offtarget.tsv : $(wildcard cnvkit/cnn/normal/$(NORMAL_SAMPLES).antitargetcoverage.cnn) $(wildcard cnvkit/cnn/tumor/$(TUMOR_SAMPLES).antitargetcoverage.cnn)
	$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 32G -m 48G,"$(RSCRIPT) modules/copy_number/cnvkitbinqc.R --normal_files '$(CNVKIT_NORMAL_OFF_TARGET)' --tumor_files '$(CNVKIT_TUMOR_OFF_TARGET)' --out_file cnvkit/qc/bin_qc_offtarget.tsv")


.PHONY: $(PHONY)
