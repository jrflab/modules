include modules/Makefile.inc

LOGDIR ?= log/cnvkit_pca.$(NOW)
PHONY += cnvkit cnvkit/pca

CNVKIT_NORMAL_ON_TARGET ?= $(wildcard $(foreach sample,$(NORMAL_SAMPLES),cnvkit/cnn/normal/$(sample).targetcoverage.cnn))
CNVKIT_NORMAL_OFF_TARGET ?= $(wildcard $(foreach sample,$(NORMAL_SAMPLES),cnvkit/cnn/normal/$(sample).antitargetcoverage.cnn))
CNVKIT_TUMOR_ON_TARGET ?= $(wildcard $(foreach sample,$(TUMOR_SAMPLES),cnvkit/cnn/tumor/$(sample).targetcoverage.cnn))
CNVKIT_TUMOR_OFF_TARGET ?= $(wildcard $(foreach sample,$(TUMOR_SAMPLES),cnvkit/cnn/tumor/$(sample).antitargetcoverage.cnn))

cnvkit : cnvkit/pca/normal_samples_ontarget.pdf cnvkit/pca/normal_samples_offtarget.pdf cnvkit/pca/tumor_samples_ontarget.pdf cnvkit/pca/tumor_samples_offtarget.pdf

cnvkit/pca/normal_samples_ontarget.pdf cnvkit/pca/tumor_samples_ontarget.pdf : $(wildcard cnvkit/cnn/tumor/$(NORMAL_SAMPLES).targetcoverage.cnn) $(wildcard cnvkit/cnn/tumor/$(TUMOR_SAMPLES).targetcoverage.cnn)
	$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 32G -m 48G,"$(RSCRIPT) modules/copy_number/cnvkitprcomp.R --normal_files '$(CNVKIT_NORMAL_ON_TARGET)' --tumor_files '$(CNVKIT_TUMOR_ON_TARGET)' --out_file_normal cnvkit/pca/normal_samples_ontarget.pdf --out_file_tumor cnvkit/pca/tumor_samples_ontarget.pdf")
	
cnvkit/pca/normal_samples_offtarget.pdf cnvkit/pca/tumor_samples_offtarget.pdf : $(wildcard cnvkit/cnn/tumor/$(NORMAL_SAMPLES).antitargetcoverage.cnn) $(wildcard cnvkit/cnn/tumor/$(TUMOR_SAMPLES).antitargetcoverage.cnn)
	$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 32G -m 48G,"$(RSCRIPT) modules/copy_number/cnvkitprcomp.R --normal_files '$(CNVKIT_NORMAL_OFF_TARGET)' --tumor_files '$(CNVKIT_TUMOR_OFF_TARGET)' --out_file_normal cnvkit/pca/normal_samples_offtarget.pdf --out_file_tumor cnvkit/pca/tumor_samples_offtarget.pdf")
				
.PHONY: $(PHONY)
