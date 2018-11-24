include modules/Makefile.inc

LOGDIR ?= log/cnvkit_prcomp.$(NOW)
PHONY += cnvkit cnvkit/component

CNVKIT_NORMAL_ON_TARGET ?= $(wildcard $(foreach sample,$(NORMAL_SAMPLES),cnvkit/cnn/normal/$(sample).targetcoverage.cnn))
CNVKIT_NORMAL_OFF_TARGET ?= $(wildcard $(foreach sample,$(NORMAL_SAMPLES),cnvkit/cnn/normal/$(sample).antitargetcoverage.cnn))
CNVKIT_TUMOR_ON_TARGET ?= $(wildcard $(foreach sample,$(TUMOR_SAMPLES),cnvkit/cnn/tumor/$(sample).targetcoverage.cnn))
CNVKIT_TUMOR_OFF_TARGET ?= $(wildcard $(foreach sample,$(TUMOR_SAMPLES),cnvkit/cnn/tumor/$(sample).antitargetcoverage.cnn))

cnvkit : cnvkit/component/normal_samples_ontarget.pdf cnvkit/component/normal_samples_offtarget.pdf cnvkit/component/tumor_samples_ontarget.pdf cnvkit/component/tumor_samples_offtarget.pdf cnvkit/component/all_samples_ontarget.pdf cnvkit/component/all_samples_offtarget.pdf

cnvkit/component/normal_samples_ontarget.pdf : $(wildcard cnvkit/cnn/normal/$(NORMAL_SAMPLES).targetcoverage.cnn)
	$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 32G -m 48G,"$(RSCRIPT) modules/copy_number/cnvkitprcomp.R --in_file '$(CNVKIT_NORMAL_ON_TARGET)' --out_file cnvkit/component/normal_samples_ontarget.pdf")
	
cnvkit/component/normal_samples_offtarget.pdf : $(wildcard cnvkit/cnn/normal/$(NORMAL_SAMPLES).antitargetcoverage.cnn)
	$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 32G -m 48G,"$(RSCRIPT) modules/copy_number/cnvkitprcomp.R --in_file '$(CNVKIT_NORMAL_OFF_TARGET)' --out_file cnvkit/component/normal_samples_offtarget.pdf")
	
cnvkit/component/tumor_samples_ontarget.pdf : $(wildcard cnvkit/cnn/tumor/$(TUMOR_SAMPLES).targetcoverage.cnn)
	$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 32G -m 48G,"$(RSCRIPT) modules/copy_number/cnvkitprcomp.R --in_file '$(CNVKIT_TUMOR_ON_TARGET)' --out_file cnvkit/component/tumor_samples_ontarget.pdf")
	
cnvkit/component/tumor_samples_offtarget.pdf : $(wildcard cnvkit/cnn/tumor/$(TUMOR_SAMPLES).antitargetcoverage.cnn)
	$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 32G -m 48G,"$(RSCRIPT) modules/copy_number/cnvkitprcomp.R --in_file '$(CNVKIT_TUMOR_OFF_TARGET)' --out_file cnvkit/component/tumor_samples_offtarget.pdf")

cnvkit/component/all_samples_ontarget.pdf : $(wildcard cnvkit/cnn/tumor/$(NORMAL_SAMPLES).targetcoverage.cnn) $(wildcard cnvkit/cnn/tumor/$(TUMOR_SAMPLES).targetcoverage.cnn)
	$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 32G -m 48G,"$(RSCRIPT) modules/copy_number/cnvkitcompca.R --normal_files '$(CNVKIT_NORMAL_ON_TARGET)' --tumor_files '$(CNVKIT_TUMOR_ON_TARGET)' --out_file cnvkit/component/all_samples_ontarget.pdf")
	
cnvkit/component/all_samples_offtarget.pdf : $(wildcard cnvkit/cnn/tumor/$(NORMAL_SAMPLES).antitargetcoverage.cnn) $(wildcard cnvkit/cnn/tumor/$(TUMOR_SAMPLES).antitargetcoverage.cnn)
	$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 32G -m 48G,"$(RSCRIPT) modules/copy_number/cnvkitcompca.R --normal_files '$(CNVKIT_NORMAL_OFF_TARGET)' --tumor_files '$(CNVKIT_TUMOR_OFF_TARGET)' --out_file cnvkit/component/all_samples_offtarget.pdf")

				
.PHONY: $(PHONY)
