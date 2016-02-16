include modules/Makefile.inc
include modules/variant_callers/somatic/somaticVariantCaller.inc

LOGDIR = log/recurrent_mutations.$(NOW)

EXCEL_NONSYNONYMOUS_MUTECT ?= summary/tsv/mutect_nonsynonymous.tsv
EXCEL_NONSYNONYMOUS_STRELKA_VARSCAN ?= summary/tsv/strelka_varscan_nonsynonymous.tsv

# sufam plot parameters
SUFAM_PLOT_MIN_NR_SAMPLES_WITH_MUTATION?=2
SUFAM_PLOT_SAMPLE_ORDER?=$(SAMPLES)


recurrent_mutations: recurrent_mutations/recurrent_mutations.tsv recurrent_mutations/sufam/all_mutations.vcf recurrent_mutations/sufam/all_sufam.txt recurrent_mutations/sufam/sufam.ipynb recurrent_mutations/sufam/sufam.html

recurrent_mutations/recurrent_mutations.tsv: $(EXCEL_NONSYNONYMOUS_MUTECT) $(EXCEL_NONSYNONYMOUS_STRELKA_VARSCAN)
	$(INIT) unset PYTHONPATH && \
	source $(ANACONDA_27_ENV)/bin/activate $(ANACONDA_27_ENV) && \
	python modules/scripts/recurrent_mutations_plot.py $^ $(@D)

recurrent_mutations/sufam/all_mutations.vcf: $(EXCEL_NONSYNONYMOUS_MUTECT) $(EXCEL_NONSYNONYMOUS_STRELKA_VARSCAN)
	$(INIT) unset PYTHONPATH && \
	source $(ANACONDA_27_ENV)/bin/activate $(ANACONDA_27_ENV) && \
	(csvcut -tc CHROM,POS,ID,REF,ALT,ANN[*].GENE,ANN[*].HGVS_P,ANN[*].IMPACT,ANN[*].EFFECT,GMAF $< | head -1 | sed 's/CHROM/#CHROM/'; \
	 csvcut -tc CHROM,POS,ID,REF,ALT,ANN[*].GENE,ANN[*].HGVS_P,ANN[*].IMPACT,ANN[*].EFFECT,GMAF $< | tail -n +2 | sort -u; \
	 csvcut -tc CHROM,POS,ID,REF,ALT,ANN[*].GENE,ANN[*].HGVS_P,ANN[*].IMPACT,ANN[*].EFFECT,GMAF $(<<) | tail -n +2 | sort -u) | csvformat -T > $@

recurrent_mutations/sufam/%_validated_sufam.txt: recurrent_mutations/sufam/all_mutations.vcf bam/%.bam
	$(INIT) unset PYTHONPATH && \
		source $(SUFAM_ENV)/bin/activate \
			   $(SUFAM_ENV) && \
		sufam --sample_name $* $(REF_FASTA) $< $(<<) > $@

recurrent_mutations/sufam/all_sufam.txt: $(foreach sample,$(SAMPLES),recurrent_mutations/sufam/$(sample)_validated_sufam.txt)
	$(INIT) (head -1 $<; tail -qn +2 $^) > $@

recurrent_mutations/sufam/sufam.ipynb: recurrent_mutations/sufam/all_sufam.txt recurrent_mutations/sufam/all_mutations.vcf
	$(INIT) unset PYTHONPATH && \
		source $(ANACONDA_27_ENV)/bin/activate $(ANACONDA_27_ENV) && \
		cd $(@D) && cat ../../modules/scripts/recurrent_mutations_sufam.ipynb | \
		sed "s:ALL_SUFAM:`basename $<`:" | sed "s:SUFAM_ANNOTATIONS_VCF:`basename $(<<)`:" | sed 's:MIN_NR_SAMPLES_WITH_MUTATION:$(SUFAM_PLOT_MIN_NR_SAMPLES_WITH_MUTATION):' | \
		sed 's:SAMPLE_ORDER:$(SUFAM_PLOT_SAMPLE_ORDER):' | \
		runipy --no-chdir - `basename $@`

recurrent_mutations/sufam/sufam.html: recurrent_mutations/sufam/sufam.ipynb
	$(INIT) unset PYTHONPATH && \
		source $(ANACONDA_27_ENV)/bin/activate $(ANACONDA_27_ENV) && \
		ipython nbconvert $< --to html --stdout > $@

.PHONY: recurrent_mutations
