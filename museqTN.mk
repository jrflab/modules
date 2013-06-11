# vim: set ft=make :
# museq module for use by jsm.mk
SAMPLE_PAIR_FILE = sample_pairs.txt

TUMOR_SAMPLES ?= $(shell cut -f 1 $(SAMPLE_PAIR_FILE))
NORMAL_SAMPLES ?= $(shell cut -f 2 $(SAMPLE_PAIR_FILE))
SAMPLES ?= $(TUMOR_SAMPLES) $(NORMAL_SAMPLES)
NSAMPLES ?= $(words $(TUMOR_SAMPLES))

MUSEQ_THRESHOLD = 0.5
NCHUNKS = 100
MUSEQ_MODEL = /genesis/scratch/sohrab_temp/jknaggs_tmp/models/models/all/model_nov26_1000.pk
EXTRACT_FEATURES = PYTHONPATH=/genesis/scratch/sohrab_temp/jknaggs_tmp/lib/python LD_LIBRARY_PATH=/genesis/scratch/sohrab_temp/jknaggs_tmp/lib /genesis/scratch/sohrab_temp/jknaggs_tmp/alan/mutationSeq2_test_parameters/extract_features
MUSEQ = PYTHONPATH=/genesis/scratch/sohrab_temp/jknaggs_tmp/lib/python/ LD_LIBRARY_PATH=/genesis/scratch/sohrab_temp/jknaggs_tmp/lib /genesis/scratch/sohrab_temp/jknaggs_tmp/alan/mutationSeq2_test_parameters/muSeq


# chunk-ify filtered jsm as input
define museq-chunks
museq/chunks/%.chunk_$1.txt : jsm/tables/%.jsm.filtered.txt
	$$(INIT) sed '1d' $$< | awk "NR % $$(NCHUNKS) == $1 - 1 { print }" > $$@
endef
$(foreach i,$(shell seq 1 $(NCHUNKS)),$(eval $(call museq-chunks,$i)))

define museq-features-chunks
museq/chunk_features/$1_$2.features.chunk_$3.txt : $1.bam $2.bam museq/chunks/$1_$2.chunk_$3.txt
	$$(call INIT_MEM,5G,7G) $$(EXTRACT_FEATURES) $$(word 1,$$^) $$(word 2,$$^) $$(REF_FASTA) --labels normal tumour reference --counts tumour --outfile $$@ --featureset features_classic -p $$(word 3,$$^) &> $$(LOGDIR)/$$(@F).log
endef
$(foreach i,$(shell seq 1 $(NSAMPLES)),$(foreach j,$(shell seq 1 $(NCHUNKS)),$(eval $(call museq-features-chunks,$(word $i,$(NORMAL_SAMPLES)),$(word $i,$(TUMOR_SAMPLES)),$(j)))))

define museq-classify-chunks
museq/chunk_vcf/%.museq.chunk_$1.vcf : museq/chunk_features/%.features.chunk_$1.txt
	$$(call INIT_MEM,5G,7G) $$(MUSEQ) classify $$< --model $$(MUSEQ_MODEL) --out $$@ &> $$(LOGDIR)/$$(@F).log
endef
$(foreach i,$(shell seq 1 $(NCHUNKS)),$(eval $(call museq-classify-chunks,$(i))))

# merge museq chunks
museq/vcf/%.museq.vcf : $(foreach i,$(shell seq 1 $(NCHUNKS)),museq/chunk_vcf/%.museq.chunk_$(i).vcf)
	$(INIT) grep '^#' $< | cut -f -8 > $@.unsorted; cut -f -8 $^ | sed '/^#/d; /^$$/d; s/;$$//' >> $@.unsorted && $(VCF_SORT) $(UCSC_REF_DICT) $@.unsorted > $@ && $(RM) $@.unsorted

museq/vcf/%.museq.filtered.vcf : museq/vcf/%.museq.vcf 
	$(INIT) perl -lane 'if (m/^#/) { print; } else { m/PR=([^;]+)/; $$F[6] = ($$1 > $(MUSEQ_THRESHOLD))? "PASS" : "PR"; print join "\t", @F }' $< > $@

museq/tables/all.%.txt : $(foreach i,$(shell seq 1 $(NSAMPLES)),museq/tables/$(word $i,$(NORMAL_SAMPLES))_$(word $i,$(TUMOR_SAMPLES)).%.txt)
	$(INIT) head -n 1 $< | sed 's/^/TUMOR_SAMPLE\tNORMAL_SAMPLE\t/; s/[^\t]\+\.//g' > $@; \
	for txt in $^; do \
		nsample=`echo $$txt | sed 's/.*\///; s/\..*//; s/_.*//'`; \
		tsample=`echo $$txt | sed 's/.*\///; s/\..*//; s/.*_//'`; \
		sed "1d; s/^/$$tsample\t$$nsample\t/" $$txt >> $@; \
	done

museq/tables/%.txt : museq/vcf/%.vcf
	$(call INIT_MEM,2G,3G) $(VCF_TO_TABLE) $< > $@

include ~/share/modules/vcftools.mk
