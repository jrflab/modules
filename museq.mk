# vim: set ft=make :
# apply mutation seq to GATK SNP call results (single sample)

include ~/share/modules/Makefile.inc

SAMPLE_FILE ?= samples.txt
SAMPLES ?= $(shell cat $(SAMPLE_FILE))

MUSEQ_DIR = /genesis/scratch/sohrab_temp/jknaggs_tmp/mutationSeq2
EXTRACT_FEATURES = PYTHONPATH=/genesis/scratch/sohrab_temp/jknaggs_tmp/lib/python LD_LIBRARY_PATH=/genesis/scratch/sohrab_temp/jknaggs_tmp/lib $(MUSEQ_DIR)/extract_features
MUSEQ = PYTHONPATH=/genesis/scratch/sohrab_temp/jknaggs_tmp/lib/python LD_LIBRARY_PATH=/genesis/scratch/sohrab_temp/jknaggs_tmp/lib $(MUSEQ_DIR)/muSeq
MUSEQ_MODEL = /genesis/scratch/sohrab_temp/jknaggs_tmp/all_single/model.pk
MUSEQ_THRESHOLD = 0.5
NCHUNKS = 100
SNP_FILTERS = --filterName "muSeq" --filterExpression "PROB < $(MUSEQ_THRESHOLD)"

VPATH ?= bam

LOGDIR = log/museq.$(NOW)


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: all

#all : $(foreach sample,$(SAMPLES),gatk/vcf/$(sample).snps.annotated.museq.filtered.vcf)
#all : $(foreach sample,$(SAMPLES),gatk/tables/$(sample).snps.annotated.museq.filtered.txt)
all : gatk/tables/all.snps.annotated.museq.filtered.txt

define museq-chunks
museq/chunks/%.chunk_$1.txt : gatk/vcf/%.snps.annotated.vcf
	$$(INIT) sed '/^#/d' $$< | awk "NR % $$(NCHUNKS) == $1 - 1 { print }"  | cut -f1,2 > $$@
endef
$(foreach i,$(shell seq 1 $(NCHUNKS)),$(eval $(call museq-chunks,$i)))

define museq-features-chunks
museq/chunk_features/%.chunk_$1.features : %.bam museq/chunks/%.chunk_$1.txt
	$$(call INIT_MEM,5G,7G) $$(EXTRACT_FEATURES) normal:$$< reference:$$(REF_FASTA) --outfile $$@ --featureset features_single -p $$(word 2,$$^) &> $$(LOG)
endef
$(foreach i,$(shell seq 1 $(NCHUNKS)),$(eval $(call museq-features-chunks,$i)))

define museq-classify-chunks
museq/chunk_vcf/%.museq.chunk_$1.vcf : museq/chunk_features/%.chunk_$1.features
	$$(call INIT_MEM,5G,7G) $$(MUSEQ) classify $$< --model $$(MUSEQ_MODEL) --out $$@ &> $$(LOGDIR)/$$(@F).log
endef
$(foreach i,$(shell seq 1 $(NCHUNKS)),$(eval $(call museq-classify-chunks,$(i))))

# merge museq chunks
museq/vcf/%.museq.vcf : $(foreach i,$(shell seq 1 $(NCHUNKS)),museq/chunk_vcf/%.museq.chunk_$(i).vcf)
	$(INIT) grep '^#' $< | cut -f -8 > $@.unsorted; cut -f -8 $^ | sed '/^#/d; /^$$/d; s/;$$//' >> $@.unsorted && $(VCF_SORT) $(UCSC_REF_DICT) $@.unsorted > $@ && $(RM) $@.unsorted

museq/vcf/%.museq.filtered.vcf : museq/vcf/%.museq.vcf 
	$(INIT) perl -lane 'if (m/^#/) { print; } else { m/PR=([^;]+)/; $$F[6] = ($$1 > $(MUSEQ_THRESHOLD))? "PASS" : "PR"; print join "\t", @F }' $< > $@

gatk/tables/all.%.txt : $(foreach sample,$(SAMPLES),museq/tables/$(sample).%.txt)
	$(call INIT_MEM,2G,3G) $(RBIND) $^ > $@

# add in museq probablility and re-filter
gatk/vcf/%.snps.annotated.museq.vcf : gatk/vcf/%.snps.annotated.vcf museq/vcf/%.museq.vcf
	$(INIT) grep '^##' $< > $@; echo '##INFO=<ID=PROB,Number=1,Type=Float,Description="MutationSeq Probability">' >> $@; grep '^#' $< | sed '/^##/d' >> $@; perl -e 'open IN, $$ARGV[0]; open IN2, $$ARGV[1]; while (<IN2>) { chomp; next if /^#/; @F = split /\t/;  $$M{$$F[0]}{$$F[1]} = $$F[7]; } $$i = 0; while (<IN>) { next if /^#/; @F = split /\t/;  $$F[7] = $$M{$$F[0]}{$$F[1]} . ";" . $$F[7] if exists $$M{$$F[0]}{$$F[1]} ; print join "\t", @F; } close IN; close IN2;' $^ >> $@

gatk/vcf/%.snps.annotated.museq.filtered.vcf : gatk/vcf/%.snps.annotated.museq.vcf gatk/vcf/%.snps.annotated.museq.vcf.idx
	$(call INIT_MEM,5G,6G) $(call GATK_MEM,5G) -T VariantFiltration -R $(REF_FASTA) $(SNP_FILTERS) -o $@ --variant $< &> $(LOG)

gatk/tables/%.txt : gatk/vcf/%.vcf
	$(call INIT_MEM,2G,3G) $(VCF_TO_TABLE) $< > $@

%.vcf.idx : %.vcf
	$(call INIT_MEM,2G,3G) $(IGVTOOLS) index $< &> $(LOGDIR)/$(@F).log

gatk/tables/all.%.txt : $(foreach sample,$(SAMPLES),gatk/tables/$(sample).%.txt)
	$(call INIT_MEM,2G,3G) $(RBIND) $^ > $@

