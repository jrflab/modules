# Run cloneHD on tumour-normal samples
# reconstruct clonal composition
##### DEFAULTS ######
include ~/share/modules/Makefile.inc
include ~/share/modules/variant_callers/gatk.inc

CLONEHD = $(HOME)/share/usr/bin/cloneHD
FILTERHD = $(HOME)/share/usr/bin/filterHD
PREFILTER = $(HOME)/share/usr/bin/pre-filter

TABLE_TO_CLONEHD = $(PERL) $(HOME)/share/scripts/tableToCloneHDFormat.pl

MAX_TOTAL_COPY_NUM ?= 4
MAX_SUBCLONE_NUM ?= 4
NUM_TRIALS ?= 2
NUM_RESTARTS ?= 10

LOGDIR = log/clonehd.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all

all : $(foreach s,$(SAMPLE_SETS),clonehd/results/$s.summary.txt) $(foreach s,$(SAMPLE_SETS),clonehd/results/$s.snv.summary.txt)

include ~/share/modules/variant_callers/gatk.mk
include ~/share/modules/variant_callers/somatic/mutect.mk

clonehd/cov/%.cov.txt : bam/%.bam
	$(call LSCRIPT_MEM,4G,7G,"$(SAMTOOLS) bedcov $(TARGETS_FILE) $< > $@")

clonehd/cna/%.cna.txt : clonehd/cov/%.cov.txt
	$(call LSCRIPT_MEM,2G,3G,"awk 'BEGIN { OFS = \"\t\" } {print \$$1$(,)\$$3$(,)int(0.5+\$$4/1000.0)$(,)1}' $< > $@")

clonehd/cna/%.cna.pref.txt : clonehd/cna/%.cna.txt
	$(call LSCRIPT_MEM,2G,4G,"$(PREFILTER) --data $< --pre clonehd/cna/$*.cna --print-tracks 1")

clonehd/cna/%.cna.posterior-1.txt : clonehd/cna/%.cna.txt
	$(call LSCRIPT_MEM,4G,6G,"$(FILTERHD) --data $< --mode 3 --pre clonehd/cna/$*.cna --filter-shortSeg 1")


# we need depth for positions in tumours that are heterozygous in the normal
define clonehd-set-tumors-normal
vcf/$1.gatk_het.vcf : vcf/$3.gatk_snps.het_ft.pass.vcf $$(foreach tumor,$2,bam/$$(tumor).bam) bam/$3.bam
	$$(call LSCRIPT_PARALLEL_MEM,4,2.5G,3G,"$$(call GATK_MEM,8G) -T UnifiedGenotyper -nt 4 -R $$(REF_FASTA) --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^), -I $$(bam) ) -L $$< -o $$@ --output_mode EMIT_ALL_SITES")

clonehd/baf/$1.baf.txt : tables/$1.gatk_het.tab.txt
	$$(INIT) $$(TABLE_TO_CLONEHD) $2 < $$< > $$@

clonehd/cna/$1.cna.txt : $$(foreach tumor,$2,clonehd/cna/$$(tumor).cna.txt)
	$$(INIT) cut -f1,2 $$< > $$@; \
		for x in $$^; do \
			paste $$@ <(cut -f3,4 $$$$x) > $$@.tmp && mv $$@.tmp $$@; \
		done

clonehd/cna/$1.cna.bias.jumps.txt : clonehd/cna/$1.cna.txt clonehd/cna/$3.cna.posterior-1.txt
	$$(call LSCRIPT_MEM,4G,7G,"$$(FILTERHD) --data $$< --mode 3 --pre clonehd/cna/$1.cna.bias --bias $$(<<) --sigma 0 --jump 1")

clonehd/baf/$1.baf.jumps.txt : clonehd/baf/$1.baf.txt
	$$(call LSCRIPT_MEM,8G,17G,"$$(FILTERHD) --data $$< --mode 1 --pre clonehd/baf/$1.baf --sigma 0 --jump 1 --reflect 1 --dist 1")

clonehd/results/$1.summary.txt : clonehd/baf/$1.baf.txt clonehd/cna/$1.cna.txt clonehd/baf/$1.baf.jumps.txt clonehd/cna/$1.cna.bias.jumps.txt clonehd/cna/$3.cna.posterior-1.txt
	$$(call LSCRIPT_MEM,4G,8G,"$$(CLONEHD) --cna $$(<<) --baf $$< \
		--pre clonehd/results/$1 --bias $$(5<) --seed 123 --trials $$(NUM_TRIALS) \
		--nmax $$(MAX_SUBCLONE_NUM) --force --max-tcn $$(MAX_TOTAL_COPY_NUM) --cna-jumps $$(4<) --baf-jumps $$(<<<) \
		--min-jump 0.01 --restarts $$(NUM_RESTARTS) --mass-gauging 1")

clonehd/results/$1.snv.summary.txt : clonehd/snv/$1.snv.txt clonehd/results/$1.summary.txt
	$$(call LSCRIPT_MEM,4G,8G,"$$(CLONEHD) --snv $$< --pre clonehd/results/$1.snv \
		--seed 123 --trials 2 --nmax $$(MAX_SUBCLONE_NUM) --force --max-tcn $$(MAX_TOTAL_COPY_NUM) \
		--mean-tcn clonehd/results/$1.mean-tcn.txt \
		--avail-cn clonehd/results/$1.avail-cn.txt")

vcf/$1.gatk_som.vcf : $$(foreach tumor,$2,bam/$$(tumor).bam vcf/$$(tumor)_$3.mutect.$$(MUTECT_FILTER_SUFFIX).vcf) bam/$3.bam
	$$(call LSCRIPT_PARALLEL_MEM,4,2.5G,3G,"$$(call GATK_MEM,8G) -T UnifiedGenotyper -nt 4 -R $$(REF_FASTA) --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^), -I $$(bam) ) $$(foreach vcf,$$(filter %.vcf,$$^), -L $$(vcf) ) -o $$@ --output_mode EMIT_ALL_SITES")

clonehd/snv/$1.snv.txt : tables/$1.gatk_som.tab.txt
	$$(INIT) $$(TABLE_TO_CLONEHD) $2 < $$< > $$@
endef
$(foreach s,$(SAMPLE_SETS),$(eval $(call clonehd-set-tumors-normal,$s,$(tumor.$s),$(normal.$s))))



