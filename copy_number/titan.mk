# titan module
include ~/share/modules/Makefile.inc

LOGDIR = log/titan.$(NOW)

#EXTRACT_ALLELE_READ_COUNTS = $(RSCRIPT) $(HOME)/share/scripts/extractTitanAlleleReadCounts.R
EXTRACT_ALLELE_READ_COUNTS = $(ANACONDA_PYTHON) $(HOME)/share/usr/TITANRunner-0.0.3/scripts/count.py
TITAN = $(RSCRIPT) $(HOME)/share/scripts/runTitan.R
TITAN_SEG = $(PERL) $(HOME)/share/usr/TITANRunner-0.0.3/scripts/createTITANsegmentfiles.pl
NUM_CLUSTERS ?= $(shell seq 1 5)
BQ_THRESHOLD ?= 20
MQ_THRESHOLD ?= 20


TITAN_OPTS ?= 
ifdef TARGETS_FILE
TITAN_OPTS += --targetBed $(TARGETS_FILE)
endif

READ_COUNTER = $(HOME)/share/usr/bin/readCounter
MAP_COUNTER = $(HOME)/share/usr/bin/mapCounter
GC_COUNTER = $(HOME)/share/usr/bin/gcCounter

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : all results seg

all : results seg
results : $(foreach i,$(NUM_CLUSTERS),$(foreach pair,$(SAMPLE_PAIRS),titan/results/$(pair).titan_$i.txt))
seg : $(foreach i,$(NUM_CLUSTERS),$(foreach pair,$(SAMPLE_PAIRS),titan/seg/$(pair).titan_$i.seg))

titan/wig/%.wig : bam/%.bam
	$(call LSCRIPT_MEM,6G,8G,"$(READ_COUNTER) -c $(subst $( ),$(,),$(strip $(CHROMOSOMES))) $< > $@")

define titan-tumor-normal-numcluster
vcf/$1_$2.gatk_het.vcf : vcf/$1_$2.gatk_snps.het_ft.pass.vcf bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_PARALLEL_MEM,4,2.5G,3G,"$$(call GATK_MEM,8G) -T UnifiedGenotyper -nt 4 -R $$(REF_FASTA) --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^), -I $$(bam) ) -L $$< -o $$@ --output_mode EMIT_ALL_SITES")

titan/allele_count/$1_$2.ac.txt : bam/$1.bam vcf/$1_$2.gatk_het.vcf
	$$(call LSCRIPT_MEM,4G,6G,"$$(EXTRACT_ALLELE_READ_COUNTS) $$(<<) $$< $$(REF_FASTA) $$(BQ_THRESHOLD) $$(MQ_THRESHOLD) > $$@")

titan/results/$1_$2.titan_$3.txt : titan/wig/$1.wig titan/wig/$2.wig titan/allele_count/$1_$2.ac.txt
	$$(call LSCRIPT_PARALLEL_MEM,8,1G,1.5G,"$$(TITAN) $$(TITAN_OPTS) --gcWig $$(HMMCOPY_GC_WIG) --mapWig $$(HMMCOPY_MAP_WIG)  --numClusters $3 --tumorWig $$< --normalWig $$(<<) --numCores 8 --outPrefix titan/results/$1_$2 $$(<<<)")

titan/seg/%.titan_$3.seg : titan/results/%.titan_$3.txt
	$$(call LSCRIPT_MEM,4G,6G,"$$(TITAN_SEG) -id=$$*_cluster$3 -infile=$$< -outfile=$$(@:.seg=.txt) -outIGV=$$@")
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(foreach i,$(NUM_CLUSTERS), \
	$(eval $(call titan-tumor-normal-numcluster,$(tumor.$(pair)),$(normal.$(pair)),$i))))


include ~/share/modules/variant_callers/gatk.mk
