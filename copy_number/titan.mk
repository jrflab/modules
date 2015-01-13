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
TITAN_WINDOW_SIZE ?= 10000

TITAN_SELF_TRANSITION ?= 1e15
TITAN_CLONAL_CLUSTER_TRANSITION ?= 5e3

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

include ~/share/modules/variant_callers/gatk.mk
include ~/share/modules/variant_callers/samtoolsHet.mk

titan/wig/%.w$(TITAN_WINDOW_SIZE).wig : bam/%.bam
	$(call LSCRIPT_MEM,6G,8G,"$(READ_COUNTER) -w $(TITAN_WINDOW_SIZE) -c $(subst $( ),$(,),$(strip $(CHROMOSOMES))) $< > $@")

titan/wig/gc.w$(TITAN_WINDOW_SIZE).wig :
	$(call LSCRIPT_MEM,6G,8G,"$(GC_COUNTER) -w $(TITAN_WINDOW_SIZE) -c $(subst $( ),$(,),$(strip $(CHROMOSOMES))) $(REF_FASTA) > $@")

titan/wig/map.w$(TITAN_WINDOW_SIZE).wig :
	$(call LSCRIPT_MEM,6G,8G,"$(MAP_COUNTER) -w $(TITAN_WINDOW_SIZE) -c $(subst $( ),$(,),$(strip $(CHROMOSOMES))) $(MAP_BIGWIG) > $@")

define titan-tumor-normal
titan/vcf/$1_$2.gatk_het.vcf : vcf/$2.het_snp.vcf bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_PARALLEL_MEM,4,2.5G,3G,"$$(call GATK_MEM2,8G) -T UnifiedGenotyper -nt 4 -R $$(REF_FASTA) --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^), -I $$(bam) ) -L $$< -o $$@ --output_mode EMIT_ALL_SITES")

titan/allele_count/$1_$2.ac.txt : bam/$1.bam titan/vcf/$1_$2.gatk_het.vcf
	$$(call LSCRIPT_MEM,4G,6G,"$$(EXTRACT_ALLELE_READ_COUNTS) $$(<<) $$< $$(REF_FASTA) $$(BQ_THRESHOLD) $$(MQ_THRESHOLD) > $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call titan-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

define titan-tumor-normal-numcluster
titan/results/$1_$2.titan_$3.txt : titan/wig/$1.w$4.wig titan/wig/$2.w$4.wig titan/allele_count/$1_$2.ac.txt titan/wig/gc.w$4.wig titan/wig/map.w$4.wig
	$$(call LSCRIPT_PARALLEL_MEM,8,1G,1.5G,"$$(TITAN) $$(TITAN_OPTS) --gcWig $$(4<) --mapWig $$(5<) --numClusters $3 --tumorWig $$< --normalWig $$(<<) --txnZstrength $$(TITAN_CLONAL_CLUSTER_TRANSITION) --txnExpLen $$(TITAN_SELF_TRANSITION) --numCores 8 --outPrefix titan/results/$1_$2 --plotPrefix titan/results/$1_$2 $$(<<<)")
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(foreach i,$(NUM_CLUSTERS), \
	$(eval $(call titan-tumor-normal-numcluster,$(tumor.$(pair)),$(normal.$(pair)),$i,$(TITAN_WINDOW_SIZE)))))

define titan-numcluster
titan/seg/%.titan_$1.seg : titan/results/%.titan_$1.txt
	$$(call LSCRIPT_MEM,4G,6G,"$$(TITAN_SEG) -id=$$*_cluster$3 -infile=$$< -outfile=$$(@:.seg=.txt) -outIGV=$$@")
endef
$(foreach i,$(NUM_CLUSTERS),$(eval $(call titan-numcluster,$i)))

