# titan module
include ~/share/modules/Makefile.inc

LOGDIR = log/titan.$(NOW)

#EXTRACT_ALLELE_READ_COUNTS = $(RSCRIPT) $(HOME)/share/scripts/extractTitanAlleleReadCounts.R
EXTRACT_ALLELE_READ_COUNTS = $(ANACONDA_PYTHON) $(HOME)/share/usr/TITANRunner-0.0.3/scripts/count.py
TITAN = $(RSCRIPT) $(HOME)/share/scripts/runTitan.R
TITAN_SEG = $(PERL) $(HOME)/share/usr/TITANRunner-0.0.3/scripts/createTITANsegmentfiles.pl
NUM_CLUSTERS ?= $(shell seq 1 5)
PLOIDY_PRIORS = 2 4 6 8

BQ_THRESHOLD ?= 20
MQ_THRESHOLD ?= 20
TITAN_WINDOW_SIZE ?= 10000

TITAN_SELF_TRANSITION ?= 1e15
TITAN_CLONAL_CLUSTER_TRANSITION ?= 5e5

TITAN_OPTS ?= 
ifdef TARGETS_FILE
TITAN_OPTS += --targetBed $(TARGETS_FILE)
endif

READ_COUNTER = $(HOME)/share/usr/bin/readCounter
MAP_COUNTER = $(HOME)/share/usr/bin/mapCounter
GC_COUNTER = $(HOME)/share/usr/bin/gcCounter

HET_FILTER_SUFFIX := dbsnp.dbsnp_ft
ifdef TARGETS_FILE
HET_FILTER_SUFFIX := $(HET_FILTER_SUFFIX).target_ft
endif

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : titan results seg

titan : results seg
RESULT_FILES := $(foreach i,$(NUM_CLUSTERS),\
	$(foreach j,$(PLOIDY_PRIORS),\
	$(foreach pair,$(SAMPLE_PAIRS),titan/results_w$(TITAN_WINDOW_SIZE)_p$j/$(pair).z$i.titan.txt)))
results : $(RESULT_FILES)
summary : $(foreach j,$(PLOIDY_PRIORS),titan/optclust_results_w$(TITAN_WINDOW_SIZE)_p$j/titan_summary.txt)
seg : $(RESULT_FILES:.txt=.seg)

include ~/share/modules/variant_callers/gatk.mk
include ~/share/modules/variant_callers/samtoolsHet.mk

titan/wig/%.w$(TITAN_WINDOW_SIZE).wig : bam/%.bam
	$(call LSCRIPT_MEM,6G,8G,"$(READ_COUNTER) -w $(TITAN_WINDOW_SIZE) -c $(subst $( ),$(,),$(strip $(CHROMOSOMES))) $< > $@")

titan/wig/gc.w$(TITAN_WINDOW_SIZE).wig :
	$(call LSCRIPT_MEM,6G,8G,"$(GC_COUNTER) -w $(TITAN_WINDOW_SIZE) -c $(subst $( ),$(,),$(strip $(CHROMOSOMES))) $(REF_FASTA) > $@")

titan/wig/map.w$(TITAN_WINDOW_SIZE).wig :
	$(call LSCRIPT_MEM,6G,8G,"$(MAP_COUNTER) -w $(TITAN_WINDOW_SIZE) -c $(subst $( ),$(,),$(strip $(CHROMOSOMES))) $(MAP_BIGWIG) > $@")

define titan-tumor-normal
titan/vcf/$1_$2.gatk_het.vcf : vcf/$2.het_snp.$(HET_FILTER_SUFFIX).pass.vcf bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_PARALLEL_MEM,8,1.5G,3G,"$$(call GATK_MEM2,12G) -T UnifiedGenotyper -nt 8 -R $$(REF_FASTA) --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^), -I $$(bam) ) -L $$< -o $$@ --output_mode EMIT_ALL_SITES")

titan/allele_count/$1_$2.ac.txt : bam/$1.bam titan/vcf/$1_$2.gatk_het.vcf
	$$(call LSCRIPT_MEM,4G,6G,"$$(EXTRACT_ALLELE_READ_COUNTS) $$(<<) $$< $$(REF_FASTA) $$(BQ_THRESHOLD) $$(MQ_THRESHOLD) > $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call titan-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

define titan-tumor-normal-numcluster-ploidy-windowsize
titan/results_w$5_p$4/$1_$2.z$3.titan.txt : titan/wig/$1.w$5.wig titan/wig/$2.w$5.wig titan/allele_count/$1_$2.ac.txt titan/wig/gc.w$5.wig titan/wig/map.w$5.wig
	$$(call LSCRIPT_PARALLEL_MEM,8,1G,1.5G,"$$(TITAN) $$(TITAN_OPTS) --gcWig $$(4<) --mapWig $$(5<) --numClusters $3 --tumorWig $$< --normalWig $$(<<) --ploidyPrior $4 --txnZstrength $$(TITAN_CLONAL_CLUSTER_TRANSITION) --txnExpLen $$(TITAN_SELF_TRANSITION) --numCores 8 --outPrefix titan/results_w$5_p$4/$1_$2.z$3 --plotPrefix titan/results_w$5_p$4/$1_$2.z$3 $$(<<<)")
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(foreach i,$(NUM_CLUSTERS), \
		$(foreach j,$(PLOIDY_PRIORS), \
			$(eval $(call titan-tumor-normal-numcluster-ploidy-windowsize,$(tumor.$(pair)),$(normal.$(pair)),$i,$j,$(TITAN_WINDOW_SIZE))))))

titan/optclust_results_%/titan_summary.txt : $(foreach pair,$(SAMPLE_PAIRS),$(foreach i,$(NUM_CLUSTERS),titan/results_%/$(pair).z$i.titan.txt))
	$(SUMMARIZE_TITAN) --outDir $(@D) $^ 

%.titan.seg %.titan_seg.txt : %.titan.txt
	$(call LSCRIPT_MEM,4G,6G,"$(TITAN_SEG) -id=$(notdir $*) -infile=$< -outfile=$(@:.seg=_seg.txt) -outIGV=$@")

