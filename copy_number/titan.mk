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


define hetsnp-chr
chr_vcf/%.$1.het_snp.vcf : bam/%.bam
	$$(call LSCRIPT_MEM,6G,8G,"$$(SAMTOOLS2) mpileup -r $1 -f $$(REF_FASTA) -g -I $$< | $$(BCFTOOLS2) call -c | $$(BCFTOOLS2) view -g het | $$(VCFUTILS) varFilter -d 10 -a 5 - > $$@")
endef
$(foreach chr,$(CHROMOSOMES),$(eval $(call hetsnp-chr,$(chr))))

vcf/%.het_snp.vcf : $(foreach chr,$(CHROMOSOMES),chr_vcf/%.$(chr).het_snp.vcf)
	$(INIT) grep -P '^#' $< > $@ && sed '/^#/d' $^ >> $@

define titan-tumor-normal
vcf/$1_$2.gatk_het.vcf : vcf/$2.het_snp.vcf bam/$1.bam bam/$2.bam
	$$(call LSCRIPT_PARALLEL_MEM,4,2.5G,3G,"$$(call GATK_MEM,8G) -T UnifiedGenotyper -nt 4 -R $$(REF_FASTA) --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^), -I $$(bam) ) -L $$< -o $$@ --output_mode EMIT_ALL_SITES")

titan/allele_count/$1_$2.ac.txt : bam/$1.bam vcf/$1_$2.gatk_het.vcf
	$$(call LSCRIPT_MEM,4G,6G,"$$(EXTRACT_ALLELE_READ_COUNTS) $$(<<) $$< $$(REF_FASTA) $$(BQ_THRESHOLD) $$(MQ_THRESHOLD) > $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(eval $(call titan-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

define titan-tumor-normal-numcluster
titan/results/$1_$2.titan_$3.txt : titan/wig/$1.wig titan/wig/$2.wig titan/allele_count/$1_$2.ac.txt
	$$(call LSCRIPT_PARALLEL_MEM,8,1G,1.5G,"$$(TITAN) $$(TITAN_OPTS) --gcWig $$(HMMCOPY_GC_WIG) --mapWig $$(HMMCOPY_MAP_WIG)  --numClusters $3 --tumorWig $$< --normalWig $$(<<) --numCores 8 --outPrefix titan/results/$1_$2 --plotPrefix titan/results/$1_$2 $$(<<<)")
endef
$(foreach pair,$(SAMPLE_PAIRS), \
	$(foreach i,$(NUM_CLUSTERS), \
	$(eval $(call titan-tumor-normal-numcluster,$(tumor.$(pair)),$(normal.$(pair)),$i))))

define titan-numcluster
titan/seg/%.titan_$1.seg : titan/results/%.titan_$1.txt
	$$(call LSCRIPT_MEM,4G,6G,"$$(TITAN_SEG) -id=$$*_cluster$3 -infile=$$< -outfile=$$(@:.seg=.txt) -outIGV=$$@")
endef
$(foreach i,$(NUM_CLUSTERS),$(eval $(call titan-numcluster,$i)))

include ~/share/modules/variant_callers/gatk.mk
