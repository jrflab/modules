# apply oncosnpseq
#
# how to create dbsnp file: 
# java -jar ~/share/usr/snpEff/SnpSift.jar filter 'exists GMAF & GMAF < 0.5' < dbsnp_137.b37.vcf | awk '{OFS="\t"; if (!/^#/){print $1,$2-1,$2,$4"/"$5,"+"}}' > dbsnp_gmafLT50p.bed
# perl -ne 'print if (rand() < .01)' dbsnp_gmafLT50p.bed > dbsnp_gmafLT50p.rand.bed

include ~/share/modules/Makefile.inc

LOGDIR = log/oncoseq.$(NOW)

export LD_LIBRARY_PATH = /home/limr/usr/MATLAB/v82/runtime/glnxa64:/home/limr/usr/MATLAB/v82/bin/glnxa64:/home/limr/usr/MATLAB/v82/sys/os/glnxa64:/home/limr/usr/MATLAB/v82/sys/java/jre/glnxa64/jre/lib/amd64/native_threads:/home/limr/usr/MATLAB/v82/sys/java/jre/glnxa64/jre/lib/amd64/server:/home/limr/usr/MATLAB/v82/sys/java/jre/glnxa64/jre/lib/amd64
export XAPPLRESDIR = /home/limr/usr/MATLAB/v82/X11/app-defaults

MCR_DIR = $(HOME)/share/usr/MATLAB

ONCOSEQ_PATH = $(HOME)/share/usr/oncosnpseq
ONCOSEQ = $(ONCOSEQ_PATH)/executables/oncoseq
PROCESS_PILEUP = $(PERL) $(ONCOSEQ_PATH)/scripts/process_pileup.pl
DBSNP_BED = $(HOME)/share/reference/dbsnp_gmafLT50p.rand.bed
TUMOR_STATES_TABLE = $(ONCOSEQ_PATH)/config/tumourStates.txt
HG_TABLE = $(ONCOSEQ_PATH)/config/hgTables_b37.txt

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all

all: $(foreach pair,$(SAMPLE_PAIRS),oncoseq/$(pair).oncoseq_timestamp)

oncoseq/infile/%.oncoseq.txt.gz : bam/%.bam
	$(call LSCRIPT_MEM,4G,6G,"$(PROCESS_PILEUP) --infile <($(SAMTOOLS) mpileup $<) --outfile >(gzip -c > $@) --snpfile $(DBSNP_BED)")

define oncoseq-tumor-normal
oncoseq/$1_$2.oncoseq_timestamp : oncoseq/infile/$1.oncoseq.txt.gz oncoseq/infile/$2.oncoseq.txt.gz
	$$(call LSCRIPT_MEM,8G,12G,"$$(ONCOSEQ) $$(MCR_DIR) \
		--maxcopy 10 \
		--read_depth_range "[10:40]" \
		--chr_range "[1:22]" \
		--n_train 30000 \
		--maxploidy 4.5 \
		--minploidy 1.5 \
		--normalcontamination \
		--tumourheterogeneity \
		--tumourstatestable $$(TUMOR_STATES_TABLE) \
		--maxnormalcontamination 0.6 \
		--hgtable $$(HG_TABLE) \
		--samplename $1 \
		--infile $$< \
		--outdir $$(@D)/$1_$2 \
		&& touch $@")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call oncoseq-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

