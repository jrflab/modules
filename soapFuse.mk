# soapfuse
# vim: set ft=make :

include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

LOGDIR = log/soapfuse.$(NOW)

SOAPFUSE = $(HOME)/usr/SOAPfuse-v1.26/SOAPfuse-RUN.pl 
SOAPFUSE_CONFIG = $(HOME)/share/usr/SOAPfuse-v1.26/config/config.txt
PREPARE_SOAPFUSE = $(HOME)/share/scripts/prepareSoapFuse.pl

SOAPFUSE_NORMAL_FILTER = $(PERL) $(HOME)/share/scripts/normalFilterSoapFuse.pl
SOAPFUSE_NORMAL_FILTER_OPTS = -w 1000

ONCOFUSE_MEM = $(JAVA) -Xmx$1 -jar $(HOME)/share/usr/oncofuse-v1.0.3/Oncofuse.jar
ONCOFUSE_TISSUE_TYPE ?= EPI

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all sample_tables soapfuse

ifdef NORMAL_SOAPFUSE_RESULTS
ALL_SUFFIX = nft.oncofuse.merged
else
ALL_SUFFIX = oncofuse.merged
endif
all : soapfuse/alltables/all.sfuse.$(ALL_SUFFIX).txt soapfuse/alltables/all.isoform_sfuse.$(ALL_SUFFIX).txt
sample_tables : $(foreach sample,$(SAMPLES),soapfuse/tables/$(sample).sfuse.txt soapfuse/tables/$(sample).isoform_sfuse.txt)
soapfuse : $(foreach sample,$(SAMPLES),soapfuse/$(sample).timestamp)

soapfuse/%.timestamp : soapfuse/sample_lists/%.txt
	$(call LSCRIPT_NAMED_PARALLEL_MEM,$*_soapfuse,4,3G,5G,"$(PERL) $(SOAPFUSE) -c $(SOAPFUSE_CONFIG) -fd $(@D) -l $< -o $(@D)/$* && touch $@ && $(RM) -rf soapfuse/$*/alignWG")

soapfuse/sample_lists/%.txt : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(INIT) SAMPLEDIR=soapfuse/$*/$*; \
	mkdir -p $$SAMPLEDIR; ln -f -t $$SAMPLEDIR $^; \
	rename .1.fastq.gz _1.fastq.gz $$SAMPLEDIR/*.gz; \
	rename .2.fastq.gz _2.fastq.gz $$SAMPLEDIR/*.gz; \
	READLEN=`zcat $< | head -2 | sed 1d | wc -c`; \
	echo -e "$*\t$*\t$*\t$$READLEN" > $@

soapfuse/tables/%.sfuse.txt : soapfuse/%.timestamp
	$(INIT) cp -f soapfuse/$*/final_fusion_genes/$*/$*.final.Fusion.specific.for.genes $@

soapfuse/tables/%.isoform_sfuse.txt : soapfuse/%.timestamp
	$(INIT) cp -f soapfuse/$*/final_fusion_genes/$*/$*.final.Fusion.specific.for.trans $@

soapfuse/alltables/all.%.txt : $(foreach sample,$(SAMPLES),soapfuse/tables/$(sample).%.txt)
	$(INIT) head -1 $< | sed 's/^/Sample\t/;' > $@ && for i in $^; do sed "1d; s/^/$$(basename $${i%%.$*.txt})\t/" $$i >> $@; done

soapfuse/alltables/all.sfuse%coord.txt : soapfuse/alltables/all.sfuse%txt
	$(INIT) cut -f3,5,8,10 $< | awk 'BEGIN { OFS = "\t" } { print $$0, "$(ONCOFUSE_TISSUE_TYPE)" }' | sed '1d' > $@

soapfuse/alltables/all.isoform_sfuse%coord.txt : soapfuse/alltables/all.isoform_sfuse%txt
	$(INIT) cut -f4,7,11,14 $< | awk 'BEGIN { OFS = "\t" } { print $$0, "$(ONCOFUSE_TISSUE_TYPE)" }' | sed '1d' > $@

soapfuse/alltables/all.%.nft.txt : soapfuse/alltables/all.%.txt
	$(INIT) $(SOAPFUSE_NORMAL_FILTER) $(SOAPFUSE_NORMAL_FILTER_OPTS) $(NORMAL_SOAPFUSE_RESULTS) $< > $@

%.oncofuse.txt : %.coord.txt
	$(call LSCRIPT_MEM,8G,12G,"$(call ONCOFUSE_MEM,7G) $< coord $(ONCOFUSE_TISSUE_TYPE) $@")

%.oncofuse.merged.txt : %.txt %.oncofuse.txt 
	$(INIT) head -1 $< | sed 's/^/RowID\t/' > $<.tmp && awk 'BEGIN {OFS = "\t" } NR > 1 { print NR-1, $$0 }' $< >> $<.tmp ;\
		$(RSCRIPT) $(MERGE) -X --byColX 1 --byColY 1 -H $<.tmp $(<<) > $@ && rm -f $<.tmp

