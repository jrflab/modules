# soapfuse
# vim: set ft=make :

include ~/share/modules/Makefile.inc
include ~/share/modules/gatk.inc

LOGDIR = log/soapfuse.$(NOW)

SOAPFUSE = $(HOME)/usr/SOAPfuse-v1.26/SOAPfuse-RUN.pl 
SOAPFUSE_CONFIG = $(HOME)/share/usr/SOAPfuse-v1.26/config/config.txt
PREPARE_SOAPFUSE = $(HOME)/share/scripts/prepareSoapFuse.pl

ONCOFUSE_MEM = $(JAVA) -Xmx$1 -jar $(HOME)/share/usr/oncofuse-v1.0.3/Oncofuse.jar
ONCOFUSE_TISSUE_TYPE ?= EPI

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: all sample_tables soapfuse

all : soapfuse/alltables/all.sfuse.oncofuse.merged.txt soapfuse/alltables/all.sfuse_isoforms.txt
sample_tables : $(foreach sample,$(SAMPLES),soapfuse/tables/$(sample).sfuse.txt soapfuse/tables/$(sample).sfuse_isoforms.txt)
soapfuse : $(foreach sample,$(SAMPLES),soapfuse/$(sample).timestamp)

soapfuse/%.timestamp : soapfuse/sample_lists/%.txt
	$(call LSCRIPT_NAMED_PARALLEL_MEM,$*_soapfuse,4,2G,2.5G,"$(SOAPFUSE) -c $(SOAPFUSE_CONFIG) -fd $(@D) -l $< -o $(@D)/$* && touch $@ && $(RM) -rf soapfuse/$*/alignWG")

soapfuse/sample_lists/%.txt : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(INIT) SAMPLEDIR=soapfuse/$*/$*; \
	mkdir -p $$SAMPLEDIR; ln -f -t $$SAMPLEDIR $^; \
	rename .1.fastq.gz _1.fastq.gz $$SAMPLEDIR/*.gz; \
	rename .2.fastq.gz _2.fastq.gz $$SAMPLEDIR/*.gz; \
	READLEN=`zcat $< | head -2 | sed 1d | wc -c`; \
	echo -e "$*\t$*\t$*\t$$READLEN" > $@

soapfuse/tables/%.sfuse.txt : soapfuse/%.timestamp
	$(INIT) ln -f soapfuse/$*/final_fusion_genes/$*/$*.final.Fusion.specific.for.genes $@

soapfuse/tables/%.sfuse_isoforms.txt : soapfuse/%.timestamp
	$(INIT) ln -f soapfuse/$*/final_fusion_genes/$*/$*.final.Fusion.specific.for.trans $@

soapfuse/alltables/all.%.txt : $(foreach sample,$(SAMPLES),soapfuse/tables/$(sample).%.txt)
	$(INIT) head -1 $< | sed 's/^/Sample\t/;' > $@ && for i in $^; do sed "1d; s/^/$$(basename $${i%%.sfuse.$*.txt})\t/" $$i >> $@; done

soapfuse/alltables/all.sfuse.coord.txt : soapfuse/alltables/all.sfuse.txt
	$(INIT) cut -f3,5,8,10 $< | awk 'BEGIN { OFS = "\t" } { print $$0, "$(ONCOFUSE_TISSUE_TYPE)" }' | sed '1d' > $@

soapfuse/alltables/all.sfuse_isoforms.coord.txt : soapfuse/alltables/all.sfuse_isforms.txt
	$(INIT) cut -f4,6,11,13 $< | awk 'BEGIN { OFS = "\t" } { print $$0, "$(ONCOFUSE_TISSUE_TYPE)" }' | sed '1d' > $@

%.oncofuse.txt : %.coord.txt
	$(call LSCRIPT_MEM,8G,12G,"$(call ONCOFUSE_MEM,7G) $< coord $(ONCOFUSE_TISSUE_TYPE) $@")

%.oncofuse.merged.txt : %.oncofuse.txt %.txt
	$(INIT) head -1 $< | sed 's/^/RowID\t/' > $<.tmp && awk 'BEGIN {OFS = "\t" } NR > 1 { print NR-1, $$0 }' $< >> $<.tmp ;\
		$(RSCRIPT) $(MERGE) -X --byColX 1 --byColY 1 -H $<.tmp $(word 2,$^) > $@

