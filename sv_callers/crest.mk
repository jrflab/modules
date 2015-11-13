# run crest

include modules/Makefile.inc

LOGDIR = log/crest.$(NOW)

CREST_DIR = $(HOME)/share/usr/crest
CREST = PERL5LIB=/home/limr/share/usr/crest:/home/limr/share/usr/src/bioperl-live:/home/limr/share/usr/perl5/lib/perl5:/home/limr/share/usr/src/ensembl/modules:/home/limr/share/usr/src/ensembl-compara/modules:/home/limr/share/usr/src/ensembl-variation/modules:/home/limr/share/usr/src/ensembl-funcgen/modules $(PERL) $(CREST_DIR)/CREST.pl
CREST_OPTS = --blat $(BLAT) --cap3 $(CAP3) --blatclient $(GFCLIENT) --blatserver 140.163.153.48 --blatport 88878 -t $(REF_2BIT) --ref_genome $(REF_FASTA) 
EXTRACT_SCLIP = PERL5LIB=$(CREST_DIR):$(PERL5LIB) $(CREST_DIR)/extractSClip.pl
EXTRACT_SCLIP_OPTS = --ref_genome $(REF_FASTA)

.SECONDARY:
.DELETE_ON_ERROR:


ifdef SAMPLE_PAIRS
PHONY += crestTN
crestTN : $(foreach pair,$(SAMPLE_PAIRS),crest/$(pair).crest_timestamp)
else
PHONY += crest
crest : $(foreach sample,$(SAMPLES),crest/$(sample).crest_timestamp)
endif

crest/%.read_len : bam/%.bam
	$(call LSCRIPT_MEM,3G,5G,"$(SAMTOOLS) view $< | tail -n+100000 | head -1 | awk '{ print length(\$$10) }' > $@")

crest/%.sclip.txt : bam/%.bam
	$(call LSCRIPT_MEM,6G,8G,"$(EXTRACT_SCLIP) $(EXTRACT_SCLIP_OPTS) -p $(@D)/$* -i $<")

crest/%.crest_timestamp : bam/%.bam crest/%.sclip.txt crest/%.read_len
	$(call LSCRIPT_MEM,15G,60G,"$(CREST) $(CREST_OPTS) -f $(<<) -d $< -p $(@D)/$* -l `cat $(<<<)` && touch $@")

define crest-tumor-normal
crest/$1_$2.crest_timestamp : bam/$1.bam bam/$2.bam crest/$1.sclip.txt crest/$1.read_len
	$$(call LSCRIPT_MEM,15G,60G,"$$(CREST) $$(CREST_OPTS) -f $$(<<<) -d $$< -g $$(<<) -p $$(@D)/$1_$2 -l `cat $$(<<<<)` && touch $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call crest-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

.PHONY: $(PHONY)
