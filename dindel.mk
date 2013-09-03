include ~/share/modules/Makefile.inc
include ~/share/modules/hg19.inc

SAMPLE_FILE ?= samples.txt
SAMPLES ?= $(shell cat $(SAMPLE_FILE))

DINDEL = $(HOME)/share/usr/bin/dindel
MAKE_WINDOWS = $(HOME)/share/usr/bin/makeWindows.py
MERGE_OUTPUT = $(HOME)/share/usr/bin/mergeOutputDiploid.py
VCF_SORT = $(PERL) $(HOME)/share/usr/bin/vcfsorter.pl $(UCSC_REF_DICT)
NUM_WINDOWS_PER_FILE = 1000

.PHONY : windows vcf
.SECONDARY:
.DELETE_ON_ERROR:

windows : $(foreach sample,$(SAMPLES),dindel/$(sample).realn_windows_timestamp)
vcf : $(foreach sample,$(SAMPLES),dindel/vcf/$(sample).dindel.sorted.annotated.vcf)

# extract candidate indels and insert size dist
dindel/candidate/%.variants.txt dindel/candidate/%.libraries.txt : %.bam %.bam.bai
	$(call LSCRIPT_MEM,8G,12G,"$(DINDEL) --analysis getCIGARindels --bamFile $< --outputFile dindel/candidate/$* --ref $(REF_FASTA) &> $(LOGDIR)/$(@F).log")

# make realignment windows
dindel/%.realn_windows_timestamp : dindel/candidate/%.variants.txt
	$(call LSCRIPT_MEM,8G,12G,"$(MKDIR) dindel/windows; $(MAKE_WINDOWS) --inputVarFile $< --windowFilePrefix dindel/windows/$*.realn_windows --numWindowsPerFile $(NUM_WINDOWS_PER_FILE) &> $(LOGDIR)/$*.realn_windows.log && time > $@")

# realignment
define find-indels-window
dindel/windows/$1.variants.%.glf.txt : $1.bam dindel/windows/$1.realn_windows.%.txt dindel/candidate/$1.libraries.txt $1.bam.bai
	$$(call INIT_MEM,5G,10G) $$(DINDEL) --analysis indels --doDiploid --bamFile $$< --ref $$(REF_FASTA) --varFile $$(word 2,$$^) --libFile $$(word 3,$$^) --outputFile $$(@:.glf.txt=) &> $$(LOGDIR)/$$(@F).log
endef
$(foreach sample,$(SAMPLES),$(eval $(call find-indels-window,$(sample))))

define file-list
dindel/$1.variants_filelist.txt : $(subst .txt,.glf.txt,$(subst realn_windows,variants,$(wildcard dindel/windows/$1.realn_windows.*.txt)))
	$$(INIT) ls dindel/windows/$1.variants.*.txt > $$(@)
endef
$(foreach sample,$(SAMPLES),$(eval $(call file-list,$(sample))))

dindel/vcf/%.dindel.vcf : dindel/%.variants_filelist.txt
	$(call INIT_MEM,5G,10G) $(MERGE_OUTPUT) --inputFiles $< --outputFile $@ --ref $(REF_FASTA) > $(LOGDIR)/$(@F).log

%.sorted.vcf : %.vcf
	$$(INIT) $(VCF_SORT) $< | grep -v '1/2' > $@

%.vcf.idx : %.vcf
	$$(INIT) $(VCFTOOLS) index $<

