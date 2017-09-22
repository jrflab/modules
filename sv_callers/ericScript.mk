# run EricScript
# author: Raymond Lim

ERICSCRIPT_ENV = $(HOME)/share/usr/anaconda-envs/ericscript-0.5.5/
ERICSCRIPT = ericscript.pl
ERICSCRIPT_OPTS ?= --refid $(ERICSCRIPT_SPECIES)  -db $(ERICSCRIPT_DB) --remove
ERICSCRIPT_TO_USV = python modules/sv_callers/ericscript2usv.py

.PHONY: ericscript
.SECONDARY:
.DELETE_ON_ERROR:

ericscript: $(foreach sample,$(SAMPLES),usv/$(sample).ericscript.tsv)

ericscript/%_ericscript.timestamp : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$(call RUN,-s 14G -m 14G -n 7 -N $*_ericscript,"$(ERICSCRIPT) $(ERICSCRIPT_OPTS) -p 7 -name $* -o $(@D)/$* $^ && \
		touch $@")

usv/%.ericscript.tsv : ericscript/%_ericscript.timestamp
	$(call RUN,,"$(ERICSCRIPT_TO_USV) < eriscript/$*/$*.results.filtered.tsv > $@")
