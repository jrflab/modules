include modules/Makefile.inc

LOGDIR ?= log/pratchet.$(NOW)
PHONY += pratchet

pratchet : $(foreach set,$(SAMPLE_SETS),pratchet/$(set)/tree_final.RData) $(foreach set,$(SAMPLE_SETS),pratchet/$(set)/tree_final.pdf)

define parsimony-ratchet
pratchet/%/tree_final.RData : sufam/%.tsv
	$$(call RUN,-c -s 8G -m 12G -v $(PHANGORN_ENV),"mkdir -p pratchet && \
								 					mkdir -p pratchet/$$* && \
								 					$(RSCRIPT) modules/test/phylogeny/pratchet.R --sample_set $$* --normal_samples '$(NORMAL_SAMPLES)'")

pratchet/%/tree_final.pdf : pratchet/%/tree_final.RData
	$$(call RUN,-c -s 4G -m 6G -v $(PHYLO_ENV),"$(RSCRIPT) modules/test/phylogeny/plotratchet.R --sample_set $$(*)")

endef
$(foreach set,$(SAMPLE_SETS),\
		$(eval $(call parsimony-ratchet,$(set))))
		

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
