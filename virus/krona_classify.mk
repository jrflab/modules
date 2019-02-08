include modules/Makefile.inc

LOGDIR ?= log/krona_classify.$(NOW)
PHONY += unmapped_reads

krona_classify : $(foreach sample,$(SAMPLES),unmapped_reads/$(sample).html)

define krona-classify
unmapped_reads/%.html : unmapped_reads/%.blast
	$(call RUN,-n 1 -s 4G -m 9G,"ktClassifyBLAST -s $$< -o unmapped_reads/$$*.tax && ktImportTaxonomy -m 1 unmapped_reads/$$*.tax -o unmapped_reads/$$*.html")
endef
$(foreach sample,$(SAMPLES),\
		$(eval $(call krona-classify,$(sample))))


.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
