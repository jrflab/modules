include modules/Makefile.inc

LOGDIR ?= log/krona_classify.$(NOW)

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: krona_classify

VPATH = unmapped_reads

krona_classify : $(foreach sample,$(SAMPLES),unmapped_reads/$(sample).html)

unmapped_reads/%.html : unmapped_reads/%.blast
	$(call RUN,-n 1 -s 4G -m 9G,"ktClassifyBLAST -s $< -o unmapped_reads/$*.tax && ktImportTaxonomy -m 1 unmapped_reads/$*.tax -o unmapped_reads/$*.html")
