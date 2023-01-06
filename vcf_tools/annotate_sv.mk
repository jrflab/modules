include modules/Makefile.inc

LOGDIR ?= log/anotate_sv.$(NOW)

SV_CALLERS = svaba manta gridss merged
ANNOTATE_SV ?= $(HOME)/share/usr/env/annot_sv-3.1.3/opt/AnnotSV/bin/AnnotSV

annotate_sv :  $(foreach pair,$(SAMPLE_PAIRS), \
			$(foreach caller,$(SV_CALLERS),annotate_sv/$(pair)/$(pair).$(caller)_sv.tsv)) \
	       $(foreach pair,$(SAMPLE_PAIRS), \
			$(foreach caller,$(SV_CALLERS),annotate_sv/$(pair)/$(pair).$(caller)_sv.maf))
			
define annotate-sv
annotate_sv/$1/$2/$1.$2_sv.tsv : vcf/$1.$2_sv.vcf
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(ANNOTATE_SV_ENV),"set -o pipefail && \
							       mkdir -p annotate_sv/$1/$2 && \
							       $$(ANNOTATE_SV) \
							       -SVinputFile $$(<) \
							       -outputFile ./annotate_sv/$1/$2/$1.$2_sv.tsv \
							       -genomeBuild GRCh37")
							       
annotate_sv/$1/$1.$2_sv.tsv : annotate_sv/$1/$2/$1.$2_sv.tsv
	$$(INIT) cat $$(<) > $$(@)
	
annotate_sv/$1/$1.$2_sv.maf : vcf/$1.$2_sv.vcf
	$$(call RUN,-c -n 12 -s 1G -m 2G -v $(VEP_ENV),"set -o pipefail && \
							$$(VCF2MAF) \
							--input-vcf $$(<) \
							--tumor-id $1 \
							--filter-vcf $$(EXAC_NONTCGA) \
							--ref-fasta $$(REF_FASTA) \
							--vep-path $$(VEP_PATH) \
							--vep-data $$(VEP_DATA) \
							--tmp-dir `mktemp -d` \
							--output-maf $$(@)")

endef
$(foreach pair,$(SAMPLE_PAIRS),\
	$(foreach caller,$(SV_CALLERS), \
		$(eval $(call annotate-sv,$(pair),$(caller)))))
		
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: annotate_sv
