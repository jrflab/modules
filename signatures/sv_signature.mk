include modules/Makefile.inc

LOGDIR ?= log/sv_signature.$(NOW)

MIN_SIZE = 1
MAX_SIZE = 100000000000000000000

signature_sv :  $(foreach pair,$(SAMPLE_PAIRS),sv_signature/$(pair)/$(pair).merged.bed) \
		$(foreach pair,$(SAMPLE_PAIRS),sv_signature/$(pair)/$(pair).merged.bedpe) \
		$(foreach pair,$(SAMPLE_PAIRS),sv_signature/$(pair)/$(pair).merged_exposures.txt) \
		sv_signature/exposures.txt
		
define signature-sv
sv_signature/$1_$2/$1_$2.merged.bed : vcf/$1_$2.merged_sv.vcf
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(SURVIVOR_ENV),"set -o pipefail && \
							    SURVIVOR vcftobed \
							    $$(<) \
							    $(MIN_SIZE) \
							    $(MAX_SIZE) \
							    $$(@)")
							    
sv_signature/$1_$2/$1_$2.merged.bedpe : sv_signature/$1_$2/$1_$2.merged.bed
	$$(call RUN,-c -n 1 -s 4G -m 8G,"set -o pipefail && \
					 echo \"chrom1	start1	end1	chrom2	start2	end2	sv_id	pe_support	strand1	strand2	svclass\" > \
					 $$(@) && \
					 cat $$(<) >> $$(@)")
					 
sv_signature/$1_$2/$1_$2.merged_exposures.txt : sv_signature/$1_$2/$1_$2.merged.bedpe
	$$(call RUN,-c -n 4 -s 2G -m 4G -v $(SIGNATURE_TOOLS_ENV),"set -o pipefail && \
								  $(RSCRIPT) $(SCRIPTS_DIR)/sv_signature.R \
								  --option 1 \
								  --sample_name $1_$2 \
								  --input_file $$(<) \
								  --output_file sv_signature/$1_$2/$1_$2.merged")

endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call signature-sv,$(tumor.$(pair)),$(normal.$(pair)))))
		
sv_signature/exposures.txt : $(foreach pair,$(SAMPLE_PAIRS),sv_signature/$(pair)/$(pair).merged_exposures.txt)
	$(call RUN, -c -n 1 -s 8G -m 12G -v $(SIGNATURE_TOOLS_ENV),"set -o pipefail && \
					  			    $(RSCRIPT) $(SCRIPTS_DIR)/sv_signature.R --option 2 --sample_name '$(SAMPLE_PAIRS)' --output_file $(@)")


..DUMMY := $(shell mkdir -p version; \
	     $(SURVIVOR_ENV)/bin/SURVIVOR --version &> version/sv_signature.txt;)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: signature_sv
