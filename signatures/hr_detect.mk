include modules/Makefile.inc

LOGDIR ?= log/hr_detect.$(NOW)

MIN_SIZE = 1
MAX_SIZE = 100000000000000000000

hr_detect :  $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).merged.bed) \
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).merged.bedpe) \
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).snv.vcf) \
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).indel.vcf) \
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).cn.txt)

define hr-detect
hr_detect/$1_$2/$1_$2.merged.bed : vcf/$1_$2.merged_sv.vcf
	$$(call RUN,-c -n 1 -s 4G -m 8G -v $(SURVIVOR_ENV),"set -o pipefail && \
							    SURVIVOR vcftobed \
							    $$(<) \
							    $(MIN_SIZE) \
							    $(MAX_SIZE) \
							    $$(@)")
							    
hr_detect/$1_$2/$1_$2.merged.bedpe : hr_detect/$1_$2/$1_$2.merged.bed
	$$(call RUN,-c -n 1 -s 4G -m 8G,"set -o pipefail && \
					 echo \"chrom1	start1	end1	chrom2	start2	end2	sv_id	pe_support	strand1	strand2	svclass\" > \
					 $$(@) && \
					 cat $$(<) >> $$(@)")
					 
hr_detect/$1_$2/$1_$2.snv.vcf : summary/tsv/all.tsv
	$$(call RUN,-c -n 1 -s 12G -m 16G,"set -o pipefail && \
					   $(RSCRIPT) modules/scripts/hr_detect.R \
					   --option 1 \
					   --sample_name $1_$2")

hr_detect/$1_$2/$1_$2.indel.vcf : summary/tsv/all.tsv
	$$(call RUN,-c -n 1 -s 12G -m 16G,"set -o pipefail && \
					   $(RSCRIPT) modules/scripts/hr_detect.R \
					   --option 2 \
					   --sample_name $1_$2")

hr_detect/$1_$2/$1_$2.cn.txt : facets/cncf/$1_$2.txt
	$$(call RUN,-c -n 1 -s 12G -m 16G,"set -o pipefail && \
					   $(RSCRIPT) modules/scripts/hr_detect.R \
					   --option 3 \
					   --sample_name $1_$2")

endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call hr-detect,$(tumor.$(pair)),$(normal.$(pair)))))
		
..DUMMY := $(shell mkdir -p version; \
	     R --version &> version/hr_detect.txt;)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: hr_detect
