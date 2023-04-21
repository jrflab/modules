include modules/Makefile.inc

LOGDIR ?= log/hr_detect.$(NOW)

MIN_SIZE = 1
MAX_SIZE = 100000000000000000000

hr_detect :  $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).merged.bed) \
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).merged.bedpe) \
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).snv.vcf) \
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).snv.vcf.bgz) \
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).snv.vcf.bgz.tbi) \
	     
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).snv_repaired.vcf) \
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).snv_repaired.vcf.bgz) \
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).snv_repaired.vcf.bgz.tbi) \
	     
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).indel.vcf) \
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).indel.vcf.bgz) \
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).indel.vcf.bgz.tbi) \
	     
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).indel_repaired.vcf) \
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).indel_repaired.vcf.bgz) \
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).indel_repaired.vcf.bgz.tbi) \
	     
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).cn.txt) \
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).sv.bedpe) \
	     hr_detect/hrdetect.txt \
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).png) \
	     $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).svg)

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
	$$(call RUN,-c -n 1 -s 12G -m 16G -v $(SIGNATURE_TOOLS_ENV),"set -o pipefail && \
					   			     $(RSCRIPT) modules/scripts/hr_detect.R \
								     --option 1 \
								     --sample_name $1_$2")

hr_detect/$1_$2/$1_$2.snv.vcf.bgz : hr_detect/$1_$2/$1_$2.snv.vcf
	$$(call RUN,-c -n 1 -s 12G -m 16G -v $(INNOVATION_ENV),"set -o pipefail && \
								bgzip -c $$(<) > $$(@)")

hr_detect/$1_$2/$1_$2.snv.vcf.bgz.tbi : hr_detect/$1_$2/$1_$2.snv.vcf.bgz
	$$(call RUN,-c -n 1 -s 12G -m 16G -v $(INNOVATION_ENV),"set -o pipefail && \
								tabix -p vcf $$(<)")
								
hr_detect/$1_$2/$1_$2.indel.vcf : summary/tsv/all.tsv
	$$(call RUN,-c -n 1 -s 12G -m 16G -v $(SIGNATURE_TOOLS_ENV),"set -o pipefail && \
								     $(RSCRIPT) modules/scripts/hr_detect.R \
								     --option 2 \
								     --sample_name $1_$2")
								     
hr_detect/$1_$2/$1_$2.indel.vcf.bgz : hr_detect/$1_$2/$1_$2.indel.vcf
	$$(call RUN,-c -n 1 -s 12G -m 16G -v $(INNOVATION_ENV),"set -o pipefail && \
								bgzip -c $$(<) > $$(@)")

hr_detect/$1_$2/$1_$2.indel.vcf.bgz.tbi : hr_detect/$1_$2/$1_$2.indel.vcf.bgz
	$$(call RUN,-c -n 1 -s 12G -m 16G -v $(INNOVATION_ENV),"set -o pipefail && \
								tabix -p vcf $$(<)")
								
hr_detect/$1_$2/$1_$2.snv_repaired.vcf : hr_detect/$1_$2/$1_$2.snv.vcf.bgz
	$$(call RUN,-c -n 1 -s 12G -m 16G -v $(INNOVATION_ENV),"set -o pipefail && \
								bcftools view $$(<) > $$(@)")
								
hr_detect/$1_$2/$1_$2.snv_repaired.vcf.bgz : hr_detect/$1_$2/$1_$2.snv_repaired.vcf
	$$(call RUN,-c -n 1 -s 12G -m 16G -v $(INNOVATION_ENV),"set -o pipefail && \
								bgzip -c $$(<) > $$(@)")

hr_detect/$1_$2/$1_$2.snv_repaired.vcf.bgz.tbi : hr_detect/$1_$2/$1_$2.snv_repaired.vcf.bgz
	$$(call RUN,-c -n 1 -s 12G -m 16G -v $(INNOVATION_ENV),"set -o pipefail && \
								tabix -p vcf $$(<)")

hr_detect/$1_$2/$1_$2.indel_repaired.vcf : hr_detect/$1_$2/$1_$2.indel.vcf.bgz
	$$(call RUN,-c -n 1 -s 12G -m 16G -v $(INNOVATION_ENV),"set -o pipefail && \
								bcftools view $$(<) > $$(@)")
								

hr_detect/$1_$2/$1_$2.indel_repaired.vcf.bgz : hr_detect/$1_$2/$1_$2.indel_repaired.vcf
	$$(call RUN,-c -n 1 -s 12G -m 16G -v $(INNOVATION_ENV),"set -o pipefail && \
								bgzip -c $$(<) > $$(@)")

hr_detect/$1_$2/$1_$2.indel_repaired.vcf.bgz.tbi : hr_detect/$1_$2/$1_$2.indel_repaired.vcf.bgz
	$$(call RUN,-c -n 1 -s 12G -m 16G -v $(INNOVATION_ENV),"set -o pipefail && \
								tabix -p vcf $$(<)")

hr_detect/$1_$2/$1_$2.cn.txt : facets/cncf/$1_$2.txt
	$$(call RUN,-c -n 1 -s 12G -m 16G -v $(SIGNATURE_TOOLS_ENV),"set -o pipefail && \
								     $(RSCRIPT) modules/scripts/hr_detect.R \
								     --option 3 \
								     --sample_name $1_$2")
								     
hr_detect/$1_$2/$1_$2.sv.bedpe : hr_detect/$1_$2/$1_$2.merged.bedpe
	$$(call RUN,-c -n 1 -s 12G -m 16G -v $(SIGNATURE_TOOLS_ENV),"set -o pipefail && \
								     $(RSCRIPT) modules/scripts/hr_detect.R \
								     --option 4 \
								     --sample_name $1_$2")
								     
hr_detect/$1_$2/$1_$2.png : hr_detect/$1_$2/$1_$2.snv_repaired.vcf.bgz.tbi hr_detect/$1_$2/$1_$2.indel_repaired.vcf.bgz.tbi hr_detect/$1_$2/$1_$2.cn.txt hr_detect/$1_$2/$1_$2.sv.bedpe
	$$(call RUN,-c -n 1 -s 12G -m 16G -v $(SIGNATURE_TOOLS_ENV),"set -o pipefail && \
								    $(RSCRIPT) modules/scripts/hr_detect.R \
								    --option 5 \
								    --sample_name $1_$2 && \
								    mv hr_detect/$1_$2/$1_$2.genomePlot.png $$(@)")

hr_detect/$1_$2/$1_$2.svg : hr_detect/$1_$2/$1_$2.snv_repaired.vcf.bgz.tbi hr_detect/$1_$2/$1_$2.indel_repaired.vcf.bgz.tbi hr_detect/$1_$2/$1_$2.cn.txt hr_detect/$1_$2/$1_$2.sv.bedpe
	$$(call RUN,-c -n 1 -s 12G -m 16G -v $(SIGNATURE_TOOLS_ENV),"set -o pipefail && \
								     $(RSCRIPT) modules/scripts/hr_detect.R \
								     --option 6 \
								     --sample_name $1_$2 && \
								     mv hr_detect/$1_$2/$1_$2.genomePlot.svg $$(@)")

endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call hr-detect,$(tumor.$(pair)),$(normal.$(pair)))))
		
hr_detect/hrdetect.txt : $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).sv.bedpe) $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).snv_repaired.vcf) $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).indel_repaired.vcf.bgz) $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).indel_repaired.vcf.bgz.tbi) $(foreach pair,$(SAMPLE_PAIRS),hr_detect/$(pair)/$(pair).cn.txt)
	$(call RUN, -c -n 4 -s 6G -m 9G -v $(SIGNATURE_TOOLS_ENV),"set -o pipefail && \
					  			     $(RSCRIPT) modules/scripts/hr_detect.R --option 7 --sample_name '$(SAMPLE_PAIRS)'")

		
..DUMMY := $(shell mkdir -p version; \
	     R --version &> version/hr_detect.txt;)
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: hr_detect
