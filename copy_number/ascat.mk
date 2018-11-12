include modules/Makefile.inc

LOGDIR ?= log/ascat.$(NOW)
PHONY += ascat ascat/log2 ascat/bafall ascat/bafhet ascat/mad ascat/log2nbaf ascat/ascat ascat/total

ascat : $(foreach pair,$(SAMPLE_PAIRS),ascat/log2/$(pair).pdf ascat/bafall/$(pair).pdf ascat/bafhet/$(pair).pdf ascat/mad/$(pair).RData ascat/log2nbaf/$(pair).pdf ascat/ascat/$(pair).pdf ascat/total/$(pair).pdf)

define ascat-plot-log2
ascat/log2/$1_$2.pdf : facets/cncf/$1_$2.Rdata
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 1G -m 2G,"$(RSCRIPT) modules/copy_number/ascat.R --type log2 --file_in $$< --file_out ascat/log2/$1_$2.pdf")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-plot-log2,$(tumor.$(pair)),$(normal.$(pair)))))
		
define ascat-plot-bafall
ascat/bafall/$1_$2.pdf : facets/cncf/$1_$2.Rdata
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 1G -m 2G,"$(RSCRIPT) modules/copy_number/ascat.R --type bafall --file_in $$< --file_out ascat/bafall/$1_$2.pdf")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-plot-bafall,$(tumor.$(pair)),$(normal.$(pair)))))

define ascat-plot-bafhet
ascat/bafhet/$1_$2.pdf : facets/cncf/$1_$2.Rdata
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 1G -m 2G,"$(RSCRIPT) modules/copy_number/ascat.R --type bafhet --file_in $$< --file_out ascat/bafhet/$1_$2.pdf")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-plot-bafhet,$(tumor.$(pair)),$(normal.$(pair)))))
		
define ascat-aspcf
ascat/mad/$1_$2.RData : facets/cncf/$1_$2.Rdata
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 3G -m 6G,"$(RSCRIPT) modules/copy_number/ascat.R --type aspcf --file_in $$< --file_out ascat/mad/$1_$2.RData --gamma '$${aspcf_gamma.$1}'")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-aspcf,$(tumor.$(pair)),$(normal.$(pair)))))

define ascat-plot-aspcf
ascat/log2nbaf/$1_$2.pdf : ascat/mad/$1_$2.RData
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 3G -m 6G,"$(RSCRIPT) modules/copy_number/ascat.R --type plot-aspcf --file_in $$< --file_out ascat/log2nbaf/$1_$2.pdf --nlog2 '$${aspacf_nlog2.$1}' --nbaf '$${aspacf_nbaf.$1}'")
	
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-plot-aspcf,$(tumor.$(pair)),$(normal.$(pair)))))
		
define ascat-run-ascat
ascat/ascat/$1_$2.pdf : ascat/mad/$1_$2.RData
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 3G -m 6G,"$(RSCRIPT) modules/copy_number/ascat.R --type run-ascat --file_in $$< --file_out ascat/ascat/$1_$2.pdf --rho '$${ascat_rho.$1}' --psi '$${ascat_psi.$1}' --nlog2 '$${aspacf_nlog2.$1}' --nbaf '$${aspacf_nbaf.$1}'")
	
ascat/total/$1_$2.pdf : facets/cncf/$1_$2.Rdata ascat/ascat/$1_$2.pdf
	$$(call RUN,-c -v ~/share/usr/anaconda-envs/ascat -s 6G -m 12G,"$(RSCRIPT) modules/copy_number/ascat.R --type total-copy --file_in $$< --file_out ascat/total/$1_$2.pdf")	

endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call ascat-run-ascat,$(tumor.$(pair)),$(normal.$(pair)))))
		
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
