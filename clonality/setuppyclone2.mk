include modules/Makefile.inc

LOGDIR ?= log/setup_pyclone2.$(NOW)
PHONY += pyclone

setup_pyclone2 : $(foreach pair,$(SAMPLE_PAIRS),pyclone/$(pair)/$(pair).yaml)

define make-input-pyclone
pyclone/$1_$2/$1_$2.tsv : summary/tsv/mutation_summary.tsv ascat/ascat/$1_$2.RData
	$$(call RUN,-c -s 4G -m 6G,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
								if [ ! -d pyclone/$1_$2 ]; then mkdir pyclone/$1_$2; fi && \
								$(RSCRIPT) modules/clonality/tsvforpyclone2.R --tumor_name $1 --normal_name $2")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call make-input-pyclone,$(tumor.$(pair)),$(normal.$(pair)))))

define build-mutations-file
pyclone/$1_$2/$1_$2.yaml : pyclone/$1_$2/$1_$2.tsv
	$$(call RUN,-c -s 4G -m 6G,"source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate /home/${USER}/share/usr/anaconda-envs/PyClone-0.13.1 && \
								PyClone build_mutations_file --in_file pyclone/$1_$2/$1_$2.tsv --out_file pyclone/$1_$2/$1_$2.yaml --prior total_copy_number")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call build-mutations-file,$(tumor.$(pair)),$(normal.$(pair)))))

