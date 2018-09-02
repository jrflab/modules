include modules/Makefile.inc

LOGDIR ?= log/hla_optitype.$(NOW)
PHONY += hla_optitype 

hla_optitype : $(foreach pair,$(SAMPLE_PAIRS),hla_optitype/$(pair)/$(pair).taskcomplete)

define hla-optitype
hla_optitype/$1_$2/winners.hla.txt : bam/$1.bam bam/$2.bam
	$$(call RUN,-n 8 -s 12G -m 24G, "source /home/${USER}/share/usr/anaconda-envs/jrflab-modules-0.1.5/bin/activate \
									 /home/${USER}/share/usr/anaconda-envs/optitype && \
								 	 if [ ! -d hla_optitype/$1_$2 ]; then mkdir hla_optitype/$1_$2; fi && \
								 	 ")
								 	  
hla_optitype/$1_$2/$1_$2.taskcomplete : hla_optitype/$1_$2/
	$$(call RUN,-s 1G -m 1G,"touch hla_optitype/$1_$2/$1_$2.taskcomplete")
endef
$(foreach pair,$(SAMPLE_PAIRS),\
		$(eval $(call hla-optitype,$(tumor.$(pair)),$(normal.$(pair)))))

.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
