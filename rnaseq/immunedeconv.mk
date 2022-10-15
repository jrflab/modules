include modules/Makefile.inc

LOGDIR = log/immunedeconv.$(NOW)

immunedeconv : immunedeconv/quantiseq.txt \
	       immunedeconv/mcpcounter.txt \
	       immunedeconv/cibersort.txt

immunedeconv/quantiseq.txt : kallisto/tpm_bygene.txt
	$(call RUN, -c -n 1 -s 8G -m 16G -v $(IMMUNE_ENV),"set -o pipefail && \
							   $(RSCRIPT) $(SCRIPTS_DIR)/immunedeconv.R --option 1 --input_file $(<) --output_file $(@)")

immunedeconv/mcpcounter.txt : kallisto/tpm_bygene.txt
	$(call RUN, -c -n 1 -s 8G -m 16G -v $(IMMUNE_ENV),"set -o pipefail && \
							   $(RSCRIPT) $(SCRIPTS_DIR)/immunedeconv.R --option 2 --input_file $(<) --output_file $(@)")

immunedeconv/cibersort.txt : kallisto/tpm_bygene.txt
	$(call RUN, -c -n 1 -s 8G -m 16G -v $(IMMUNE_ENV),"set -o pipefail && \
							   $(RSCRIPT) $(SCRIPTS_DIR)/immunedeconv.R --option 3 --input_file $(<) --output_file $(@)")

..DUMMY := $(shell mkdir -p version; \
	     ~/share/usr/env/r-immunedeconv-2.1.0/bin/R --version >> version/immunedeconv.txt;)
.SECONDARY:
.DELETE_ON_ERROR:
.PHONY: kallisto
