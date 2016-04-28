# annotations that have pre-req annotations and are run last or should be run after the 2nd round of filtering

CLASSIFY_PATHOGENICITY = $(PYTHON) modules/vcf_tools/classify_pathogenicity_vcf.py
vcf/ann3/pathogen/%.pathogen.vcf : vcf/ann2/%.vcf
	$(INIT) $(CLASSIFY_PATHOGENICITY) $< $@

MUTATION_TASTER = $(PYTHON) modules/vcf_tools/mutation_taster_vcf.py
vcf/ann2/mut_taste/%.mut_taste.vcf : vcf/ft2/%.vcf
	$(INIT) $(MUTATION_TASTER) $< $@ &> $(LOG)

PROVEAN = $(RSCRIPT) modules/vcf_tools/proveanVcf.R
AA_TABLE = $(HOME)/share/reference/aa_table.tsv

PROVEAN_OPTS = --genome $(REF) --aaTable $(AA_TABLE) --ensemblTxdb $(ENSEMBL_TXDB) --mysqlHost $(EMBL_MYSQLDB_HOST) \
			   --mysqlPort $(EMBL_MYSQLDB_PORT) --mysqlUser $(EMBL_MYSQLDB_USER) --mysqlPassword $(EMBL_MYSQLDB_PW) \
			   --mysqlDb $(EMBL_MYSQLDB_DB) --numThreads 8 --memPerThread 1G --queue $(QUEUE) --qsubPriority $(QSUB_PRIORITY)

vcf/ann2/provean/%.provean.vcf : vcf/ft2/%.vcf
	$(call LSCRIPT_MEM,8G,10G,"$(PROVEAN) $(PROVEAN_OPTS) --outFile $@ $<")

CHASM = $(RSCRIPT) modules/vcf_tools/chasmVcf.R 
#CHASM_DIR = /ifs/opt/common/CHASM/CHASMDL.1.0.7
CHASM_DIR = $(HOME)/share/usr/CHASM
CHASM_PYTHON_ENV = $(HOME)/share/usr/anaconda-envs/pyenv27-chasm
CHASM_CLASSIFIER ?= Breast
vcf/ann2/chasm/%.chasm.vcf : vcf/ft2/%.vcf
	$(call CHECK_VCF,$(call LSCRIPT_CHECK_MEM,8G,17G,"unset PYTHONPATH && source $(CHASM_PYTHON_ENV)/bin/activate $(CHASM_PYTHON_ENV) && $(CHASM) --genome $(REF) --classifier $(subst $( ),$(,),$(CHASM_CLASSIFIER)) --chasmDir $(CHASM_DIR) --python $(shell which python) --outFile $@ $<"))

FATHMM = $(MY_RSCRIPT) modules/vcf_tools/fathmmVcf.R 
FATHMM_DIR = $(HOME)/share/usr/fathmm
FATHMM_PYTHON = $(HOME)/share/usr/bin/python
FATHMM_PYTHONPATH = $(HOME)/share/usr/lib/python:$(HOME)/share/usr/lib/python2.7
FATHMM_OPTS = --genome $(REF) --ensemblTxdb $(ENSEMBL_TXDB) --ref $(REF_FASTA) --fathmmDir $(FATHMM_DIR) --python $(FATHMM_PYTHON) \
			  --mysqlHost $(EMBL_MYSQLDB_HOST) --mysqlPort $(EMBL_MYSQLDB_PORT) --mysqlUser $(EMBL_MYSQLDB_USER) --mysqlPassword $(EMBL_MYSQLDB_PW) --mysqlDb $(EMBL_MYSQLDB_DB)

vcf/ann2/fathmm/%.fathmm.vcf : vcf/ft2/%.vcf
	$(call CHECK_VCF,$(call LSCRIPT_CHECK_MEM,8G,10G,"PYTHONPATH=$(FATHMM_PYTHONPATH) $(FATHMM) $(FATHMM_OPTS) \
		--outFile $@ $<"))

ANNOTATE_FACETS_VCF = $(RSCRIPT) modules/copy_number/annotateFacets2Vcf.R
define annotate-facets-pair
vcf/ann2/facets/$1.%.facets.vcf : vcf/ft2/$1.%.vcf facets/$1.cncf.txt
	$$(call LSCRIPT_MEM,4G,6G,"$$(ANNOTATE_FACETS_VCF) --facetsFile $$(<<) --outFile $$@ $$<")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call annotate-facets-pair,$(pair))))

define hrun-tumor-normal
vcf/ann2/hrun/$1_$2.%.hrun.vcf : vcf/ft2/$1_$2.%.vcf bam/$1.bam bam/$2.bam bam/$1.bai bam/$2.bai
	$$(call LSCRIPT_CHECK_PARALLEL_MEM,4,2G,3G,"$$(call GATK_MEM,8G) -T VariantAnnotator -nt 4 -R $$(REF_FASTA) -A HomopolymerRun --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -L $$< -o $$@")
endef
$(foreach pair,$(SAMPLE_PAIRS),$(eval $(call hrun-tumor-normal,$(tumor.$(pair)),$(normal.$(pair)))))

define hrun-sample
vcf/ann2/hrun/$1.%.hrun.vcf : vcf/ft2/$1.%.vcf bam/$1.bam bam/$1.bai
	$$(call LSCRIPT_CHECK_PARALLEL_MEM,4,2G,3G,"$$(call GATK_MEM,8G) -T VariantAnnotator -nt 4 -R $$(REF_FASTA) -L $$< -A HomopolymerRun --dbsnp $$(DBSNP) $$(foreach bam,$$(filter %.bam,$$^),-I $$(bam) ) -V $$< -o $$@")
endef
$(foreach sample,$(SAMPLES),$(eval $(call hrun-sample,$(sample))))

