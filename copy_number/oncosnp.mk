# OncoSNP module
include modules/Makefile.inc

SAMPLE_FILE = samples.txt
SAMPLE_SEX_FILE = samples.sex.txt
SAMPLES = $(shell cat $(SAMPLE_FILE))

LOGDIR = oncosnp/log
GENERATE_GENO_CLUSTER ?= true

SNP6_CDF = $(HOME)/share/references/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.Full.cdf
SNP6_MODELS = $(HOME)/share/references/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.birdseed.models
SNP6_SPECIAL = $(HOME)/share/references/GenomeWideSNP_6/LibFiles/GenomeWideSNP_6.Full.specialSNPs
APT_PROBESET_GENOTYPE = $(HOME)/share/usr/apt/bin/apt-probeset-genotype
GENOTYPE_PATHWAY = birdseed
APT_PROBESET_SUMMARIZE = $(HOME)/share/usr/apt/bin/apt-probeset-summarize
SUMMARIZE_PATHWAY = quant-norm.sketch=50000,pm-only,med-polish,expr.genotype=true
SUMMARIZE_TARGET = $(HOME)/share/usr/penncnv/gw6/lib/hapmap.quant-norm.normalization-target.txt
PENNCNV_NORMALIZE = $(HOME)/share/usr/penncnv/gw6/bin/normalize_affy_geno_cluster.pl
PENNCNV_GENERATE_AFFY_GENO_CLUSTER = $(HOME)/share/usr/penncnv/gw6/bin/generate_affy_geno_cluster.pl

ifeq ($(GENERATE_GENO_CLUSTER),true) 
HAPMAP_GENOCLUSTER = oncosnp/penncnv/gw6.genocluster
else
HAPMAP_GENOCLUSTER = $(HOME)/share/usr/penncnv/gw6/lib/hapmap.genocluster
endif

#GW6_LOC = $(HOME)/share/usr/penncnv/gw6/lib/affygw6.hg18.pfb
GW6_LOC = $(HOME)/share/usr/penncnv/gw6/lib/affygw6.hg19.pfb
COLUMN_EXTRACT = $(HOME)/share/usr/penncnv/kcolumn.pl
ONCOSNP = $(HOME)/share/usr/oncosnp/run_oncosnp.sh $(HOME)/share/usr/matlab/mcr_v714/v714
ONCOSNP_GC = $(HOME)/share/references/quantisnp2/b37
ONCOSNP_PARAMS = $(HOME)/share/usr/oncosnp/configuration/hyperparameters-affy.dat
ONCOSNP_LEVELS = $(HOME)/share/usr/oncosnp/configuration/levels-affy.dat
ONCOSNP_TRAINING_STATES = $(HOME)/share/usr/oncosnp/configuration/trainingStates.dat
ONCOSNP_TUMOUR_STATES = $(HOME)/share/usr/oncosnp/configuration/tumourStates.dat

.SECONDARY:
.DELETE_ON_ERROR: 
.PHONY : all

all : $(foreach sample,${SAMPLES},oncosnp/oncosnp.$(sample).timestamp)

# Generate genotyping calls from CEL files
oncosnp/apt/birdseed.calls.txt : $(foreach sample,$(SAMPLES),cel/$(sample).CEL)
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,8G,10G)" $(MKDIR) $(@D) $(LOGDIR); \
	$(APT_PROBESET_GENOTYPE) -c $(SNP6_CDF) -a $(GENOTYPE_PATHWAY) --read-models-birdseed $(SNP6_MODELS) --special-snps $(SNP6_SPECIAL) --out-dir $(@D) $^ &> $(LOGDIR)/$(@F).log

# Allele-specific signal extraction from CEL files
oncosnp/apt/quant-norm.pm-only.med-polish.expr.summary.txt : $(foreach sample,$(SAMPLES),cel/$(sample).CEL)
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,8G,10G)" $(MKDIR) $(@D) $(LOGDIR); \
	$(APT_PROBESET_SUMMARIZE) -d $(SNP6_CDF) -a $(SUMMARIZE_PATHWAY) --target-sketch $(SUMMARIZE_TARGET) --out-dir $(@D) $^ &> $(LOGDIR)/$(@F).log

# Generate canonical genotype clustering file
# This step is not necessary because you can use the default geno-cluster file provided by penncnv. 
# But you can can generate more accurate calls if you generate it, but this should only be done if you have more than a few dozen .CEL files
oncosnp/penncnv/gw6.genocluster : oncosnp/apt/birdseed.calls.txt oncosnp/apt/quant-norm.pm-only.med-polish.expr.summary.txt
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,8G,10G)" $(MKDIR) $(@D) $(LOGDIR); \
	${PENNCNV_GENERATE_AFFY_GENO_CLUSTER} $(word 1,$^) oncosnp/apt/birdseed.confidences.txt $(word 2,$^) -locfile ${GW6_LOC} -sexfile ${SAMPLE_SEX_FILE} -out $@ &> $(LOGDIR)/$(@F).log

# Log R and B allelic freq calculation
oncosnp/penncnv/logR_BAF.txt : oncosnp/apt/quant-norm.pm-only.med-polish.expr.summary.txt
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,8G,10G)" $(MKDIR) $(@D) $(LOGDIR); \
	$(PENNCNV_NORMALIZE) $(HAPMAP_GENOCLUSTER) $< -locfile $(GW6_LOC) -out $@ &> $(LOGDIR)/$(@F).log

# split log R and B AFs
oncosnp/penncnv/logR_BAF_split.timestamp : oncosnp/penncnv/logR_BAF.txt
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,8G,10G)" rm -fr oncosnp/penncnv/logR_BAF; $(MKDIR) oncosnp/penncnv/logR_BAF; $(COLUMN_EXTRACT) $< split 2 -tab -head 3 -name -beforestring .CEL -out oncosnp/penncnv/logR_BAF/gw6 &> $(LOGDIR)/logR_BAF_split.log && touch $@

# create batchfile for oncosnp
oncosnp/%.batchfile : oncosnp/penncnv/logR_BAF_split.timestamp
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,1G,2G)" \
	echo -e "sampleid\ttumourfile\tnormalfile" > $@; \
	echo -e "$*\toncosnp/penncnv/logR_BAF/gw6.$*\t" >> $@
#	sed 's;^;oncosnp/penncnv/logR_BAF/gw6.;' $(SAMPLE_FILE) | paste $(SAMPLE_FILE) - >> $@

# write oncosnp results to results dir
oncosnp/oncosnp.%.timestamp : oncosnp/%.batchfile
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,8G,10G)" \
	$(ONCOSNP) --batch-file $< --output-dir oncosnp/results --gcdir $(ONCOSNP_GC) --paramsfile $(ONCOSNP_PARAMS) --levelsfile $(ONCOSNP_LEVELS) --plot --fulloutput --subsample 2 --emiters 15 --stromal --intratumor --trainingstatesfile $(ONCOSNP_TRAINING_STATES) --tumourstatesfile $(ONCOSNP_TUMOUR_STATES) &> $(LOGDIR)/oncosnp.$*.log && touch $@
