ifndef TITAN_INC
EXTRACT_ALLELE_READ_COUNTS = $(ANACONDA_PYTHON) $(HOME)/share/usr/TITANRunner-0.0.3/scripts/count.py
TITAN = $(RSCRIPT) modules/copy_number/runTitan.R
TITAN_SEG = $(PERL) $(HOME)/share/usr/TITANRunner-0.0.3/scripts/createTITANsegmentfiles.pl
SUMMARIZE_TITAN = $(RSCRIPT) modules/copy_number/summarizeTitan.R
ANNOTATE_TITAN_LOH_VCF = $(RSCRIPT) modules/copy_number/annotateTitanLOHVcf.R
NUM_CLUSTERS ?= $(shell seq 1 5)
PLOIDY_PRIORS = 2 3 4
DEFAULT_PLOIDY_PRIOR ?= 2

BQ_THRESHOLD ?= 20
MQ_THRESHOLD ?= 20
TITAN_WINDOW_SIZE ?= 10000

TITAN_SELF_TRANSITION ?= 1e15
TITAN_CLONAL_CLUSTER_TRANSITION ?= 5e5


override TITAN_OPTS := $(if $(UCSC_REF),--genomeStyle UCSC,--genomeStyle NCBI) $(if $(TARGETS_FILE),--targetBed $(TARGETS_FILE)) 
READ_COUNTER = $(HOME)/share/usr/bin/readCounter
MAP_COUNTER = $(HOME)/share/usr/bin/mapCounter
GC_COUNTER = $(HOME)/share/usr/bin/gcCounter

#VCF_FIELDS += titanCN titanMinorCN titanMajorCN titanCall titanMedianRatio titanMedianLogR
endif
TITAN_INC = true
