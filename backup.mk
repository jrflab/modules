# backup data/results to lustre

include ~/share/modules/Makefile.inc

SCRATCH_DATA = /genesis/scratch/fcchan/data
SCRATCH_PROJECTS = /genesis/scratch/fcchan/projects
PROJECTS = /projects/fcchan/projects

RELPATH = $(shell echo $(PWD) | sed 's;.*share/;;')
DATA_DEST = beast:/share/lustre/clc
RESULTS_DEST = beast:/share/lustre/projects/clc

RESULT_FILETYPES = txt vcf vcf.gz gz bz2
DATA_FILETYPES = bam bai

RSYNC = ssh thost05 rsync --verbose --progress -H -a $(foreach filetype,$(1),--include='*.$(filetype)') --include='*/' --exclude='log*' --exclude='*' 
RSYNC_DATA = $(call RSYNC,$(DATA_FILETYPES))
RSYNC_RESULTS = $(call RSYNC,$(RESULT_FILETYPES))
.PHONY: all data results current_dir

all : data results


# backup the current directory
current_dir :
	RELPATH=`echo $$PWD | sed 's;.*share/;;'`;
	$(RSYNC_DATA) $$PWD $(DATA_DEST)/$$RELPATH/*; \
	$(RSYNC_DATA) $$PWD $(RESULTS_DEST)/$$RELPATH

# backup data (bam files) to /share/lustre/clc
data :
	$(RSYNC_DATA) $(SCRATCH_DATA) $(DATA_DEST); \
	$(RSYNC_DATA) $(SCRATCH_PROJECTS) $(DATA_DEST); \
	$(RSYNC_DATA) $(PROJECTS) $(DATA_DEST)

# backup results to /share/lustre/projects/clc
results :
	$(RSYNC_RESULTS) $(SCRATCH_DATA) $(RESULTS_DEST); \
	$(RSYNC_RESULTS) $(SCRATCH_PROJECTS) $(RESULTS_DEST); \
	$(RSYNC_RESULTS) $(PROJECTS) $(RESULTS_DEST)
