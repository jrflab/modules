# find and transfer bams from the GSC projects directory
# requires SSH setup for no password

include ~/share/modules/Makefile.inc

SAMPLE_FILE = samples.txt
SAMPLES = $(shell cat $(SAMPLE_FILE))

SAMPLE_DIRS = $(HOME)/share/references/sample_dirs.txt

RSYNC = rsync -avr --progress --partial

TRANSFER_HOST = xhost09
FIND_BAM = sh $(HOME)/share/scripts/findBam.sh $(SAMPLE_DIRS)
SSH = ssh apollo LOADEDMODULES="default-manpath/1.0.1:mpiJava/1.2.5:lam-oscar-7.1.2:mpi/lam-7.1.2:switcher/1.0.13:pvm/3.4.5+6:sge/6.2u5:oscar-modules/1.0.5" COMMD_HOST=apollo SGE_CELL=apollo SGE_ROOT=/opt/sge SGE_CLUSTER_NAME=apollo 
QSUB = qsub -q thosts.q -now no -P transfer -sync y -N $(@F) -e $(TRANSFER_LOGDIR)/$(@F).err.log -o $(TRANSFER_LOGDIR)/$(@F).out.log
#TRANSFER = ssh $(TRANSFER_HOST) ${RSYNC}
CP = cp -v -p

TRANSFER_LOGDIR = ~/transfer_log

all : $(foreach sample,$(SAMPLES),gsc_bam/$(sample).bam)
transfer : transfer.sh

gsc_bam/%.bam :
	$(MKDIR) $(@D); \
	DEST=$(CURDIR)/$@; if [[ "$$DEST" != /genesis/* ]]; then DEST=/genesis/$$DEST; fi; \
	BAM=`ssh ${TRANSFER_HOST} $(FIND_BAM) $*`; \
	$(SSH) $(MKDIR) $(TRANSFER_LOGDIR); \
	echo "$(CP) $$BAM $$DEST" | $(SSH) $(QSUB)

# Builds a transfer script that you can then use on the apollo cluster
transfer.sh : 
	for sample in ${SAMPLES}; do \
		DEST=$(CURDIR)/gsc_bam/$${sample}.bam; if [[ "$$DEST" != /genesis/* ]]; then DEST=/genesis/$$DEST; fi; \
		BAM=`${FIND_BAM} $${sample}`; \
		echo "${RSYNC} $$BAM $$DEST" >> $@; \
	done;

# BAM=`$(FIND_BAM) $*`; ssh apollo "SGE_ROOT=/opt/sge SGE_CELL=apollo qrsh -q thosts.q -P transfer -now no /bin/cp $$BAM /genesis/$(CURDIR)/$@"

$(SAMPLE_DIRS) :
	find -L /projects/analysis /projects/seq_analysis /archive/analysis* /archive/solexa* -maxdepth 2 -name 'HS*' -or -name 'A*' > $@;
