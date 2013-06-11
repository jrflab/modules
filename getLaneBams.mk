# transfer lane bams for a set of samples
# requires a lane sample file
# 	format is tab-delimited "sample lane_id bam_file_path"
# 	Use scripts/findLanes.sh to generate from a samples.txt file
include ~/share/modules/Makefile.inc

LANE_SAMPLE_FILE ?= samples.lane.txt
SAMPLES = $(shell cut -f1 $(LANE_SAMPLE_FILE) | sort | uniq)
LANES = $(shell cut -f2 $(LANE_SAMPLE_FILE))
LANE_BAMS = $(shell cut -f3 $(LANE_SAMPLE_FILE))

$(foreach i,$(shell seq 1 $(words $(LANES))),$(eval lane_lookup.$(word $i,$(LANES)) := $(word $i,$(LANE_BAMS))))

SAMPLE_DIRS = $(HOME)/share/references/sample_dirs.txt

RSYNC = rsync -avr --progress --partial

TRANSFER_HOST = xhost09
APOLLO_SSH = ssh apollo LOADEDMODULES="default-manpath/1.0.1:mpiJava/1.2.5:lam-oscar-7.1.2:mpi/lam-7.1.2:switcher/1.0.13:pvm/3.4.5+6:sge/6.2u5:oscar-modules/1.0.5" COMMD_HOST=apollo SGE_CELL=apollo SGE_ROOT=/opt/sge SGE_CLUSTER_NAME=apollo 
QSUB = qsub -q thosts.q -now no -P transfer -sync y -N $(@F) -e $(TRANSFER_LOGDIR)/$(@F).err.log -o $(TRANSFER_LOGDIR)/$(@F).out.log
#TRANSFER = ssh $(TRANSFER_HOST) ${RSYNC}
CP = cp -v -p

TRANSFER_LOGDIR = ~/transfer_log
LOGDIR = log/get_lane_bams.$(NOW)

.PHONY : all

all : $(foreach lane,$(LANES),gsc_bam/$(lane).bam)

gsc_bam/%.bam : 
	$(INIT) \
	DEST=$(CURDIR)/$@; if [[ "$$DEST" != /genesis/* ]]; then DEST=/genesis/$$DEST; fi; \
	ssh $(TRANSFER_HOST) $(MKDIR) $(TRANSFER_LOGDIR); \
	echo "$(CP) $(lane_lookup.$*) $$DEST" | $(APOLLO_SSH) $(QSUB)

