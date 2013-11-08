NUM_ATTEMPTS ?= 3
NOW := $(shell date +"%F")
MAKELOG = log/$(@).$(NOW).log

QMAKE_BINARY = /common/sge/bin/lx-amd64/qmake
QMAKE = ~/share/scripts/qmake.pl -n $@.$(NOW) -r $(NUM_ATTEMPTS) -m -- $(QMAKE_BINARY)
MAKE = ~/share/scripts/qmake.pl -n $@.$(NOW) -r $(NUM_ATTEMPTS) -m -- make
QMAKEFLAGS = -cwd -v -inherit -q jrf.q
FLAGS = -j 25

.PHONY : all alignment variants qc fusions sum_reads

all : alignment variants qc fusions sum_reads

alignment :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/bowtieAlignerMD5.mk $(FLAGS)

variants : alignment
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/gatkVariantCaller.mk $(FLAGS)
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/mutect.mk $(FLAGS)
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/museqTN.mk $(FLAGS)

fusions :
	$(MAKE) -e -f ~/share/modules/defuse.mk -j15 -k
	$(MAKE) -e -f ~/share/modules/chimerascan.mk $(FLAGS)
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/tophatFusion.mk $(FLAGS)

sum_reads :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/sumRNASeqReads.mk $(FLAGS)

qc : alignment
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/bamIntervalMetrics.mk $(FLAGS)
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/fastqc.mk $(FLAGS)
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/bamMetrics.mk $(FLAGS)
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/rnaseqMetrics.mk $(FLAGS)
