NUM_ATTEMPTS ?= 3
NOW := $(shell date +"%F")
MAKELOG = log/$(@).$(NOW).log

QMAKE_BINARY = /common/sge/bin/lx-amd64/qmake
QMAKE = ~/share/scripts/qmake.pl -n $@.$(NOW) -r $(NUM_ATTEMPTS) -m -- $(QMAKE_BINARY)
MAKE = ~/share/scripts/qmake.pl -n $@.$(NOW) -r $(NUM_ATTEMPTS) -m -- make
QMAKEFLAGS = -cwd -v -inherit -q jrf.q
FLAGS = -j 25

.PHONY : all

all : variants qc cnv qc

alignment.timestamp :
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/bwaAlignerMD5.mk $(FLAGS) && touch $@

variants: alignment.timestamp
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/gatkVariantCaller.mk $(FLAGS)
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/mutect.mk $(FLAGS)
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/museqTN.mk $(FLAGS)
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/pindel.mk $(FLAGS)

cnv : alignment.timestamp
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/controlFreeCTN.mk $(FLAGS)
	$(MAKE) $(MAKEFLAGS) -f ~/share/modules/varscanTN.mk $(FLAGS) cnv

qc : alignment.timestamp
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/bamIntervalMetrics.mk $(FLAGS)
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/fastqc.mk $(FLAGS)
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/bamMetrics.mk $(FLAGS)
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/rnaseqMetrics.mk $(FLAGS)
	$(MAKE) $(MAKEFLAGS) -e -f ~/share/modules/qualimap.mk $(FLAGS) $(TARGET)
