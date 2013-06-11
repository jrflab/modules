# run miso events analysis on a set of bams
include ~/share/modules/Makefile.inc

SAMPLE_FILE ?= samples.txt
SAMPLES ?= $(shell cat $(SAMPLE_FILE))

LOGDIR = log/miso.$(NOW)

VPATH = bam

MISO = $(PYTHON) $(HOME)/share/usr/bin/run_miso.py
MISO_INDEX_DIR = $(HOME)/share/references/miso/hg19
#EVENTS_FILE = $(MISO_INDEX_DIR)/exon_events_list.txt
EVENT_TYPES = A3SS A5SS AFE ALE MXE SE TandemUTR RI
#EVENT_TYPES = RI
#GENES_FILE = $(MISO_INDEX_DIR)/gene_list.txt
#GENES = $(shell cat $(GENES_FILE))

PE_UTILS = $(PYTHON) $(HOME)/share/usr/bin/pe_utils.py
CONST_EXONS_GFF = $(HOME)/share/references/miso/exons/Homo_sapiens.GRCh37.65.with_chr.min_800.const_exons.gff

NUM_CHUNKS = 100
CHUNK_SEQ = $(shell seq 0 $(shell expr $(NUM_CHUNKS) - 1))

.SECONDARY : 
.DELETE_ON_ERROR :
.PHONY : metrics events_analyses summaries process_bams all chunk

all : process_bams metrics summaries chunk tar

# define CHUNK and SAMPLE to process a chunk
ifdef CHUNK_FILE
CHUNK_EVENTS = $(shell cat $(CHUNK_FILE))
chunk : $(addprefix miso/$(SAMPLE).events_analysis/,$(addsuffix .miso,$(CHUNK_EVENTS)))
endif
	
summaries : $(foreach event_type,$(EVENT_TYPES),$(foreach sample,$(SAMPLES),miso/summary_output/$(sample)/summary/$(event_type).miso_summary))

events_analyses : $(foreach sample,$(SAMPLES),miso/$(sample).events_analysis_timestamp)

process_bams : $(foreach sample,$(SAMPLES),bam/$(sample).bam) $(foreach sample,$(SAMPLES),bam/$(sample).bam.bai)


metrics : $(foreach sample,$(SAMPLES),miso/metrics/$(sample).read_len) $(foreach sample,$(SAMPLES),miso/metrics/$(sample).bam.insert_len) $(foreach chunk,$(CHUNK_SEQ),miso/chunks/$(chunk).txt)

tar : $(foreach sample,$(SAMPLES),miso/$(sample).events_analysis.tar.bz2)

miso/metrics/%.read_len : %.bam
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,4G,5G)" $(MKDIR) $(@D); $(SAMTOOLS) view $< | awk '{ print length($$10) }' | sort -n | uniq -c | sort -rn | sed 's/^ \+//' > $@

miso/metrics/%.bam.insert_len : %.bam
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,3G,4G)" $(MKDIR) $(@D) $(LOGDIR); $(PE_UTILS) --compute-insert-len $< $(CONST_EXONS_GFF) --output-dir $(@D) &> $(LOGDIR)/$(@F).log

miso/chunks/%.txt : miso/exon_events_list.txt
	$(MKDIR) $(@D); awk 'NR % $(NUM_CHUNKS) == $*' $< > $@

miso/exon_events_list.txt :
	( cd $(MISO_INDEX_DIR) && find $(EVENT_TYPES) -name '*.pickle' ) | sed 's/\.pickle$$//; s/:/\\:/g; s/|/\\|/g;' > $@

# $(call events-analysis,sample)
define events-analysis
miso/$1.events_analysis/.%.timestamp : miso/chunks/%.txt
	SGE_RREQ="-q all.q -N $1_$$*.miso_events_analysis -now no $$(call MEM_FREE,2G,3G) -pe $$(PARALLEL_ENV) 5" \
	 echo "$$$$HOSTNAME" >>  $$(LOGDIR)/$1.$$*.make_events_analysis.log; \
	$$(MKDIR) $$(@D) $$(LOGDIR); make -e -j 5 -f ~/share/modules/miso.mk SAMPLE=$1 CHUNK_FILE=$$< chunk >> $$(LOGDIR)/$1.$$*.make_events_analysis.log 2>&1 && touch $$@

miso/$1.events_analysis_timestamp : $$(foreach chunk,$$(CHUNK_SEQ),miso/$1.events_analysis/.$$(chunk).timestamp)
	$$(MKDIR) $$(@D); touch $$@
endef
$(foreach sample,$(SAMPLES),$(eval $(call events-analysis,$(sample))))

miso/%.events_analysis.tar.bz2 : miso/%.events_analysis_timestamp
	tar -C $(@D) -cjf $@ $*.events_analysis  # && $(RM) -r miso/$*.events_analysis

# $(call compute-gene-psi,sample,event_type)
define compute-gene-psi
miso/$1.events_analysis/$2/%.miso : $1.bam miso/metrics/$1.bam.insert_len miso/metrics/$1.read_len
	SGE_RREQ="$$(SGE_RREQ) $$(call MEM_FREE,1G,2G) -l tmpfree=10G" \
	 INS=`head -1 $$(word 2,$$^) | sed 's/.*mean=\([0-9]\+\)\.\+[0-9]\+,sdev=\([0-9]\+\).*/\1 \2/'`; \
	 READ_LEN=`head -1 $$(word 3,$$^) | cut -f 2 -d' '`; \
	 $$(MKDIR) $$(@D); \
	 $$(RM) "$$@" && \
	 $$(MISO) --read-len $$$$READ_LEN --paired-end $$$$INS \
	 --compute-gene-psi "$$(subst \,,$$(notdir $$*))" "$$(MISO_INDEX_DIR)/$2/$$(subst \,,$$*).pickle" \
	 $$< miso/$1.events_analysis/$2 &> /dev/null
endef
$(foreach sample,$(SAMPLES),$(foreach event_type,$(EVENT_TYPES),$(eval $(call compute-gene-psi,$(sample),$(event_type)))))

# $(call miso-summary,event-type)
define miso-summary
miso/summary_output/%/summary/$1.miso_summary : miso/%.events_analysis_timestamp
	SGE_RREQ="-q all.q -now no -N $$*_$1.miso_summary $$(call MEM_FREE,10G,14G)" $$(MKDIR) $$(@D); $$(MISO) --summarize-samples miso/$$*.events_analysis/$1 miso/summary_output/$$* &> $$(LOGDIR)/$$(@F).log
endef
$(foreach event_type,$(EVENT_TYPES),$(eval $(call miso-summary,$(event_type))))

include ~/share/modules/processBam.mk
