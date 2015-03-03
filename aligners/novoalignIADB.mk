include modules/Makefile.inc

NOVOALIGN = $(HOME)/share/usr/bin/novoalignMPI
NOVOINDEX = $(HOME)/share/usr/bin/novoindex

REF_FASTA = $(HOME)/share/references/IADB_feb2011.fa
IADB_NIX = $(HOME)/share/references/IADB_feb2011.fa.nix

SAMPLE_FILE = samples.txt
SAMPLES = $(shell cat $(SAMPLE_FILE))

VPATH = bam

LOGDIR = iadb/log

.DELETE_ON_ERROR:
.SECONDARY: 
.PHONY : all

all : $(foreach sample,$(SAMPLES),iadb/processed_bam/$(sample).bam)

iadb/fastq/%.1.fastq.gz iadb/fastq/%.2.fastq.gz : %.bam
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,10G,12G)" $(MKDIR) $(@D) $(LOGDIR); $(BAM2FASTQ) --no-aligned -o $(@D)/$*#.fastq $< &> $(LOGDIR)/$(@F).log && mv $(@D)/$*_1.fastq $(@D)/$*.1.fastq && mv $(@D)/$*_2.fastq $(@D)/$*.2.fastq && gzip $(@D)/$*.1.fastq $(@D)/$*.2.fastq

iadb/bam/%.novoalign.bam : iadb/fastq/%.1.fastq.gz iadb/fastq/%.2.fastq.gz
	SGE_RREQ="$(SGE_RREQ) $(call MEM_FREE,14G,18G) -pe openmpi 5" $(MKDIR) $(@D) $(LOGDIR); mpiexec -np 5 $(NOVOALIGN) -i 200,50 -r A -R 0 -d $(IADB_NIX) -f $^ -o SAM $$'@RG\tID:$*\tPU:illumina\tLB:$*' 2> $(LOGDIR)/$(@F).log | $(SAMTOOLS) view -bS - > $@
	
iadb/processed_bam/%.bam : iadb/bam/%.novoalign.sorted.filtered.fixmate.markdup.bam
	$(MKDIR) $(@D); ln -v $< $@

include modules/bam_tools/processBam.mk
