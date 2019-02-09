include modules/Makefile.inc

LOGDIR ?= log/defuse.$(NOW)
.PHONY: defuse

DEFUSE_SCRIPTS = /opt/common/CentOS_7/defuse/defuse-0.8.0/scripts
CONFIG = modules/config/defuse.conf
RELEASE = Grch37.p13
G38 = ${HOME}/share/reference/defuse/homo_sapiens/Ensembl/Grch38.p5
HSA = ${HOME}/share/reference/defuse/homo_sapiens/Ensembl/${RELEASE}
H37 = ${HSA}/Sequence/WholeGenomeFasta/Homo_sapiens.GRCh37.75.dna.primary_assembly.fa
WGFASTA_DIR = ${HSA}/Sequence/WholeGenomeFasta/
BWA = ${HSA}/Sequence/BWAIndex/genome.fa
BWT2 = ${HSA}/Sequence/Bowtie2Index/genome
BWT = ${HSA}/Sequence/BowtieIndex/Homo_sapiens.GRCh37.75.dna.primary_assembly
STARG = ${HSA}/Sequence/StarGenome/
STAR2PG = ${HSA}/Sequence/Star2pass/
RSEM = ${HSA}/Sequence/RSEMIndex/Grch37.p13
RTG = ${HSA}/Sequence/RTGIndex/SDF
CTAT19 = ${HSA}/Sequence/CTAT/GRCh37_gencode_v19/
CTAT19v2 = ${HSA}/Sequence/CTAT/GRCh37_gencode_v19_CTAT_lib_July192017/ctat_genome_lib_build_dir/
CTAT25 = ${HSA}/Sequence/CTAT/GRCh37_gencode_v25/
DEFUSE69 = ${HSA}/Sequence/defuse_e69/
DEFUSE75 = ${HSA}/Sequence/defuse_e75/
GMAP_DB = ${HSA}/Sequence/GmapDB
ERIC74 = ${HSA}/Sequence/ericscript_db/e74/
ERIC83 = $G38/Sequence/ericscript_db/
VARBIN_DIR = ${HSA}/Sequence/varbin/
VARBIN_BWT = ${VARBIN_DIR}/GRCh37.p13
ENSEMBL66 = ${HSA}/Genes/Homo_sapiens.GRCh37.66.gtf
GTF = ${HSA}/Annotation/Genes/gencode.v25lift37.annotation.gtf
REFFLAT = ${HSA}/Annotation/Genes/Homo.sapiens.annot.ucsc.txt
DBSNP = ${HSA}/Annotation/Variation/Homo_sapiens.vcf.gz
DBSNP_IC = ${HSA}/Annotation/Variation/Homo_sapiens_incl_consequences.vcf.gz
DBSNP_SOMATIC = ${HSA}/Annotation/Variation/Homo_sapiens_somatic.vcf.gz
DBSNP_SOMATIC_IC = ${HSA}/Annotation/Variation/Homo_sapiens_somatic_incl_consequences.vcf.gz
INDELVCF = ${HSA}/Annotation/Variation/
INDELMO = ${HSA}/Annotation/Variation/
SANGERINDEL = ${HSA}/Annotation/Variation/
SV_ALL = ${HSA}/Annotation/Variation/
SV_NE = ${HSA}/Annotation/Variation/
PERL = /usr/bin/perl

VPATH = fastq defuse defuse/tables

defuse : $(foreach sample,$(SAMPLES),defuse/$(sample).taskcomplete)

define defuse-single-sample
fastq/%.1.fastq fastq/%.2.fastq : fastq/%.1.fastq.gz fastq/%.2.fastq.gz
	$$(call RUN,-c -s 4G -m 9G,"gzip -d $$(<) $$(<<)")

defuse/tables/%.results.filtered.tsv : fastq/%.1.fastq fastq/%.2.fastq
	$$(call RUN,-c -n 10 -s 3G -m 4G -w 540,"$$(PERL) $$(DEFUSE_SCRIPTS)/defuse_run.pl -c $$(CONFIG) -d $$(DEFUSE75) -o defuse/$$* --res defuse/tables/$$*.results.tsv -resfil defuse/tables/$$*.results.filtered.tsv -1 fastq/$$*.1.fastq -2 fastq/$$*.2.fastq -p 10 -s direct")
	
defuse/%.taskcomplete : defuse/tables/%.results.filtered.tsv
	$$(call RUN,-c -s 1G -m 3G,"echo $$* > defuse/$$*.taskcomplete")
endef

$(foreach sample,$(SAMPLES),\
		$(eval $(call defuse-single-sample,$(sample))))
		
.DELETE_ON_ERROR:
.SECONDARY:
.PHONY: $(PHONY)
