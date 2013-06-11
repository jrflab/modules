#SNV module
#Requirements:
#	1) *.filtered.rmdup.bam
#	2) Set the REF_GENOME variable to the genome fasta file that the libraries were aligned to
GET_PILEUP_RESULTS_BY_POS_SCRIPT=/share/lustre/gascoyne/scripts/getPileupResultsByPosns.pl
SNP6_POS_FILE=/share/lustre/gascoyne/references/affy6/affysnp6_positions.hg19.txt
KNOWNSNPS_FILE=/share/lustre/gascoyne/references/ucsc/common_snps.txt
SNP_EFF_DB=GRCh37.65
SNP_EFF_OPTS=-chr chr -no-downstream -no-intergenic -no-intron -no-upstream -no-utr -onlyCoding true
#CODON_FILE=/share/data/rmorin/data/codon_lookup.sort

THRESHOLD=0.5
MIN_READS=2
BASE_QUAL=15
MAPPING_QUAL=10

FILTER_SNVMIX_SCRIPT=${SCRIPTS_DIR}/filterSnvmix.pl
IDENTIFY_NONSYN_MUT_SCRIPT=${SCRIPTS_DIR}/identify_nonsynonymous_mutations.pl
ARTIFACTFINDER_SCRIPT=${SCRIPTS_DIR}/generateArtifactFinder.pl
SNVMIX_OPTS= -t MB -q ${BASE_QUAL} -Q ${MAPPING_QUAL}

#Training SNVMix - Grab only positions with SNP6 probes
snvmix/modelparams/%.snvmix2_modelparams : ${ALIGNER}/filtered_bam/%.filtered.rmdup.bam
	SGE_RREQ="-N $(@F) -l mem_free=1G -q all.q" \
	$(MKDIR) $(@D)/log;\
	$(SAMTOOLS) mpileup -s -f $(REF_GENOME) $< | perl ${GET_PILEUP_RESULTS_BY_POS_SCRIPT} -f SAMtools -p ${SNP6_POS_FILE} | $(SNVMIX) $(SNVMIX_OPTS) -T -m $@ 2> $(@D)/log/$*.log 

snvmix/results/%.snvmix2.txt : ${ALIGNER}/filtered_bam/%.filtered.rmdup.bam snvmix/modelparams/%.snvmix2_modelparams
	SGE_RREQ="-N $(@F) -l mem_free=1G -q all.q" \
	$(MKDIR) $(@D)/log;\
	$(SAMTOOLS) mpileup -s -f $(REF_GENOME) $< | $(SNVMIX) $(SNVMIX_OPTS) -p s -C -m $(word 2,$^) -o $@ 2> $(@D)/log/$*.log

snvmix/results/%.filtered.snvmix2.txt : snvmix/results/%.snvmix2.txt
	SGE_RREQ="-N $(@F) -l mem_free=1G -q all.q" \
	$(FILTER_SNVMIX_SCRIPT) $(THRESHOLD) $(MIN_READS) < $< > $@
		
snvmix/results/%.filtered.snvs.txt : snvmix/results/%.filtered.snvmix2.txt
	SGE_RREQ="-N $(@F) -l mem_free=1G -q all.q" \
	sort -k 1 $< | join -v 1 -a 1 -t\t - ${KNOWNSNPS_FILE} > $@ \
		&& rm -f $<

snvmix/snpeff/%.snpeff.tmp : snvmix/results/%.filtered.snvs.txt
	SGE_RREQ="-n $(@F) -l mem_free=4G -q all.q" \
	$(MKDIR) snvmix/snpeff/logs;\
	cat $< | perl -lane '@tmp=split(/:/,$$F[0]); print "$$tmp[0]\t$$tmp[1]\t$$F[1]\t$$F[2]\t+"' | ${SNP_EFF} -i txt -1 -c ${SNP_EFF_CONFIG} ${SNP_EFF_OPTS} ${SNP_EFF_DB} > $@ 2> $(@D)/logs/$(@F).log

#Building the SNPEff Header from the snpeff.tmp file 	
snvmix/snpeff/%.snpeff.header : snvmix/snpeff/%.snpeff.tmp
	SGE_RREQ="-n $(@F) -l mem_free=100M -q all.q" \
	grep "#" $< > $@;\
	tail -n 1 $@ | perl -F'\t' -lane 'print "$$F[0]\t$$F[1]\tRef_Count\tVariant_Count\t".join("\t",@F[2..$$#F])' > $(@D)/$*.tmp;\
	head -n 2 $@ | cat - $(@D)/$*.tmp > $(@D)/$*.tmp2;\
	mv $(@D)/$*.tmp2 $@ \
		&& rm $(@D)/$*.tmp

#Modify the snpeff output so that the chromosome and position are one column seperated by ":".  Th
snvmix/snpeff/%.snpeff.tmp2 : snvmix/snpeff/%.snpeff.tmp snvmix/snpeff/%.snpeff.header
	SGE_RREQ="-n $(@F) -l mem_free=100M -q all.q" \
	grep -v "#" $< | perl -F'\t' -lane 'print "$$F[0]:$$F[1]\t".join("\t",@F[2..$$#F])' > $@ \
		&& rm -f $<

snvmix/snpeff/%.snpeff.txt : snvmix/snpeff/%.snpeff.tmp2 snvmix/results/%.filtered.snvs.txt snvmix/snpeff/%.snpeff.header
	SGE_RREQ="-n $(@F) -l mem_free=100M -q all.q" \
	cat $(word 2,$^) | perl -lane '/.:(\d+),.:(\d+)/; print "$$F[0]\t$$1\t$$2"' | join -t$$'\t' - $< | perl -F'\t' -lane '@tmp=split(/:/,$$F[0]); print "$$tmp[0]\t$$tmp[1]\t".join("\t",@F[1..$$#F])' | cat $(word 3,$^) - > $@ \
		&& rm -f $< $(word 3, $^)
