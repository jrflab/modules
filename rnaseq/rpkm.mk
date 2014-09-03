#Module calculates RPKM values.  Depends on the sumRNASeqReads.mk
RPKM_RSCRIPT = ${RSCRIPT} ~/gascoyne/scripts/calculateRPKM.R

rpkm/%.rpkm.txt : summarized_reads/%.summarized_reads.txt
	SGE_RREQ="-N $(@F) -l mem_free=1G -q all.q -now n" \
	$(MKDIR) $(@D)/logs;\
	$(RPKM_RSCRIPT) ${TXDB_FILE} $< $@ > $(@D)/logs/$*.log 2>&1
