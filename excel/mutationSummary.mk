include modules/variant_callers/somatic/somaticVariantCaller.inc
include modules/Makefile.inc

ALLTABLES_HIGH_MODERATE_MUTECT = alltables/allTN.$(call VCF_SUFFIXES,mutect).tab.high_moderate.novel.txt
ALLTABLES_LOW_MODIFIER_MUTECT = alltables/allTN.$(call VCF_SUFFIXES,mutect).tab.low_modifier.novel.txt
ALLTABLES_HIGH_MODERATE_STRELKA_VARSCAN = alltables/allTN.strelka_varscan_indels.tab.high_moderate.novel.txt
ALLTABLES_LOW_MODIFIER_STRELKA_VARSCAN = alltables/allTN.strelka_varscan_indels.tab.low_modifier.novel.txt

mutation_summary: excel/mutation_summary.xlsx


excel/mutation_summary.xlsx : $(ALLTABLES_HIGH_MODERATE_MUTECT) $(ALLTABLES_LOW_MODIFIER_MUTECT) $(ALLTABLES_HIGH_MODERATE_STRELKA_VARSCAN) $(ALLTABLES_LOW_MODIFIER_STRELKA_VARSCAN)
	$(INIT) source $(ANACONDA_27_ENV)/bin/activate $(ANACONDA_27_ENV); \
	python modules/scripts/tsvToExcel.py $< $@ 'mutect_high_moderate' && \
	python modules/scripts/tsvToExcel.py $(<<) $@ 'mutect_low_modifier' && \
	python modules/scripts/tsvToExcel.py $(<<<) $@ 'strelka_varscan_high_moderate' && \
	python modules/scripts/tsvToExcel.py $(<<<<) $@ 'strelka_varscan_low_modifier' && \
	python modules/scripts/tsvToExcel.py -c CHROM,POS,TUMOR_SAMPLE,ANN[*].GENE,ANN[*].HGVS_P,ANN[*].EFFECT,TUMOR.AD,NORMAL.AD,TUMOR.DP,NORMAL.DP,chasm_score,dbNSFP_MutationTaster_pred,FATHMM_pred,cancer_gene_census,kandoth,lawrence,hap_insuf,REF,ALT,ANN[*].IMPACT $(<) $@ 'SNV_NONSILENT' && \
	python modules/scripts/tsvToExcel.py -c CHROM,POS,TUMOR_SAMPLE,ANN[*].GENE,ANN[*].HGVS_P,ANN[*].EFFECT,TUMOR.AD,NORMAL.AD,TUMOR.DP,NORMAL.DP,chasm_score,dbNSFP_MutationTaster_pred,FATHMM_pred,cancer_gene_census,kandoth,lawrence,hap_insuf,REF,ALT,ANN[*].IMPACT $(<<) $@ 'SNV_SILENT' && \
	python modules/scripts/tsvToExcel.py -c CHROM,POS,TUMOR_SAMPLE,ANN[*].GENE,ANN[*].HGVS_P,ANN[*].EFFECT,TUMOR.AD,NORMAL.AD,TUMOR.DP,NORMAL.DP,chasm_score,dbNSFP_MutationTaster_pred,FATHMM_pred,cancer_gene_census,kandoth,lawrence,hap_insuf,REF,ALT,ANN[*].IMPACT $(<<<) $@ 'INDEL_NONSILENT' && \
	python modules/scripts/tsvToExcel.py -c CHROM,POS,TUMOR_SAMPLE,ANN[*].GENE,ANN[*].HGVS_P,ANN[*].EFFECT,TUMOR.AD,NORMAL.AD,TUMOR.DP,NORMAL.DP,chasm_score,dbNSFP_MutationTaster_pred,FATHMM_pred,cancer_gene_census,kandoth,lawrence,hap_insuf,REF,ALT,ANN[*].IMPACT $(<<<<) $@ 'INDEL_SILENT'
