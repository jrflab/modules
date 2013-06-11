# Targeted Sequencing QC module
# needs bed file with targets and baits
#
#
# Author: Raymond Lim <raylim@mm.st>
#
#
# QC
# 
TARGETS_BED := data/baits.bed
BAITS_BED := data/baits.bed

TEQC_REPORT := /share/lustre/gascoyne/scripts/TEQCreport.R --ref=$(REF)

TEQCreports : $(foreach sample,$(SAMPLES),TEQCreports/$(sample))

TEQCreports/% : %.bam $(TARGETS_BED) $(BAITS_BED)
	SGE_RREQ="-q all.q -l mem_free=4G -now no" $(MKDIR) $@; $(TEQC_REPORT) --outputDir=$@ --sampleName=$* --targetsName=ubr5_et_al --saveWorkspace --duplicatesPlot --baitsFile=$(BAITS_BED) --covGCPlot --covUniformityPlot $< $(BAITS_BED) &> log/$*.teqc_report.log || rm -rf $@
