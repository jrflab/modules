# annotate titan module
include modules/Makefile.inc

LOGDIR = log/mutsig_report.$(NOW)

MUTSIG_REPORT = modules/mut_sigs/mutSigReport.Rmd
KNIT = $(RSCRIPT) modules/scripts/knit.R
ALEXANDROV_DATA = $(HOME)/share/reference/Alexandrov_NMF_signatures.txt

#$(eval mutsig-report-name-vcfs,name,vcfs)
# creates [name]_mutsig_report target
define mutsig-report-name-vcfs
PHONY += $1_mutsig_report 
$1_mutsig_report : mutsig_report/$1/index.html
mutsig_report/$1/index.html : $2
	$$(call LSCRIPT_NAMED_MEM,$1_mutsig_report,8G,20G,"$$(KNIT) $$(MUTSIG_REPORT) $(@D) --outDir $(@D) --name $1 --alexandrovData $$(ALEXANDROV_DATA) $$(if $$(TARGETS_FILE),--targetBed $$(TARGETS_FILE)) $$^")
endef


