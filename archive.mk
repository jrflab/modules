# archive bams and various tables to archive dir

ARCHIVE_DIR = /mount/limr/zedshared
RELATIVE_PATH = $(shell echo "$(CURDIR)" | sed 's/.*\(projects\|data\)\//\1\//')
ARCHIVE_PATH = $(ARCHIVE_DIR)/$(RELATIVE_PATH)

LOGDIR = log/archive.$(NOW)

RSYNC_OPTS = --verbose --stats --recursive -a -0 --prune-empty-dirs
RSYNC = rsync

.SECONDARY:
.DELETE_ON_ERROR:
.PHONY : archive_bams archive_tables

archive_all : archive_bams archive_tables

archive_bams :
	$(RSYNC) $(RSYNC_OPTS) bam $(ARCHIVE_PATH)

archive_tables :
	$(RSYNC) $(RSYNC_OPTS) alltables tables vcf metrics $(ARCHIVE_PATH)
