export REF = b37

ALIGNER = bwamem
DUP_TYPE = markdup

TARGETS_FILE = test.bed
EXOME = true

# gatk options
HARD_FILTER_SNPS = true

QSUB_PRIORITY = -800

VALIDATION = false
