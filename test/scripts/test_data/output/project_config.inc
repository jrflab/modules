ANN_FACETS = true
GATK_HARD_FILTER_SNPS = true
BAM_NO_REALN = false
VCF_POST_ANN_FILTER_EXPRESSION = ExAC_AF > 0.05
CLUSTER_ENGINE = SGE
BAM_NO_FILTER = false
BAM_DUP_TYPE = markdup
ALIGNER = bwamem
ANN_MUT_TASTE = false
EXOME = true
BAM_FIX_RG = false
FACETS_UNION_GATK_DBSNP = false
ANN_PATHOGEN = true
BAM_REPROCESS = false
GATK_POOL_SNP_RECAL = false
FACETS_MIN_NHET = 25
BAM_NO_SORT = false
BAM_PHRED64 = false
BAM_CHR1_BASE_RECAL = true
BAM_NO_RECAL = false
QSUB_PRIORITY = 0
INDEL_TYPE = strelka_varscan_indels
FACETS_CVAL1 = 150
USE_CLUSTER = true
FACETS_CVAL2 = 50
FACETS_PRE_CVAL = 50
REF = b37
SNV_TYPE = mutect

SAMPLES = GCT10T GCT11T GCT12T 
TUMOR_SAMPLES = GCT10T GCT11T GCT12T
NORMAL_SAMPLES = 
SAMPLE_PAIRS = 
SAMPLE_SETS = 

# tumor.normal = tumors

# normal.tumor = normal

# tumor.pair = tumor

# normal.pair = normal

# set.normal = set

# set.tumor = set

# set.pair = set

# sample_fastq_file
fq.GCT10T.0 = /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT10_T_IGO_06575_C_3/GCT10_T_IGO_06575_C_3_S7_L002_R1_001.fastq.gz /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT10_T_IGO_06575_C_3/GCT10_T_IGO_06575_C_3_S7_L002_R2_001.fastq.gz
fq.GCT10T.1 = /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT10_T_IGO_06575_C_3/GCT10_T_IGO_06575_C_3_S55_L005_R1_001.fastq.gz /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT10_T_IGO_06575_C_3/GCT10_T_IGO_06575_C_3_S55_L005_R2_001.fastq.gz
fq.GCT11T.0 = /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT11_T_IGO_06575_C_4/GCT11_T_IGO_06575_C_4_S54_L005_R1_001.fastq.gz /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT11_T_IGO_06575_C_4/GCT11_T_IGO_06575_C_4_S54_L005_R2_001.fastq.gz
fq.GCT11T.1 = /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT11_T_IGO_06575_C_4/GCT11_T_IGO_06575_C_4_S6_L002_R1_001.fastq.gz /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT11_T_IGO_06575_C_4/GCT11_T_IGO_06575_C_4_S6_L002_R2_001.fastq.gz
fq.GCT12T.0 = /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT12_T_IGO_06575_C_5/GCT12_T_IGO_06575_C_5_S56_L005_R1_001.fastq.gz /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT12_T_IGO_06575_C_5/GCT12_T_IGO_06575_C_5_S56_L005_R2_001.fastq.gz
fq.GCT12T.1 = /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT12_T_IGO_06575_C_5/GCT12_T_IGO_06575_C_5_S8_L002_R1_001.fastq.gz /home/limr/share/data/test/tseq/modules/test/scripts/test_data/rawdata/Project_06575_C/Sample_GCT12_T_IGO_06575_C_5/GCT12_T_IGO_06575_C_5_S8_L002_R2_001.fastq.gz
split.GCT10T.0 = GCT10T
split.GCT10T.1 = GCT10T
split.GCT10T = GCT10T.0 GCT10T.1

split.GCT11T.0 = GCT11T
split.GCT11T.1 = GCT11T
split.GCT11T = GCT11T.0 GCT11T.1

split.GCT12T.0 = GCT12T
split.GCT12T.1 = GCT12T
split.GCT12T = GCT12T.0 GCT12T.1

SPLIT_SAMPLES = GCT11T.0 GCT11T.1 GCT12T.1 GCT12T.0 GCT10T.1 GCT10T.0
