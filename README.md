# Title
[![Build Status](https://travis-ci.org/cBioPortal/cbioportal.svg?branch=master)](https://travis-ci.org/jrflab/modules)


	## Introduction
	This is the implementation of the jrflab pipeline.

	## Installation
	The easiest way to download this pipeline is to clone the repository.

	```
	git clone https://github.com/jrflab/modules.git
	```

	## Dependencies
	An instance of [anaconda](https://www.anaconda.com) or [miniconda](https://conda.io/en/latest/miniconda.html)
	LSG or PBS for resource management  

		### Following R Packages
		- [xxx](https://)

	## Command line execution
	
		### Conventions
		- Sample names cannot have "/" or "." in them
		- Fastq files end in ".fastq.gz"
		- Fastq files are stored in DATA_DIR (Set as Environment Variable) 

		### Whole genome, whole exome and targeted sequencing
		- QC
		- BWA
- Broad Standard Practices on bwa bam  
- Haplotype Caller, Platupys, Bam2MPG, MuTect, Strelka  
- snpEff, Annovar, SIFT, pph2, Custom Annotation  
- Coverage Plot, Circos Plot, Hotspot Coverage Box Plot  
- Create input format for oncogenomics database (Patient Level)  
- Make Actionable Classification for Germline and Somatic Mutations   
- Copy number based on the simple T/N LogRatio (N cov >=30), Corrected for Total # Reads  
- Copy number, tumor purity using sequenza   
- LRR adjusted to center  
- Contamination using 
- HLA Typing
	* [xxx](http://)
	* [xxx](http://)  

### RNASeq:
- QC
- Tophat, STAR
- Broad Standard Practices on STAR bam
- fusion-catcher, tophat-fusion, deFuse
- Cufflinks (ENS and UCSC)
- Rsubread TPM (ENS, UCSC), Gene, Transcript and Exon Level
- In-house Exon Expression (ENS and UCSC)
- Haplotype Caller
- snpEff, Annovar, SIFT, pph2, Custom Annotation
- Actionable Fusion classification

### Patient:
- Genotyping On Patient. 
	1000g sites are evaluated for every library and then compared (all vs all)
	If two libraries come from a patient the match should be pretty good >80%
- Still to develop:
	If the match is below a certain threshold, break the pipeline for patient.

## Detailed usage
[wiki](https://github.com/jrflab/modules/wiki).

## Known issues
