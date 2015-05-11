#!/usr/bin/env Rscript
# Run absolute

suppressPackageStartupMessages(library("optparse"));
suppressPackageStartupMessages(library(HAPSEG));

options(warn = -1, error = quote({ traceback(2); q('no', status = 1) }))

optList <- list(
        make_option("--disease", default = 'breastcancer', help = "disease [default %default]"),
        make_option("--platform", default = "SNP_6.0", help = "platform [default %default]"),
        make_option("--plate", default = 'Plate', help = "Plate name [default %default]"),
        make_option("--callsFile", default = NULL, help = "birdseed calls file"),
        make_option("--summaryFile", default = NULL, help = "APT summary file"),
        make_option("--clustersFile", default = NULL, help = "birdseed models file"),
        make_option("--ref", default = 'hg19', help = "genome build [default %default]"),
        make_option("--resultsDir", default = NULL, help = "results directory"),
        make_option("--phasedBGLDir", default = NULL, help = "Phased BGL directory"),
        make_option("--outFile", default = NULL, help = "output file")
        )
parser <- OptionParser(usage = "%prog tumour [normal]", option_list = optList);
arguments <- parse_args(parser, positional_arguments = T);
opt <- arguments$options;

if (is.null(opt$outFile)) {
    cat("Need output file\n");
    print_help(parser);
    stop();
} else  if (is.null(opt$resultsDir)) {
    cat("Need results dir\n");
    print_help(parser);
    stop();
} else if (is.null(opt$summaryFile)) {
    cat("Need summary file\n");
    print_help(parser);
    stop();
} else if (is.null(opt$phasedBGLDir)) {
    cat("Need phased BGL directory\n");
    print_help(parser);
    stop();
} else if (is.null(opt$callsFile)) {
    cat("Need calls file\n");
    print_help(parser);
    stop();
} else if (is.null(opt$clustersFile)) {
    cat("Need clusters file\n");
    print_help(parser);
    stop();
} else if (length(arguments$args) < 1) {
    cat("Need tumour sample\n");
    print_help(parser);
    stop();
}

tumour <- arguments$args[1]
if (length(arguments$args) > 1) {
    normal <- arguments$args[2]
    useNormal <- T
} else {
    useNormal <- F
    normal <- NULL
}


RunHapSeg(out.file = opt$outFile, plate.name = opt$plate, array.name = tumour, seg.fn = NULL, 
	snp.fn = opt$summaryFile, 
	genome.build = opt$ref, results.dir = opt$resultsDir,
	platform = opt$platform, use.pop = "CEPH", impute.gt = TRUE, plot.segfit = TRUE, 
	merge.small = TRUE, merge.close = TRUE, min.seg.size = 3, normal = FALSE, out.p = 0.05, 
	seg.merge.thresh = 0.0001, phased.bgl.dir = opt$phasedBGLDir,
	drop.x = FALSE, drop.y = TRUE, calls.fn = opt$callsFile , 
	calibrate.data = TRUE, use.normal = useNormal,
    mn.sample = normal, 
	clusters.fn = opt$clustersFile,
	snp.file.parser = AptSnpFileParser, clusters.file.parser = BirdseedClustersFileParser, 
	verbose = TRUE, adj.atten = 0)

