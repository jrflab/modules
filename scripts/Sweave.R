#!/usr/bin/env Rscript
# runs sweave a file, passes on arguments

suppressPackageStartupMessages(library("optparse"));

allArguments <- commandArgs(trailingOnly = T)

if (length(allArguments) < 1) {
    cat("Need Rnw file");
    stop();
}

rnwFile <- allArguments[1];
arguments <<- allArguments[-1];

Sweave(rnwFile);
