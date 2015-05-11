#!/usr/bin/env Rscript

library(knitr)
library(markdown)


args <- commandArgs(T)

if (length(args) < 2) stop("Need input and output script")

input <- args[1]
outPrefix <- args[2]
args <- args[c(-1,-2)]

#create output dirs
figPath <- file.path(outPrefix, 'figure/')
cachePath <- file.path(outPrefix, 'cache/')
dir.create(figPath, showWarnings = F, recursive = T)
dir.create(cachePath, showWarnings = F, recursive = T)

opts_chunk$set(dev = c("png", 'pdf'), dev.args = list(png = list(type = 'cairo-png')), cache.path = cachePath) # , fig.path = file.path('mutsig_report/figure/'))
opts_knit$set(root.dir = getwd(), base.dir = file.path(paste(outPrefix, '/', sep = '')), progress = F, verbose = T)

#options(warn = -1, error = quote({ traceback(2); q('no', status = 1) }))

knit(input, paste(outPrefix, '/index.md', sep = ''))
markdownToHTML(paste(outPrefix, '/index.md', sep = ''), paste(outPrefix, '/index.html', sep = ''))

