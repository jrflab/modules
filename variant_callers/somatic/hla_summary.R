suppressPackageStartupMessages(library("optparse"))

options(warn = -1, error = quote({ traceback(); q('no', status = 1) }))

optList <- list(
                make_option("--sample_names", default = "NA", help = "tumor normal sample pair names")
                )

parser <- OptionParser(usage = "%prog [options]", option_list = optList)
arguments <- parse_args(parser, positional_arguments = T)
opt <- arguments$options

sample_names = unlist(strsplit(opt$sample_names, split=" ", fixed=TRUE))
hla_genotypes = list()
for (i in 1:length(sample_names)) {
	data = read.csv(file=paste0("hla_polysolver/", sample_names[i], "/winners.hla.txt"), header=FALSE, sep="\t", stringsAsFactors=FALSE)
	gen_1 = t(data[,2,drop=FALSE])
	gen_2 = t(data[,3,drop=FALSE])
	colnames(gen_1) = paste0(c("HLA-A", "HLA-B", "HLA-C"), "_1")
	colnames(gen_2) = paste0(c("HLA-A", "HLA-B", "HLA-C"), "_2")
	hla_genotypes[[i]] = cbind(gen_1, gen_2)
}
hla_genotypes = do.call(rbind, hla_genotypes)
hla_genotypes = cbind("SAMPLE_NAMES"=sample_names, hla_genotypes)
write.table(hla_genotypes, file="hla_polysolver/summary/genotype_summary.txt", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)
