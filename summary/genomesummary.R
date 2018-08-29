#!/usr/bin/env Rscript

file_names = dir(path="genome_stats/", pattern=".tsv", full.names=FALSE)
summary_scores = NULL
for (i in 1:length(file_names)) {
  data = read.csv(file=paste0("genome_stats/", file_names[i]), header=FALSE, sep="\t", stringsAsFactors=FALSE)
  summary_scores = cbind(summary_scores, data[,2])
}
summary_scores = cbind(data[,1], summary_scores)
colnames(summary_scores) = c("sample_names", gsub(".tsv", "", file_names))
write.table(summary_scores, file="genome_stats/genome_summary.tsv", col.names=TRUE, row.names=FALSE, sep="\t", quote=FALSE)

warnings()
