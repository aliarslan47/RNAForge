#!/usr/bin/env Rscript
# tximport: quant.sf'leri gen-seviyesine topla (lengthScaledTPM → m06 DESeq2 offset gerektirmez).
# argv: <tx2gene.tsv> <out_gene_counts.tsv> sid1 sf1 sid2 sf2 ...
suppressMessages(library(tximport))
args <- commandArgs(trailingOnly = TRUE)
tx2gene_path <- args[1]
out_path <- args[2]
rest <- args[-(1:2)]
sids <- rest[seq(1, length(rest), by = 2)]
sfs  <- rest[seq(2, length(rest), by = 2)]
names(sfs) <- sids
tx2gene <- read.table(tx2gene_path, header = FALSE, sep = "\t",
                      stringsAsFactors = FALSE)
txi <- tximport(sfs, type = "salmon", tx2gene = tx2gene,
                countsFromAbundance = "lengthScaledTPM")
counts <- round(txi$counts)                      # uzunluk-düzeltilmiş sayım
df <- data.frame(gene = rownames(counts), counts, check.names = FALSE)
write.table(df, out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("tximport ok:", nrow(counts), "genes x", ncol(counts), "samples\n")
