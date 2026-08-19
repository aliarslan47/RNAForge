# rnaforge/scripts/deseq2.R
# Argümanlar: counts.tsv coldata.tsv design reference out_dir
# DESeq2 ile DE; temiz TSV'ler yazar. Hata → stderr + nonzero exit (sessiz kısmi çıktı yok).
args <- commandArgs(trailingOnly = TRUE)
counts_tsv <- args[1]; coldata_tsv <- args[2]; design_str <- args[3]
reference  <- args[4]; out_dir    <- args[5]
# 6. arg (opsiyonel): açık kontrastlar "test:ref;test2:ref2". Boş → varsayılan (son-vs-ilk).
contrasts_spec <- if (length(args) >= 6) args[6] else ""
suppressMessages(library(DESeq2))

counts  <- read.delim(counts_tsv, row.names = 1, check.names = FALSE)
coldata <- read.delim(coldata_tsv, row.names = 1, check.names = FALSE)
counts  <- counts[, rownames(coldata), drop = FALSE]        # örnek sırasını hizala
coldata$condition <- factor(coldata$condition)
if (nzchar(reference)) coldata$condition <- relevel(coldata$condition, ref = reference)
if ("batch" %in% colnames(coldata)) coldata$batch <- factor(coldata$batch)
if ("subject" %in% colnames(coldata)) coldata$subject <- factor(coldata$subject)

dds <- DESeqDataSetFromMatrix(countData = as.matrix(counts),
                              colData = coldata, design = as.formula(design_str))
# Dispersiyon uyumu: parametrik başarısız olursa (ör. az gen / zor veri) sessizce
# çökme yerine local'e, o da olmazsa mean'e düş. Hangi uyumun kullanıldığı stderr'e yazılır.
dds <- tryCatch(
  DESeq(dds, quiet = TRUE),
  error = function(e1) {
    message("parametric dispersion fit failed (", conditionMessage(e1),
            "); retrying with fitType='local'")
    tryCatch(
      DESeq(dds, quiet = TRUE, fitType = "local"),
      error = function(e2) {
        message("local dispersion fit failed (", conditionMessage(e2),
                "); retrying with fitType='mean'")
        tryCatch(
          DESeq(dds, quiet = TRUE, fitType = "mean"),
          error = function(e3) {
            # Dejenere durum (ör. neredeyse tekdüze sayımlar): dispersiyon trendi
            # uydurulamaz. DESeq2'nin önerdiği yol — gen-bazlı tahminleri doğrudan kullan.
            message("all dispersion trend fits failed (", conditionMessage(e3),
                    "); using gene-wise dispersion estimates directly")
            d <- estimateSizeFactors(dds)
            d <- estimateDispersionsGeneEst(d)
            dispersions(d) <- mcols(d)$dispGeneEst
            nbinomWaldTest(d)
          }
        )
      }
    )
  }
)
write_results <- function(res, path) {
  res <- as.data.frame(res)
  res <- cbind(gene = rownames(res), res)
  write.table(res, path, sep = "\t", quote = FALSE, row.names = FALSE)
}

cond_levels <- levels(coldata$condition)
# Varsayılan kontrast etiketi (son-vs-ilk); açık kontrast verilirse ilk çift ezer.
primary_contrast <- paste(rev(cond_levels)[1], "vs", cond_levels[1])

if (nzchar(contrasts_spec)) {
  # Açık kontrastlar: her biri için ayrı dosya; İLK çift birincil deseq2_results.tsv olur
  # (downstream/rapor birincili tüketir). Geriye uyumlu: spec boşsa bu blok atlanır.
  pairs <- strsplit(contrasts_spec, ";", fixed = TRUE)[[1]]
  first <- TRUE
  for (pair in pairs) {
    tr <- strsplit(pair, ":", fixed = TRUE)[[1]]
    test <- tr[1]; ref <- tr[2]
    if (!(test %in% cond_levels) || !(ref %in% cond_levels)) {
      stop(sprintf("contrast '%s vs %s': condition level not found (levels: %s)",
                   test, ref, paste(cond_levels, collapse = ", ")))
    }
    res_i <- results(dds, contrast = c("condition", test, ref))
    write_results(res_i, file.path(out_dir, sprintf("deseq2_results.%s_vs_%s.tsv", test, ref)))
    if (first) {
      write_results(res_i, file.path(out_dir, "deseq2_results.tsv"))
      primary_contrast <- paste(test, "vs", ref)
      first <- FALSE
    }
  }
} else {
  write_results(results(dds), file.path(out_dir, "deseq2_results.tsv"))
}

norm <- counts(dds, normalized = TRUE)
normdf <- cbind(gene = rownames(norm), as.data.frame(norm))
write.table(normdf, file.path(out_dir, "normalized_counts.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# Dispersiyon tahminleri (m07 dispersiyon figürü için). dispFit fallback yolunda olmayabilir -> NA.
mc <- mcols(dds)
dispfit <- if ("dispFit" %in% colnames(mc)) mc$dispFit else rep(NA_real_, nrow(dds))
disp <- data.frame(gene_id = rownames(dds), baseMean = mc$baseMean,
                   dispGeneEst = mc$dispGeneEst, dispFit = dispfit,
                   dispFinal = dispersions(dds))
write.table(disp, file.path(out_dir, "dispersions.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)

# min koşul-içi replika korelasyonu (log2 normalize sayım, Pearson)
ln <- log2(norm + 1)
mincor <- 1.0
for (lvl in levels(coldata$condition)) {
  cols <- rownames(coldata)[coldata$condition == lvl]
  if (length(cols) >= 2) {
    cm <- cor(ln[, cols])
    mincor <- min(mincor, min(cm[upper.tri(cm)]))
  }
}
writeLines(c(paste0("min_replicate_correlation\t", mincor),
             paste0("contrast\t", primary_contrast),
             paste0("n_genes\t", nrow(dds))),
           file.path(out_dir, "de_metrics.tsv"))
