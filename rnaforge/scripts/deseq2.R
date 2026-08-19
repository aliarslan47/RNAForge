# rnaforge/scripts/deseq2.R
# Argümanlar: counts.tsv coldata.tsv design reference out_dir
# DESeq2 ile DE; temiz TSV'ler yazar. Hata → stderr + nonzero exit (sessiz kısmi çıktı yok).
args <- commandArgs(trailingOnly = TRUE)
counts_tsv <- args[1]; coldata_tsv <- args[2]; design_str <- args[3]
reference  <- args[4]; out_dir    <- args[5]
# 6. arg (opsiyonel): açık kontrastlar "factor:test:ref;factor2:test2:ref2". factor=condition
# (ANA tetkik faktörü) çoğu koşuda; condition-dışı faktör de tetkik edilebilir (Faz 3).
# Boş → varsayılan (condition son-vs-ilk).
contrasts_spec <- if (length(args) >= 6) args[6] else ""
suppressMessages(library(DESeq2))

counts  <- read.delim(counts_tsv, row.names = 1, check.names = FALSE)
coldata <- read.delim(coldata_tsv, row.names = 1, check.names = FALSE)
counts  <- counts[, rownames(coldata), drop = FALSE]        # örnek sırasını hizala
# coldata'daki HER sütunu (condition/batch/subject + keyfi kovaryatlar) faktörle —
# özel-durum yerine generic döngü (Faz 3: keyfi faktör adları). condition ANA tetkik
# faktörü olduğundan reference'a relevel edilir; diğer faktörlerin test/ref'i kontrastta açık.
for (col in colnames(coldata)) coldata[[col]] <- factor(coldata[[col]])
if (nzchar(reference)) coldata$condition <- relevel(coldata$condition, ref = reference)

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
  # Açık kontrastlar: her biri "factor:test:ref"; her biri için ayrı dosya. İLK giriş
  # birincil deseq2_results.tsv olur (downstream/rapor birincili tüketir). condition
  # faktörü → dosya adı geriye uyumlu (öneksiz); condition-dışı → faktör-önekli.
  specs <- strsplit(contrasts_spec, ";", fixed = TRUE)[[1]]
  first <- TRUE
  for (spec in specs) {
    ftr <- strsplit(spec, ":", fixed = TRUE)[[1]]
    factor_name <- ftr[1]; test <- ftr[2]; ref <- ftr[3]
    if (!(factor_name %in% colnames(coldata))) {
      stop(sprintf("contrast factor '%s' not found in coldata columns (%s)",
                   factor_name, paste(colnames(coldata), collapse = ", ")))
    }
    flevels <- levels(coldata[[factor_name]])
    if (!(test %in% flevels) || !(ref %in% flevels)) {
      stop(sprintf("contrast '%s: %s vs %s': level not found (%s levels: %s)",
                   factor_name, test, ref, factor_name, paste(flevels, collapse = ", ")))
    }
    res_i <- results(dds, contrast = c(factor_name, test, ref))
    fname <- if (factor_name == "condition")
               sprintf("deseq2_results.%s_vs_%s.tsv", test, ref)
             else
               sprintf("deseq2_results.%s.%s_vs_%s.tsv", factor_name, test, ref)
    write_results(res_i, file.path(out_dir, fname))
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
