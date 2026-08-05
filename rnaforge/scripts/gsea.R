# rnaforge/scripts/gsea.R
# Argümanlar: rnk gmt gene_map out_dir collection min_size max_size title
# fgsea (preranked) + işaretli NES dot-plot. NES>0 artan tarafta, NES<0 azalan tarafta.
suppressMessages({library(fgsea); library(ggplot2)})
a <- commandArgs(trailingOnly = TRUE)
rnk_p<-a[1]; gmt_p<-a[2]; gm_p<-a[3]; out<-a[4]; coll<-a[5]
min_size<-as.integer(a[6]); max_size<-as.integer(a[7]); title<-a[8]
dir.create(out, showWarnings=FALSE, recursive=TRUE)

wrap_lab <- function(x, width=48)
  vapply(x, function(s) paste(strwrap(s, width=width), collapse="\n"), character(1))
sav <- function(p,name,w,h){ ggsave(file.path(out,paste0(name,".png")),p,width=w,height=h,dpi=300,limitsize=FALSE)
                             ggsave(file.path(out,paste0(name,".svg")),p,width=w,height=h,limitsize=FALSE) }
empty_panel <- function(title){
  ggplot()+annotate("text",x=0,y=0,label="Anlamlı gen-seti yok",size=5,color="grey40")+
    labs(title=title)+theme_void(base_size=13)+theme(plot.title=element_text(face="bold")) }

# Ranked liste -> adlandırılmış vektör (azalan)
rk <- read.delim(rnk_p, header=FALSE, col.names=c("gene","stat"))
stats <- setNames(as.numeric(rk$stat), rk$gene)
stats <- sort(stats[is.finite(stats)], decreasing=TRUE)

# GMT -> pathways listesi + id->ad
gmt_lines <- readLines(gmt_p, warn=FALSE)
gmt_lines <- gmt_lines[nzchar(gmt_lines)]
pathways <- list(); pname <- c()
for(ln in gmt_lines){
  parts <- strsplit(ln, "\t", fixed=TRUE)[[1]]
  if(length(parts) < 3) next
  id <- parts[1]; pathways[[id]] <- parts[-c(1,2)]; pname[id] <- parts[2]
}

# gen_map: locus_tag -> sembol (öncü genler için)
gm <- tryCatch(read.delim(gm_p, check.names=FALSE), error=function(e) data.frame(locus_tag=character(),gene=character()))
sym_of <- function(ids){ v <- gm$gene[match(ids, gm$locus_tag)]; ifelse(is.na(v)|v=="", ids, v) }

tsv <- file.path(out, paste0("gsea_", coll, ".tsv"))
if(length(pathways)==0){
  writeLines("pathway_id\tname\tsize\tES\tNES\tpval\tpadj\tleading_edge", tsv)
  sav(empty_panel(title), paste0("gsea_", coll), 8, 4)
  cat("gsea.R: no pathways\n"); quit(save="no", status=0)
}

res <- fgsea(pathways=pathways, stats=stats, minSize=min_size, maxSize=max_size)
if(nrow(res)==0){
  writeLines("pathway_id\tname\tsize\tES\tNES\tpval\tpadj\tleading_edge", tsv)
  sav(empty_panel(title), paste0("gsea_", coll), 8, 4)
  cat("gsea.R: no gene sets in size range\n"); quit(save="no", status=0)
}
res <- res[order(-res$NES), ]
le <- vapply(res$leadingEdge, function(g) paste(sym_of(g), collapse=";"), character(1))
out_df <- data.frame(
  pathway_id=res$pathway, name=pname[res$pathway], size=res$size,
  ES=round(res$ES,4), NES=round(res$NES,4), pval=signif(res$pval,4),
  padj=signif(res$padj,4), leading_edge=le, stringsAsFactors=FALSE)
write.table(out_df, tsv, sep="\t", quote=FALSE, row.names=FALSE)

# İşaretli NES dot-plot: anlamlı (padj<0.05); yoksa boş-durum
sig <- out_df[is.finite(out_df$padj) & out_df$padj < 0.05, , drop=FALSE]
if(nrow(sig)==0){
  sav(empty_panel(title), paste0("gsea_", coll), 8, 4)
} else {
  # her yönden en fazla 20 terim (okunurluk)
  sig <- sig[order(-abs(sig$NES)), , drop=FALSE]
  sig <- rbind(head(sig[sig$NES>0,,drop=FALSE],20), head(sig[sig$NES<0,,drop=FALSE],20))
  sig$label <- factor(wrap_lab(sig$name), levels=rev(wrap_lab(sig$name)[order(sig$NES)]))
  p <- ggplot(sig, aes(x=NES, y=label, color=padj, size=size))+
    geom_vline(xintercept=0, color="grey70")+ geom_point()+
    scale_color_gradient(low="#D55E00", high="#4575b4", name="padj")+
    scale_size_continuous(range=c(1.5,6), name="set boyutu")+
    labs(title=title, x="NES (işaretli)", y=NULL)+
    theme_minimal(base_size=11)+
    theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold"),
          axis.text.y=element_text(size=8, lineheight=0.85), legend.position="right")
  h <- max(4.5, nrow(sig)*0.34 + 1.6)
  sav(p, paste0("gsea_", coll), 10, h)
}
cat("gsea.R done:", coll, nrow(out_df), "sets\n")
