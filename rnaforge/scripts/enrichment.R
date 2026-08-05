# rnaforge/scripts/enrichment.R
# Argümanlar: enrichment_up.tsv enrichment_down.tsv out_dir top_n
# Her set için namespace başına top-N zenginleşmiş GO terim dot-plot'u. Boş -> boş-durum paneli.
# Layout: uzun terim etiketleri sarılır, panel geniş, yükseklik satır sayısına göre dinamik.
suppressMessages({library(ggplot2)})
a <- commandArgs(trailingOnly = TRUE)
up_p<-a[1]; down_p<-a[2]; out<-a[3]; top_n<-as.integer(a[4])
dir.create(out, showWarnings=FALSE, recursive=TRUE)

# Uzun etiketleri satırlara böl (base R; stringr bağımlılığı yok).
wrap_lab <- function(x, width=48)
  vapply(x, function(s) paste(strwrap(s, width=width), collapse="\n"), character(1))

sav <- function(p,name,w,h){ ggsave(file.path(out,paste0(name,".png")),p,width=w,height=h,dpi=300,limitsize=FALSE)
                             ggsave(file.path(out,paste0(name,".svg")),p,width=w,height=h,limitsize=FALSE) }
empty_panel <- function(title){
  ggplot()+annotate("text",x=0,y=0,label="Anlamlı zenginleşme yok",size=5,color="grey40")+
    labs(title=title)+theme_void(base_size=13)+theme(plot.title=element_text(face="bold"))
}

# Bir zenginleştirme TSV'sini oku, p_adj<0.05 süz, namespace başına top_n al, dot-plot çiz + kaydet.
render_set <- function(path, title, base){
  df <- tryCatch(read.delim(path, check.names=FALSE), error=function(e) NULL)
  if(is.null(df) || nrow(df)==0 || !"p_adj" %in% names(df)){ sav(empty_panel(title), base, 8, 4); return(invisible()) }
  df <- df[is.finite(df$p_adj) & df$p_adj < 0.05, , drop=FALSE]
  if(nrow(df)==0){ sav(empty_panel(title), base, 8, 4); return(invisible()) }
  # namespace başına en anlamlı top_n
  df <- do.call(rbind, lapply(split(df, df$namespace), function(g){
    g <- g[order(g$p_adj, g$p_value), , drop=FALSE]; head(g, top_n) }))
  # namespace okunur ada dönüşsün (facet şeridi); sıra BP, MF, CC
  ns_full <- c(BP="Biyolojik süreç (BP)", MF="Moleküler işlev (MF)", CC="Hücresel bileşen (CC)")
  df$ns <- factor(ifelse(df$namespace %in% names(ns_full), ns_full[df$namespace], df$namespace),
                  levels=ns_full)
  # etiketleri sar; y sırası fold'a göre (büyük üstte)
  df$label <- wrap_lab(df$term)
  df$label <- factor(df$label, levels=rev(df$label[order(df$fold_enrichment)]))
  p <- ggplot(df, aes(x=fold_enrichment, y=label, size=study_count, color=p_adj))+
    geom_point()+
    facet_grid(ns~., scales="free_y", space="free")+
    scale_x_continuous(limits=c(0, NA), expand=expansion(mult=c(0.02, 0.08)))+
    scale_color_gradient(low="#D55E00", high="#4575b4", name="padj")+
    scale_size_continuous(range=c(1.5, 6), name="gen sayısı")+
    labs(title=title, x="kat-zenginleşme (fold)", y=NULL)+
    theme_minimal(base_size=11)+
    theme(panel.grid.minor=element_blank(),
          panel.grid.major.y=element_line(color="grey92", linewidth=0.3),
          plot.title=element_text(face="bold"),
          axis.text.y=element_text(size=8, lineheight=0.85),
          strip.text.y=element_text(angle=0, face="bold", size=9),
          legend.position="right")
  # yükseklik satır sayısına göre dinamik (sarılmış etiketler için satır başına ~0.34")
  h <- max(4.5, nrow(df) * 0.34 + 1.6)
  sav(p, base, 10, h)
}

render_set(up_p,   "GO zenginleştirme — Artan (Up)",   "enrichment_up")
render_set(down_p, "GO zenginleştirme — Azalan (Down)", "enrichment_down")
cat("enrichment.R done\n")
