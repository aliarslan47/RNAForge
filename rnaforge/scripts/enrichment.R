# rnaforge/scripts/enrichment.R
# Argümanlar: enrichment_up.tsv enrichment_down.tsv out_dir top_n
# Her set için top-N zenginleşmiş GO terim dot-plot'u (namespace facet). Boş -> boş-durum paneli.
suppressMessages({library(ggplot2)})
a <- commandArgs(trailingOnly = TRUE)
up_p<-a[1]; down_p<-a[2]; out<-a[3]; top_n<-as.integer(a[4])
dir.create(out, showWarnings=FALSE, recursive=TRUE)

sav <- function(p,name,w,h){ ggsave(file.path(out,paste0(name,".png")),p,width=w,height=h,dpi=300)
                             ggsave(file.path(out,paste0(name,".svg")),p,width=w,height=h) }
empty_panel <- function(title){
  ggplot()+annotate("text",x=0,y=0,label="Anlamlı zenginleşme yok",size=5,color="grey40")+
    labs(title=title)+theme_void(base_size=13)+theme(plot.title=element_text(face="bold"))
}

# Bir zenginleştirme TSV'sini oku, p_adj<0.05 süz, namespace başına top_n al, dot-plot çiz.
plot_set <- function(path, title){
  df <- tryCatch(read.delim(path, check.names=FALSE), error=function(e) NULL)
  if(is.null(df) || nrow(df)==0 || !"p_adj" %in% names(df)) return(empty_panel(title))
  df <- df[is.finite(df$p_adj) & df$p_adj < 0.05, , drop=FALSE]
  if(nrow(df)==0) return(empty_panel(title))
  # namespace başına en anlamlı top_n
  df <- do.call(rbind, lapply(split(df, df$namespace), function(g){
    g <- g[order(g$p_adj, g$p_value), , drop=FALSE]; head(g, top_n) }))
  # terim etiketi okunur ve benzersiz olsun (aynı ad iki namespace'te olabilir)
  df$label <- factor(df$term, levels=rev(df$term[order(df$fold_enrichment)]))
  ggplot(df, aes(x=fold_enrichment, y=label, size=study_count, color=p_adj))+
    geom_point()+
    facet_grid(namespace~., scales="free_y", space="free")+
    scale_color_gradient(low="#D55E00", high="#4575b4", name="padj")+
    scale_size_continuous(name="gen")+
    labs(title=title, x="kat-zenginleşme (fold)", y=NULL)+
    theme_minimal(base_size=12)+
    theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold"),
          strip.text.y=element_text(angle=0), legend.position="right")
}

sav(plot_set(up_p,   "GO zenginleştirme — Artan (Up)"),   "enrichment_up",   7.5, 6)
sav(plot_set(down_p, "GO zenginleştirme — Azalan (Down)"), "enrichment_down", 7.5, 6)
cat("enrichment.R done\n")
