# rnaforge/scripts/operon.R
# Argümanlar: operons.tsv out_dir [top_n]
# En güçlü koordineli operonlar: yatay bar (mean log2FC, yöne göre renkli). Boş -> boş-durum paneli.
suppressMessages({library(ggplot2)})
a <- commandArgs(trailingOnly = TRUE)
op_p<-a[1]; out<-a[2]; top_n<-ifelse(length(a)>=3, as.integer(a[3]), 20L)
dir.create(out, showWarnings=FALSE, recursive=TRUE)
sav <- function(p,name,w,h){ ggsave(file.path(out,paste0(name,".png")),p,width=w,height=h,dpi=300,limitsize=FALSE)
                             ggsave(file.path(out,paste0(name,".svg")),p,width=w,height=h,limitsize=FALSE) }
empty_panel <- function(t){ ggplot()+annotate("text",x=0,y=0,label="Koordineli operon yok",size=5,color="grey40")+
  labs(title=t)+theme_void(base_size=13)+theme(plot.title=element_text(face="bold")) }

df <- tryCatch(read.delim(op_p, check.names=FALSE), error=function(e) NULL)
title <- "Koordineli operonlar (operon-düzeyi DE)"
if(is.null(df) || nrow(df)==0 || !"coordinated" %in% names(df)){ sav(empty_panel(title),"operon_coord",8,4); quit(save="no") }
df <- df[df$coordinated=="yes" & is.finite(df$mean_log2fc), , drop=FALSE]
if(nrow(df)==0){ sav(empty_panel(title),"operon_coord",8,4); quit(save="no") }
df <- head(df[order(-abs(df$mean_log2fc)), , drop=FALSE], top_n)
# etiket: ilk 4 üye + … (uzun operon adları kısalsın)
lab_of <- function(g){ v<-strsplit(g,";")[[1]]; if(length(v)>4) paste0(paste(v[1:4],collapse=";"),"…") else paste(v,collapse=";") }
df$label <- vapply(df$genes, lab_of, character(1))
df$label <- factor(df$label, levels=rev(df$label[order(df$mean_log2fc)]))
df$dir <- ifelse(df$n_up>=df$n_down,"Artan","Azalan")
p <- ggplot(df, aes(x=mean_log2fc, y=label, fill=dir))+
  geom_col()+ geom_vline(xintercept=0, color="grey70")+
  scale_fill_manual(values=c(Artan="#D55E00", Azalan="#0072B2"), name=NULL)+
  labs(title=title, x="operon ortalama log2FC", y=NULL)+
  theme_minimal(base_size=11)+
  theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold"),
        axis.text.y=element_text(size=8), legend.position="top")
h <- max(4, nrow(df)*0.32 + 1.4)
sav(p, "operon_coord", 9, h)
cat("operon.R done:", nrow(df), "operons\n")
