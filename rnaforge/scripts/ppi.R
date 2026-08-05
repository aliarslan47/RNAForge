# rnaforge/scripts/ppi.R
# Argümanlar: nodes.tsv edges.tsv out_dir
# En büyük modüllerin STRING alt-ağı: kenarlar gri çizgi, düğümler modüle göre renkli,
# hub genler (en yüksek derece) etiketli, biçim yöne (artan/azalan) göre. Boş -> boş-durum paneli.
suppressMessages({library(ggplot2)})
a <- commandArgs(trailingOnly = TRUE)
nodes_p<-a[1]; edges_p<-a[2]; out<-a[3]
dir.create(out, showWarnings=FALSE, recursive=TRUE)
sav <- function(p,name,w,h){ ggsave(file.path(out,paste0(name,".png")),p,width=w,height=h,dpi=300,limitsize=FALSE)
                             ggsave(file.path(out,paste0(name,".svg")),p,width=w,height=h,limitsize=FALSE) }
empty_panel <- function(){ ggplot()+annotate("text",x=0,y=0,label="Ağda modül yok",size=5,color="grey40")+
  labs(title="Protein etkileşim modülleri (STRING)")+theme_void(base_size=13)+theme(plot.title=element_text(face="bold")) }

nd <- tryCatch(read.delim(nodes_p, check.names=FALSE), error=function(e) NULL)
ed <- tryCatch(read.delim(edges_p, check.names=FALSE), error=function(e) NULL)
if(is.null(nd) || nrow(nd) < 2){ sav(empty_panel(), "ppi_network", 9, 7); quit(save="no") }

# hub etiketi: modül başına en yüksek dereceli 2 düğüm
nd$dir <- ifelse(nd$direction=="up","Artan",ifelse(nd$direction=="down","Azalan","ns"))
nd <- nd[order(nd$module, -nd$degree), ]
nd$lab <- ""
for(m in unique(nd$module)){ idx<-which(nd$module==m); k<-head(idx,2); nd$lab[k]<-nd$symbol[k] }

p <- ggplot()
if(!is.null(ed) && nrow(ed)>0)
  p <- p + geom_segment(data=ed, aes(x=x1,y=y1,xend=x2,yend=y2), color="grey80", linewidth=0.2, alpha=0.5)
p <- p +
  geom_point(data=nd, aes(x=x, y=y, color=module, shape=dir, size=degree), alpha=0.85)+
  ggrepel::geom_text_repel(data=nd, aes(x=x,y=y,label=lab), size=2.8, color="grey15",
                           max.overlaps=Inf, min.segment.length=0, seed=42)+
  scale_shape_manual(values=c(Artan=16, Azalan=17, ns=15), name="yön")+
  scale_size_continuous(range=c(1.5,6), guide="none")+
  scale_color_brewer(palette="Set2", name="modül")+
  labs(title="Protein etkileşim modülleri (STRING, en büyük modüller)")+
  theme_void(base_size=12)+
  theme(plot.title=element_text(face="bold"), legend.position="right",
        plot.background=element_rect(fill="white", color=NA),
        plot.margin=margin(8,8,8,8))
sav(p, "ppi_network", 10, 7.5)
cat("ppi.R done:", nrow(nd), "nodes\n")
