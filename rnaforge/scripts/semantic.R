# rnaforge/scripts/semantic.R
# Argümanlar: dist.tsv reduced.tsv out_dir basename title
# REVIGO-benzeri semantik-uzay scatter: Lin uzaklık matrisi -> cmdscale(2D) -> ggplot.
# Nokta = temsilci GO terimi; renk = namespace; boyut = temsil ettiği terim sayısı; etiket = en büyükler.
suppressMessages({library(ggplot2)})
a <- commandArgs(trailingOnly = TRUE)
dist_p<-a[1]; red_p<-a[2]; out<-a[3]; base<-a[4]; title<-a[5]
dir.create(out, showWarnings=FALSE, recursive=TRUE)
sav <- function(p,name,w,h){ ggsave(file.path(out,paste0(name,".png")),p,width=w,height=h,dpi=300,limitsize=FALSE)
                             ggsave(file.path(out,paste0(name,".svg")),p,width=w,height=h,limitsize=FALSE) }
empty_panel <- function(t){ ggplot()+annotate("text",x=0,y=0,label="Yeterli terim yok",size=5,color="grey40")+
  labs(title=t)+theme_void(base_size=13)+theme(plot.title=element_text(face="bold")) }

dm <- tryCatch(read.delim(dist_p, row.names=1, check.names=FALSE), error=function(e) NULL)
red <- tryCatch(read.delim(red_p, check.names=FALSE), error=function(e) NULL)
if(is.null(dm) || nrow(dm) < 3 || is.null(red)){ sav(empty_panel(title), base, 8, 5); quit(save="no") }

mds <- tryCatch(cmdscale(as.dist(as.matrix(dm)), k=2), error=function(e) NULL)
if(is.null(mds) || ncol(mds) < 2){ sav(empty_panel(title), base, 8, 5); quit(save="no") }
d <- data.frame(go_id=rownames(dm), x=mds[,1], y=mds[,2], stringsAsFactors=FALSE)
d <- merge(d, red[,c("go_id","namespace","term","n_collapsed")], by="go_id")
ns_full <- c(BP="Biyolojik süreç", MF="Moleküler işlev", CC="Hücresel bileşen")
d$ns <- ifelse(d$namespace %in% names(ns_full), ns_full[d$namespace], d$namespace)
# en büyük (en çok terim temsil eden) ilk 12'yi etiketle
d <- d[order(-d$n_collapsed), ]
d$lab <- ""; d$lab[1:min(12,nrow(d))] <- d$term[1:min(12,nrow(d))]

p <- ggplot(d, aes(x=x, y=y, color=ns, size=n_collapsed))+
  geom_point(alpha=0.8)+
  ggrepel::geom_text_repel(aes(label=lab), size=3, color="grey20", max.overlaps=Inf,
                           min.segment.length=0, seed=42, show.legend=FALSE)+
  scale_size_continuous(range=c(2,9), name="temsil ettiği terim")+
  scale_color_brewer(palette="Dark2", name=NULL)+
  labs(title=title, x="semantik boyut 1", y="semantik boyut 2")+
  theme_minimal(base_size=12)+
  theme(panel.grid.minor=element_blank(), plot.title=element_text(face="bold"),
        axis.text=element_blank(), legend.position="right")
sav(p, base, 9, 6.5)
cat("semantic.R done:", base, nrow(d), "terms\n")
