# rnaforge/scripts/figures.R
# Argümanlar: normalized_counts.tsv deseq2_results.tsv coldata.tsv gene_map.tsv fdr lfc out_dir
suppressMessages({library(ggplot2)})
a <- commandArgs(trailingOnly = TRUE)
nc_p<-a[1]; de_p<-a[2]; cd_p<-a[3]; gm_p<-a[4]; fdr<-as.numeric(a[5]); lfc<-as.numeric(a[6]); out<-a[7]
dir.create(out, showWarnings=FALSE, recursive=TRUE)
theme_pub <- theme_minimal(base_size=13) + theme(panel.grid.minor=element_blank(),
  plot.title=element_text(face="bold"), legend.position="top")
col_dir <- c(Up="#D55E00", Down="#0072B2", NS="grey80")
sav <- function(p,name,w,h){ ggsave(file.path(out,paste0(name,".png")),p,width=w,height=h,dpi=300)
                             ggsave(file.path(out,paste0(name,".svg")),p,width=w,height=h) }
nc <- read.delim(nc_p, row.names=1, check.names=FALSE)
de <- read.delim(de_p, check.names=FALSE); names(de)[1] <- "gene_id"
cd <- read.delim(cd_p, row.names=1, check.names=FALSE)
gm <- tryCatch(read.delim(gm_p, check.names=FALSE), error=function(e) data.frame(locus_tag=character(),gene=character()))
lab_of <- function(ids){ v <- gm$gene[match(ids, gm$locus_tag)]; ifelse(is.na(v)|v=="", ids, v) }
cond <- cd[colnames(nc), "condition"]

## PCA
lg <- log2(as.matrix(nc)+1); v <- apply(lg,1,var); top <- head(order(v,decreasing=TRUE),500)
pc <- prcomp(t(lg[top,])); ve <- round(100*pc$sdev^2/sum(pc$sdev^2),1)
d1 <- data.frame(PC1=pc$x[,1],PC2=pc$x[,2],condition=cond,sample=colnames(nc))
p1 <- ggplot(d1,aes(PC1,PC2,color=condition,label=sample))+geom_point(size=4)+
  geom_text(vjust=-1,size=3,show.legend=FALSE)+
  labs(title="PCA",x=paste0("PC1 (",ve[1],"%)"),y=paste0("PC2 (",ve[2],"%)"))+theme_pub
sav(p1,"01_pca",6.5,5)

## Volcano
de$dir<-"NS"; de$dir[!is.na(de$padj)&de$padj<fdr&de$log2FoldChange>=lfc]<-"Up"
de$dir[!is.na(de$padj)&de$padj<fdr&de$log2FoldChange<=-lfc]<-"Down"
de$dir<-factor(de$dir,levels=c("Up","Down","NS")); de$mlp<--log10(pmax(de$padj,1e-300))
sig<-de[de$dir!="NS"&!is.na(de$padj),]; sig<-sig[order(sig$padj),]; labs<-head(sig,15)
p2 <- ggplot(de,aes(log2FoldChange,mlp,color=dir))+geom_point(size=0.7,alpha=0.6)+
  scale_color_manual(values=col_dir,name=NULL)+
  geom_vline(xintercept=c(-lfc,lfc),linetype=2,color="grey50")+
  geom_hline(yintercept=-log10(fdr),linetype=2,color="grey50")+
  ggrepel::geom_text_repel(data=labs,aes(label=lab_of(gene_id)),size=3,color="black",max.overlaps=20)+
  labs(title="Volcano",x="log2 fold change",y="-log10 padj")+theme_pub
sav(p2,"02_volcano",7,5.5)

## Heatmap (top 40 DEG)
topg<-head(sig,40); m<-lg[topg$gene_id,,drop=FALSE]; rownames(m)<-lab_of(topg$gene_id)
z<-t(scale(t(m))); ord<-hclust(dist(z))$order
long<-data.frame(gene=factor(rep(rownames(z)[ord],ncol(z)),levels=rownames(z)[ord]),
  sample=factor(rep(colnames(z),each=nrow(z)),levels=colnames(z)),z=as.vector(z[ord,]))
p3 <- ggplot(long,aes(sample,gene,fill=z))+geom_tile()+
  scale_fill_gradient2(low="#0072B2",mid="white",high="#D55E00",midpoint=0,name="z")+
  labs(title="En güçlü 40 DEG (z-skor)",x=NULL,y=NULL)+
  theme_minimal(base_size=10)+theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())
sav(p3,"03_heatmap",6.5,8)

## MA
de$A<-log2(de$baseMean+1)
p4 <- ggplot(de,aes(A,log2FoldChange,color=dir))+geom_point(size=0.7,alpha=0.6)+
  scale_color_manual(values=col_dir,name=NULL)+geom_hline(yintercept=0,color="grey40")+
  geom_hline(yintercept=c(-lfc,lfc),linetype=2,color="grey60")+
  labs(title="MA plot",x="log2 baseMean",y="log2 fold change")+theme_pub
sav(p4,"04_ma",7,5)
cat("figures.R done\n")
