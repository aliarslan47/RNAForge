# rnaforge/scripts/figures.R
# Argümanlar: normalized_counts deseq2_results coldata gene_map dispersions fdr lfc out_dir
suppressMessages({library(ggplot2)})
a <- commandArgs(trailingOnly = TRUE)
nc_p<-a[1]; de_p<-a[2]; cd_p<-a[3]; gm_p<-a[4]; disp_p<-a[5]
fdr<-as.numeric(a[6]); lfc<-as.numeric(a[7]); out<-a[8]
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
lg <- log2(as.matrix(nc)+1)

## 01 PCA
v <- apply(lg,1,var); top <- head(order(v,decreasing=TRUE),500)
pc <- prcomp(t(lg[top,])); ve <- round(100*pc$sdev^2/sum(pc$sdev^2),1)
d1 <- data.frame(PC1=pc$x[,1],PC2=pc$x[,2],condition=cond,sample=colnames(nc))
p1 <- ggplot(d1,aes(PC1,PC2,color=condition,label=sample))+geom_point(size=4)+
  geom_text(vjust=-1,size=3,show.legend=FALSE)+
  labs(title="PCA",x=paste0("PC1 (",ve[1],"%)"),y=paste0("PC2 (",ve[2],"%)"))+theme_pub
sav(p1,"01_pca",6.5,5)

## 02 Sample correlation (Pearson, hclust-ordered)
if(ncol(lg)>=2){
  cm<-suppressWarnings(cor(lg, method="pearson"))
  cm[!is.finite(cm)]<-0                          # sabit/bos ornek -> NaN; kumele/goster icin 0'a indir
  ord_ok<-tryCatch({ordc<-hclust(as.dist(1-cm))$order; TRUE}, error=function(e) FALSE)
  if(ord_ok) cm<-cm[ordc,ordc]
  longc<-data.frame(s1=factor(rep(rownames(cm),ncol(cm)),levels=rownames(cm)),
    s2=factor(rep(colnames(cm),each=nrow(cm)),levels=colnames(cm)),r=as.vector(cm))
  p2<-ggplot(longc,aes(s1,s2,fill=r))+geom_tile()+
    scale_fill_gradient(low="#f7fbff",high="#08306b",limits=c(min(cm),1),name="r")+
    labs(title="Örnek korelasyonu (Pearson)",x=NULL,y=NULL)+
    theme_minimal(base_size=10)+theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())
}else{
  p2<-ggplot()+annotate("text",x=0,y=0,label="Tek örnek — korelasyon yok",size=5,color="grey40")+
    labs(title="Örnek korelasyonu (Pearson)")+theme_void(base_size=13)+theme(plot.title=element_text(face="bold"))
}
sav(p2,"02_sample_correlation",6,5)

## 03 Expression distribution (per-sample log2 boxplot)
longe<-data.frame(sample=factor(rep(colnames(lg),each=nrow(lg)),levels=colnames(lg)),
  value=as.vector(lg), condition=rep(cond,each=nrow(lg)))
p3<-ggplot(longe,aes(sample,value,fill=condition))+geom_boxplot(outlier.size=0.3)+
  labs(title="Ekspresyon dağılımı (log2 normalize)",x=NULL,y="log2(sayım+1)")+theme_pub+
  theme(axis.text.x=element_text(angle=45,hjust=1))
sav(p3,"03_expression_dist",7,5)

## 04 Dispersion (plotDispEsts style) from dispersions.tsv
dp<-read.delim(disp_p, check.names=FALSE)
dp<-dp[is.finite(dp$baseMean)&dp$baseMean>0,,drop=FALSE]
p4<-ggplot(dp,aes(baseMean))+
  geom_point(aes(y=dispGeneEst),color="grey60",size=0.5,alpha=0.5)+
  geom_point(aes(y=dispFinal),color="#D55E00",size=0.5,alpha=0.6)
if(any(is.finite(dp$dispFit))) p4<-p4+geom_line(aes(y=dispFit),color="#0072B2",linewidth=0.9)
p4<-p4+scale_x_log10()+scale_y_log10()+
  labs(title="Dispersiyon tahmini",x="ortalama normalize sayım",y="dispersiyon")+theme_pub
sav(p4,"04_dispersion",6.5,5)

## 05 p-value distribution
pv<-de$pvalue[is.finite(de$pvalue)]
p5<-ggplot(data.frame(p=pv),aes(p))+geom_histogram(bins=40,fill="#0072B2",color="white",linewidth=0.1)+
  labs(title="p-değeri dağılımı",x="p-değeri",y="gen sayısı")+theme_pub
sav(p5,"05_pval_histogram",6.5,4.5)

## 06 Volcano
de$dir<-"NS"; de$dir[!is.na(de$padj)&de$padj<fdr&de$log2FoldChange>=lfc]<-"Up"
de$dir[!is.na(de$padj)&de$padj<fdr&de$log2FoldChange<=-lfc]<-"Down"
de$dir<-factor(de$dir,levels=c("Up","Down","NS")); de$mlp<--log10(pmax(de$padj,1e-300))
sig<-de[de$dir!="NS"&!is.na(de$padj),]; sig<-sig[order(sig$padj),]; labs<-head(sig,15)
p6 <- ggplot(de,aes(log2FoldChange,mlp,color=dir))+geom_point(size=0.7,alpha=0.6)+
  scale_color_manual(values=col_dir,name=NULL)+
  geom_vline(xintercept=c(-lfc,lfc),linetype=2,color="grey50")+
  geom_hline(yintercept=-log10(fdr),linetype=2,color="grey50")+
  ggrepel::geom_text_repel(data=labs,aes(label=lab_of(gene_id)),size=3,color="black",max.overlaps=20)+
  labs(title="Volcano",x="log2 fold change",y="-log10 padj")+theme_pub
sav(p6,"06_volcano",7,5.5)

## 07 MA
de$A<-log2(de$baseMean+1)
p7 <- ggplot(de,aes(A,log2FoldChange,color=dir))+geom_point(size=0.7,alpha=0.6)+
  scale_color_manual(values=col_dir,name=NULL)+geom_hline(yintercept=0,color="grey40")+
  geom_hline(yintercept=c(-lfc,lfc),linetype=2,color="grey60")+
  labs(title="MA plot",x="log2 baseMean",y="log2 fold change")+theme_pub
sav(p7,"07_ma",7,5)

## 08 Heatmap (top 40 DEG) — az/sifir DEG veya sifir-varyansta koru
htitle <- "En güçlü 40 DEG (z-skor)"
topg<-head(sig,40); m<-lg[topg$gene_id,,drop=FALSE]; rownames(m)<-lab_of(topg$gene_id)
z<-t(scale(t(m)))
z[!is.finite(z)]<-0
z<-z[rowSums(abs(z))>0,,drop=FALSE]
if(nrow(z)>=2){ z<-z[hclust(dist(z))$order,,drop=FALSE] }
if(nrow(z)>=1){
  long<-data.frame(gene=factor(rep(rownames(z),ncol(z)),levels=rownames(z)),
    sample=factor(rep(colnames(z),each=nrow(z)),levels=colnames(z)),z=as.vector(z))
  p8 <- ggplot(long,aes(sample,gene,fill=z))+geom_tile()+
    scale_fill_gradient2(low="#0072B2",mid="white",high="#D55E00",midpoint=0,name="z")+
    labs(title=htitle,x=NULL,y=NULL)+
    theme_minimal(base_size=10)+theme(axis.text.x=element_text(angle=45,hjust=1),panel.grid=element_blank())
}else{
  p8 <- ggplot()+annotate("text",x=0,y=0,label="Anlamlı/gösterilebilir DEG yok",size=5,color="grey40")+
    labs(title=htitle)+theme_void(base_size=13)+theme(plot.title=element_text(face="bold"))
}
sav(p8,"08_heatmap",6.5,8)
cat("figures.R done\n")
