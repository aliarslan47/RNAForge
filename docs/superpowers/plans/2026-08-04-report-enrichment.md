# Rapor & Figür Zenginleştirme Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** m06/m07/m08'e 4 yeni standart figür (örnek korelasyon, ekspresyon dağılımı, dispersiyon, p-değeri histogramı), up/down ayrı DEG tabloları ve açıklayıcı içerik (figür-altı caption + bölüm-başı intro) ekle.

**Architecture:** m06 `deseq2.R` dispersiyon verisi üretir + `m06_de.py` n_up/n_down istatistiği. m07 `figures.R`+`figures.py` 4 figürü ekleyip 8'e çıkar. m08 `report_html.py` tabloyu up/down ikiye böler ve çift dilli caption/intro basar. Saf Python + R (ggplot2); yeni bağımlılık yok.

**Tech Stack:** Python 3.11 stdlib, R/ggplot2 (`rnaforge-de`), pytest.

## Global Constraints

- Kod/log İngilizce; rapor metni/caption/intro çift dilli (`tr`|`en`). Uydurma biyolojik yorum YOK — yalnız sabit/sayısal metin.
- Yeni veri-kapısı YOK; verdict m06'dan taşınır. Self-contained tek HTML; figürler base64.
- Kenar durumları çökmez (m07 dersi): <2 örnek / 0-DEG / NA → boş-durum veya atla, exit 1 yok.
- Yeni bağımlılık eklenmez. `python -m rnaforge.cli` ÇALIŞMAZ; entry point `rnaforge`.
- Test: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q`.
- Zincir: m06→m07→m08 aynı `--run-id`. m07 dispersiyon figürü m06'nın `dispersions.tsv`'sini ister.

---

### Task 1: m06 — `deseq2.R` `dispersions.tsv` üretir

**Files:**
- Modify: `rnaforge/scripts/deseq2.R`
- Test: `tests/test_deseq2.py`

**Interfaces:**
- Produces: `differential_expression/dispersions.tsv` — sütunlar `gene_id`, `baseMean`, `dispGeneEst`, `dispFit`, `dispFinal` (dispFit fallback yolunda `NA`).

- [ ] **Step 1: Write the failing test** (mevcut env-gated entegrasyon testine assertion ekle)

```python
# tests/test_deseq2.py — test_run_deseq2_detects_signal içine, son assert'ten sonra ekle:
    disp = tmp_path / "de" / "dispersions.tsv"
    assert disp.exists()
    header = disp.read_text().splitlines()[0].split("\t")
    assert header == ["gene_id", "baseMean", "dispGeneEst", "dispFit", "dispFinal"]
```

- [ ] **Step 2: Run test to verify it fails**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_deseq2.py::test_run_deseq2_detects_signal -v`
Expected: FAIL (dispersions.tsv yok) — env varsa; env yoksa skip (o zaman Task doğrulaması Task 6 smoke'ta).

- [ ] **Step 3: Implement** — `deseq2.R`'da `normalized_counts` yazımından SONRA, `de_metrics` yazımından ÖNCE ekle:

```r
# Dispersiyon tahminleri (m07 dispersiyon figürü için). dispFit fallback yolunda olmayabilir -> NA.
mc <- mcols(dds)
dispfit <- if ("dispFit" %in% colnames(mc)) mc$dispFit else rep(NA_real_, nrow(dds))
disp <- data.frame(gene_id = rownames(dds), baseMean = mc$baseMean,
                   dispGeneEst = mc$dispGeneEst, dispFit = dispfit,
                   dispFinal = dispersions(dds))
write.table(disp, file.path(out_dir, "dispersions.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
```

- [ ] **Step 4: Run test to verify it passes**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_deseq2.py -v`
Expected: PASS (env varsa) — dispersions.tsv beklenen başlıkla üretildi.

- [ ] **Step 5: Commit**

```bash
git add rnaforge/scripts/deseq2.R tests/test_deseq2.py
git commit -m "feat(m06): deseq2.R dispersions.tsv uretir (m07 dispersiyon figuru icin)"
```

---

### Task 2: m06 — `de_statistics.json`'a `n_up`/`n_down`

**Files:**
- Modify: `rnaforge/modules/m06_de.py`
- Test: `tests/test_m06_de.py`

**Interfaces:**
- Produces: `count_up_down(results: list[dict], fdr: float, lfc: float) -> tuple[int,int]` (saf); `de_statistics.json`'a `n_up`,`n_down` eklenir. `n_up+n_down == n_significant`.

- [ ] **Step 1: Write the failing test**

```python
# tests/test_m06_de.py  (append) — import satırına count_up_down ekle
from rnaforge.modules.m06_de import count_up_down

def test_count_up_down():
    res = [
        {"padj": 1e-8, "log2FoldChange": 3.0},   # up
        {"padj": 1e-6, "log2FoldChange": -2.0},  # down
        {"padj": 1e-9, "log2FoldChange": 0.2},   # |lfc|<1 -> neither
        {"padj": 0.9,  "log2FoldChange": 5.0},   # padj high -> neither
        {"padj": None, "log2FoldChange": 4.0},   # NA -> neither
    ]
    up, down = count_up_down(res, fdr=0.05, lfc=1.0)
    assert (up, down) == (1, 1)
```

- [ ] **Step 2: Run to verify fail**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m06_de.py::test_count_up_down -v`
Expected: FAIL (ImportError: count_up_down).

- [ ] **Step 3: Implement** — `m06_de.py`'da `run_de`'den ÖNCE (modül düzeyi) helper ekle:

```python
def count_up_down(results: list[dict], fdr: float, lfc: float) -> tuple[int, int]:
    up = down = 0
    for r in results:
        p = r.get("padj"); l = r.get("log2FoldChange")
        if p is None or l is None or p >= fdr:
            continue
        if l >= lfc:
            up += 1
        elif l <= -lfc:
            down += 1
    return up, down
```

`run_de` içinde, `n_sig` hesaplandıktan sonra ekle ve summary'ye koy:

```python
        n_up, n_down = count_up_down(result.results, fdr, lfc)
```

`summary` sözlüğüne (`"n_significant": n_sig,` satırından hemen sonra):

```python
            "n_up": n_up,
            "n_down": n_down,
```

- [ ] **Step 4: Run to verify pass**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_m06_de.py -v`
Expected: PASS.

- [ ] **Step 5: Commit**

```bash
git add rnaforge/modules/m06_de.py tests/test_m06_de.py
git commit -m "feat(m06): de_statistics n_up/n_down (yon ayrimi sayilari)"
```

---

### Task 3: m07 — 4 yeni figür, `FIGURE_SPECS` 8'e, runner + `figures.R`

**Files:**
- Modify: `rnaforge/figures.py`, `rnaforge/scripts/figures.R`
- Test: `tests/test_figures.py`, `tests/test_m07_figures.py`

**Interfaces:**
- Produces: `FIGURE_SPECS` 8 eleman (anlatı sırası); `run_figures_r(...)` figures.R'a `dispersions.tsv` argümanı geçer. figures.R 8 figürü PNG300+SVG üretir.

- [ ] **Step 1: Update pure tests (figures.py)** — `tests/test_figures.py`

`test_build_manifest_ok_and_missing` zaten `FIGURE_SPECS` üzerinden döngü kuruyor; sadece sabit indeks assertion'ını güncelle:

```python
    # ESKI: assert man["figures"][0]["png"] == "01_pca.png"
    assert man["figures"][0]["png"] == "01_pca.png"
    # ESKI: json.loads(p.read_text())["figures"][1]["id"] == "volcano"  -> artik sample_correlation
    p = write_manifest(fig)
    assert json.loads(p.read_text())["figures"][1]["id"] == "sample_correlation"
```

Yeni test ekle:

```python
def test_figure_specs_has_eight_in_order():
    from rnaforge.figures import FIGURE_SPECS
    ids = [s[0] for s in FIGURE_SPECS]
    assert ids == ["pca", "sample_correlation", "expression_dist", "dispersion",
                   "pval_histogram", "volcano", "ma", "heatmap"]
```

- [ ] **Step 2: Run to verify fail**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_figures.py -v`
Expected: FAIL (FIGURE_SPECS hâlâ 4).

- [ ] **Step 3: Implement figures.py** — `FIGURE_SPECS`'i değiştir:

```python
FIGURE_SPECS: list[tuple[str, str, str]] = [
    ("pca", "01_pca", "PCA"),
    ("sample_correlation", "02_sample_correlation", "Örnek korelasyonu"),
    ("expression_dist", "03_expression_dist", "Ekspresyon dağılımı"),
    ("dispersion", "04_dispersion", "Dispersiyon"),
    ("pval_histogram", "05_pval_histogram", "p-değeri dağılımı"),
    ("volcano", "06_volcano", "Volcano"),
    ("ma", "07_ma", "MA plot"),
    ("heatmap", "08_heatmap", "Heatmap"),
]
```

`run_figures_r`'da cmd'ye `dispersions.tsv`'yi ekle (gene_map'ten SONRA, fdr'den ÖNCE):

```python
    cmd = ["conda", "run", "-n", env, "Rscript", str(_SCRIPT),
           str(de_dir / "normalized_counts.tsv"), str(de_dir / "deseq2_results.tsv"),
           str(de_dir / "coldata.tsv"), str(gene_map), str(de_dir / "dispersions.tsv"),
           str(fdr), str(lfc), str(out_dir)]
```

- [ ] **Step 4: Implement figures.R** — TÜM dosyayı aşağıdakiyle değiştir (arg sırası: nc de cd gm **disp** fdr lfc out; `lg` bir kez hesaplanır; 8 figür):

```r
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
  cm<-cor(lg, method="pearson"); ordc<-hclust(as.dist(1-cm))$order; cm<-cm[ordc,ordc]
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
```

- [ ] **Step 5: Update m07 integration test** — `tests/test_m07_figures.py` `_write_de` helper'ına dispersions.tsv ekle ve figür sayısını 8'e çıkar:

`_write_de` fonksiyonunun sonuna (coldata yazımından sonra) ekle:

```python
    with (de/"dispersions.tsv").open("w") as f:
        f.write("gene_id\tbaseMean\tdispGeneEst\tdispFit\tdispFinal\n")
        for i in range(1,61):
            f.write(f"LT_{i}\t200\t0.05\t0.04\t0.045\n")
```

`test_figures_r_renders_all` içinde `assert len(man["figures"]) == 4` → `== 8`.

- [ ] **Step 6: Run tests**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_figures.py tests/test_m07_figures.py -v`
Expected: PASS (env varsa 8 figür PNG+SVG üretilir; 0-DEG parametrize dahil çökmez).

- [ ] **Step 7: Commit**

```bash
git add rnaforge/figures.py rnaforge/scripts/figures.R tests/test_figures.py tests/test_m07_figures.py
git commit -m "feat(m07): 4 yeni figur (korelasyon/ekspresyon/dispersiyon/p-hist), FIGURE_SPECS 8, runner disp arg"
```

---

### Task 4: m08 — up/down ayrı DEG tabloları

**Files:**
- Modify: `rnaforge/report_html.py`
- Test: `tests/test_report_html.py`

**Interfaces:**
- Produces: `top_degs_by_direction(de, gene_map, fdr, lfc, direction, n=25) -> list[dict]` (direction `"Up"`/`"Down"`); `section_table` iki alt-tablo (Up top-25 + Down top-25) basar; `LABELS`'e `up_table`/`down_table` anahtarları eklenir.

- [ ] **Step 1: Write failing test** — `tests/test_report_html.py` (append)

```python
from rnaforge.report_html import top_degs_by_direction

def test_top_degs_by_direction_filters():
    de = [
        {"gene": "U1", "baseMean": 100.0, "log2FoldChange": 3.0, "padj": 1e-8},
        {"gene": "U2", "baseMean": 100.0, "log2FoldChange": 2.0, "padj": 1e-4},
        {"gene": "D1", "baseMean": 100.0, "log2FoldChange": -4.0, "padj": 1e-9},
    ]
    up = top_degs_by_direction(de, {}, 0.05, 1.0, "Up", n=25)
    down = top_degs_by_direction(de, {}, 0.05, 1.0, "Down", n=25)
    assert [r["gene"] for r in up] == ["U1", "U2"]      # padj asc, Up only
    assert [r["gene"] for r in down] == ["D1"]

def test_section_table_has_up_and_down(tmp_path):
    de = [
        {"gene": "U1", "baseMean": 100.0, "log2FoldChange": 3.0, "padj": 1e-8},
        {"gene": "D1", "baseMean": 100.0, "log2FoldChange": -4.0, "padj": 1e-9},
    ]
    h = section_table(de, {"LT": "x"}, 0.05, 1.0, LABELS["tr"])
    assert LABELS["tr"]["up_table"] in h and LABELS["tr"]["down_table"] in h
    assert "U1" in h and "D1" in h

def test_section_table_empty_both(tmp_path):
    h = section_table([], {}, 0.05, 1.0, LABELS["tr"])
    assert LABELS["tr"]["no_degs"] in h
```

Not: `section_table` imzası DEĞİŞİYOR — artık ham `de_results` + gene_map + eşikleri alır (içeride yön ayırır). `render_report` çağrısı da güncellenecek (Step 3).

- [ ] **Step 2: Run to verify fail**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_report_html.py -k "direction or up_and_down or empty_both" -v`
Expected: FAIL (import/imza yok).

- [ ] **Step 3: Implement** — `report_html.py`

`LABELS` her iki dile ekle (mevcut `"table"` anahtarının yanına):

```python
        # tr:
        "up_table": "En Güçlü 25 Artan (Up)", "down_table": "En Güçlü 25 Azalan (Down)",
        # en:
        "up_table": "Top 25 Up-regulated", "down_table": "Top 25 Down-regulated",
```

`top_degs`'ten sonra ekle:

```python
def top_degs_by_direction(de: list[dict], gene_map: dict, fdr: float, lfc: float,
                          direction: str, n: int = 25) -> list[dict]:
    rows = [r for r in top_degs(de, gene_map, fdr, lfc, n=10**9) if r["direction"] == direction]
    return rows[:n]
```

`section_table`'ı TÜMÜYLE değiştir (yeni imza + iki alt-tablo):

```python
def _deg_table(rows: list[dict], L: dict) -> str:
    body = [[r["gene"],
             f'{r["log2fc"]:.2f}' if r["log2fc"] is not None else "—",
             f'{r["padj"]:.2e}' if r["padj"] is not None else "—",
             f'{r["base_mean"]:.1f}' if r["base_mean"] is not None else "—"] for r in rows]
    return _table([L["gene"], L["log2fc"], L["padj"], L["base_mean"]], body)


def section_table(de_results: list, gene_map: dict, fdr: float, lfc: float, L: dict) -> str:
    up = top_degs_by_direction(de_results, gene_map, fdr, lfc, "Up", n=25)
    down = top_degs_by_direction(de_results, gene_map, fdr, lfc, "Down", n=25)
    if not up and not down:
        return f'<section id="table"><h2>{_esc(L["table"])}</h2><p>{_esc(L["no_degs"])}</p></section>'
    up_html = (f'<h3>{_esc(L["up_table"])}</h3>{_deg_table(up, L)}' if up
               else f'<h3>{_esc(L["up_table"])}</h3><p>{_esc(L["no_degs"])}</p>')
    down_html = (f'<h3>{_esc(L["down_table"])}</h3>{_deg_table(down, L)}' if down
                 else f'<h3>{_esc(L["down_table"])}</h3><p>{_esc(L["no_degs"])}</p>')
    return (f'<section id="table"><h2>{_esc(L["table"])}</h2>{up_html}{down_html}'
            f'<p class="note">{_esc(L["full_table_note"])}</p></section>')
```

`render_report` içinde eski `top = top_degs(...)` + `section_table(top, L)` çağrısını değiştir:

```python
    # ESKI iki satiri kaldir (top = top_degs(...) ve section_table(top, L))
    # section listesindeki section_table cagrisini sununla degistir:
        section_table(inputs["de_results"], inputs["gene_map"],
                      config.de.fdr_threshold, config.de.log2fc_threshold, L),
```

- [ ] **Step 4: Run to verify pass**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_report_html.py -v`
Expected: PASS (mevcut `test_render_report_full` de geçer — "pspA" hâlâ Up tablosunda).

- [ ] **Step 5: Commit**

```bash
git add rnaforge/report_html.py tests/test_report_html.py
git commit -m "feat(m08): DEG tablosu up/down ayri (her biri top-25)"
```

---

### Task 5: m08 — figür-altı caption + bölüm-başı intro

**Files:**
- Modify: `rnaforge/report_html.py`
- Test: `tests/test_report_html.py`

**Interfaces:**
- Produces: `FIGURE_CAPTIONS: dict[str,dict[str,str]]` (`{"tr":{id:...},"en":{...}}`); `SECTION_INTRO: dict[str,dict[str,str]]`; `section_figures` her figür altına caption basar; her `section_*` başlığından sonra intro `<p class="intro">`.

- [ ] **Step 1: Write failing test** — `tests/test_report_html.py` (append)

```python
from rnaforge.report_html import FIGURE_CAPTIONS, SECTION_INTRO

def test_caption_and_intro_bilingual():
    assert set(FIGURE_CAPTIONS["tr"]) == set(FIGURE_CAPTIONS["en"])
    assert FIGURE_CAPTIONS["tr"]["pca"] != FIGURE_CAPTIONS["en"]["pca"]
    assert set(SECTION_INTRO["tr"]) == set(SECTION_INTRO["en"])

def test_section_figures_shows_caption(tmp_path):
    fig = tmp_path / "figures"; fig.mkdir()
    (fig / "01_pca.png").write_bytes(b"\x89PNG")
    manifest = {"figures": [{"id": "pca", "title": "PCA", "png": "01_pca.png", "svg": None}]}
    h = section_figures(manifest, fig, LABELS["tr"])
    assert FIGURE_CAPTIONS["tr"]["pca"][:12] in h    # caption metni basildi

def test_render_report_has_section_intro(tmp_path):
    doc = render_report(_full_inputs(tmp_path), _cfg("tr"), version="0.1.0")
    assert 'class="intro"' in doc and SECTION_INTRO["tr"]["confidence"][:12] in doc
```

Not: `section_figures` imzası DEĞİŞMEZ ama caption için dile ihtiyacı var — `L` zaten geliyor; dili `L` üzerinden değil, ayrı `lang` ile eşlemek yerine caption'ı `L` ile aynı dilde seçmek için `section_figures`'a `lang` parametresi eklenir (varsayılan "tr"), `render_report` `lang` geçer.

- [ ] **Step 2: Run to verify fail**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest tests/test_report_html.py -k "caption or intro" -v`
Expected: FAIL.

- [ ] **Step 3: Implement** — `report_html.py`

`LABELS` tanımından sonra ekle:

```python
FIGURE_CAPTIONS: dict[str, dict[str, str]] = {
    "tr": {
        "pca": "En değişken 500 genin ana bileşen izdüşümü. Aynı koşulun replikaları kümelenmeli; koşullar ayrışmalı.",
        "sample_correlation": "Örnekler arası Pearson korelasyonu (log2 normalize). Yüksek blok = tutarlı replikalar; sapan örnek burada görünür.",
        "expression_dist": "Örnek başına log2 normalize ekspresyon dağılımı. Kutular benzer olmalı; normalizasyonun dengeli olduğunu gösterir.",
        "dispersion": "DESeq2 dispersiyon tahmini: gen-bazlı (gri) vs uyum (mavi) vs son (turuncu). Nokta bulutu uyum eğrisine çökmeli.",
        "pval_histogram": "Ham p-değerlerinin dağılımı. Düz + 0'a yakın tepe sağlıklıdır; anormal biçim model/veri sorununa işaret eder.",
        "volcano": "log2 kat değişimi vs -log10 padj. Sağ üst = anlamlı artan, sol üst = anlamlı azalan genler.",
        "ma": "Ortalama ekspresyon vs log2 kat değişimi. Anlamlı genler renkli; düşük sayımda dağılım genişler.",
        "heatmap": "En güçlü 40 DEG'in örnek-başı z-skoru. Koşullar arası zıt renk blokları beklenir.",
    },
    "en": {
        "pca": "Principal-component projection of the 500 most variable genes. Replicates should cluster; conditions should separate.",
        "sample_correlation": "Between-sample Pearson correlation (log2 normalized). High blocks = consistent replicates; an outlier stands out here.",
        "expression_dist": "Per-sample log2 normalized expression distribution. Boxes should be similar, indicating balanced normalization.",
        "dispersion": "DESeq2 dispersion estimates: gene-wise (grey) vs fit (blue) vs final (orange). The cloud should shrink toward the fit.",
        "pval_histogram": "Distribution of raw p-values. Flat with a peak near 0 is healthy; an odd shape signals a model/data issue.",
        "volcano": "log2 fold change vs -log10 padj. Top-right = significant up, top-left = significant down genes.",
        "ma": "Mean expression vs log2 fold change. Significant genes coloured; spread widens at low counts.",
        "heatmap": "Per-sample z-scores of the top 40 DEGs. Contrasting colour blocks between conditions are expected.",
    },
}

SECTION_INTRO: dict[str, dict[str, str]] = {
    "tr": {
        "confidence": "Bu koşulun kalite kapılarının özeti. FAIL varsa sonuç geçersizdir; WARN varsa sonuç şüpheli damgalıdır.",
        "dataset": "Analiz edilen organizma, platform, deney tasarımı ve örnekler.",
        "quality": "Okuma işleme, hizalama ve gene atama oranları — verinin analize uygunluğu.",
        "de": "Koşullar arası diferansiyel ekspresyon özeti (DESeq2).",
        "figures": "Kalite, model tanısı ve sonuç görselleri. Her figürün altında nasıl okunacağı açıklanmıştır.",
        "table": "İstatistiksel eşiği geçen en güçlü artan ve azalan genler.",
        "methods": "Kullanılan araçlar ve parametreler.",
        "references": "Yöntemlerin dayandığı yayınlar.",
    },
    "en": {
        "confidence": "Summary of this run's quality gates. A FAIL invalidates the result; a WARN stamps it as suspect.",
        "dataset": "Organism, platform, experimental design and samples analysed.",
        "quality": "Read processing, alignment and gene-assignment rates — the data's fitness for analysis.",
        "de": "Summary of differential expression between conditions (DESeq2).",
        "figures": "Quality, model-diagnostic and result figures. Each figure includes how to read it.",
        "table": "The strongest up- and down-regulated genes passing the statistical threshold.",
        "methods": "Tools and parameters used.",
        "references": "Publications the methods are based on.",
    },
}


def _intro(section_id: str, L: dict) -> str:
    # L is a language dict; find which language by identity to pick the intro set.
    lang = "en" if L is LABELS["en"] else "tr"
    text = SECTION_INTRO[lang].get(section_id, "")
    return f'<p class="intro">{_esc(text)}</p>' if text else ""
```

Her `section_*` fonksiyonunda `<h2>...</h2>` hemen ardına `{_intro("<id>", L)}` ekle. Örnek (`section_confidence`):

```python
        f'<section id="confidence"><h2>{_esc(L["confidence"])}</h2>{_intro("confidence", L)}'
```

Aynısını `section_dataset`("dataset"), `section_quality`("quality"), `section_de`("de"), `section_figures`("figures"), `section_table`("table"), `section_methods`("methods"), `section_references`("references") için yap.

`section_figures`'ı caption basacak şekilde güncelle (imzaya `lang="tr"` ekle):

```python
def section_figures(figures_manifest: dict, figures_dir: Path, L: dict, lang: str = "tr") -> str:
    figures_dir = Path(figures_dir)
    caps = FIGURE_CAPTIONS.get(lang, FIGURE_CAPTIONS["tr"])
    blocks = []
    for fig in figures_manifest.get("figures", []):
        png = figures_dir / fig["png"]
        if not png.exists():
            raise FileNotFoundError(f"m08: figure PNG missing for report: {png}")
        cap = caps.get(fig.get("id"), "")
        cap_html = f'<figcaption><strong>{_esc(fig.get("title"))}</strong> — {_esc(cap)}</figcaption>' \
                   if cap else f'<figcaption>{_esc(fig.get("title"))}</figcaption>'
        blocks.append(f'<figure><img src="{embed_png(png)}" alt="{_esc(fig.get("title"))}"/>{cap_html}</figure>')
    return f'<section id="figures"><h2>{_esc(L["figures"])}</h2>{_intro("figures", L)}{"".join(blocks)}</section>'
```

`render_report` içindeki `section_figures(...)` çağrısına `lang` geçir:

```python
        section_figures(inputs["figures"], inputs["figures_dir"], L, lang),
```

`_CSS`'e intro stili ekle (mevcut `.note` satırının yanına):

```python
.intro{color:#444;font-size:.95rem;margin:.2rem 0 .6rem}
```

- [ ] **Step 4: Run to verify pass + full suite**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q`
Expected: tüm testler PASS.

- [ ] **Step 5: Commit**

```bash
git add rnaforge/report_html.py tests/test_report_html.py
git commit -m "feat(m08): figur-alti caption + bolum-basi intro (cift dilli, sabit metin)"
```

---

### Task 6: Canlı smoke + branş doğrulaması

**Files:** (yalnız doğrulama)

- [ ] **Step 1: Gerçek GSE300731 zincirini yeniden koş** (m06 dispersions.tsv için --force, sonra m07, m08)

```bash
cd /home/ali/rnaforge-pipeline
conda run -n rnaforge-core rnaforge de       --config raw/GSE300731/config.yaml --metadata raw/GSE300731/metadata.tsv --run-id GSE300731 --runs-dir runs --force
conda run -n rnaforge-core rnaforge figures  --config raw/GSE300731/config.yaml --metadata raw/GSE300731/metadata.tsv --run-id GSE300731 --runs-dir runs --force
conda run -n rnaforge-core rnaforge report   --config raw/GSE300731/config.yaml --metadata raw/GSE300731/metadata.tsv --run-id GSE300731 --runs-dir runs --force
```
Expected: her adım OK; `report OK: .../report/report.html`; verdict SUSPECT (değişmez).

- [ ] **Step 2: Çıktıyı gözle doğrula**

```bash
ls runs/*/figures/            # 8 png + 8 svg + manifest.json + gene_map.tsv
grep -o 'data:image/png;base64' runs/*/report/report.html | wc -l    # 8
grep -o 'dispersions.tsv\|dispGeneEst' runs/*/differential_expression/dispersions.tsv | head -1
python3 -c "import json;d=json.load(open([p for p in __import__('glob').glob('runs/*/statistics/de_statistics.json')][0]));print('n_up',d['n_up'],'n_down',d['n_down'],'n_sig',d['n_significant'])"
```
Expected: 8 gömülü figür; dispersions.tsv var; n_up+n_down==n_significant. Chrome headless ile PNG görüntüsü alıp up/down tabloları + caption/intro + 4 yeni figürü gözle doğrula.

- [ ] **Step 3: Tüm test sayısı + branş bitir**

Run: `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q` → yeşil sayı.
`superpowers:requesting-code-review` whole-branch; temizse `superpowers:finishing-a-development-branch` ile `feat/report-enrichment` → `main` merge, DURUM güncelle.

---

## Self-Review

- **Spec coverage:** dispersions.tsv → Task 1. n_up/n_down → Task 2. 4 yeni figür + FIGURE_SPECS 8 + runner arg → Task 3. up/down ayrı tablo → Task 4. figür caption + bölüm intro → Task 5. Smoke/merge → Task 6. Kenar (0-DEG/<2 örnek/NA) → Task 3 (figures.R korumaları). Uydurma yorum yok → Task 5 (sabit metin). ✓
- **Placeholder scan:** Tüm adımlarda gerçek kod. ✓
- **Type consistency:** `run_figures_r` yeni arg sırası (nc,de,cd,gm,disp,fdr,lfc,out) figures.R `a[1..8]` ile birebir. `FIGURE_SPECS` id'leri manifest + FIGURE_CAPTIONS anahtarlarıyla birebir (pca/sample_correlation/expression_dist/dispersion/pval_histogram/volcano/ma/heatmap). `section_table` yeni imzası (de_results,gene_map,fdr,lfc,L) render_report çağrısıyla birebir. `section_figures` `lang` parametresi render_report'tan geçer. `count_up_down` → m06 summary. `top_degs_by_direction` → section_table. ✓
- **Not:** `section_table`/`section_figures` imza değişiyor → mevcut testler (Task 4/5 Step 1'de) güncelleniyor; eski `test_section_table_empty_note`/`test_section_figures_embeds` çağrıları yeni imzaya uyarlanmalı (Task 4/5 aynı dosyada düzeltir).
