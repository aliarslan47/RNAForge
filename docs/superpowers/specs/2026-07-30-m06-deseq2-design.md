# m06 — Differential Expression (DESeq2) — Tasarım

**Tarih:** 2026-07-30
**Durum:** Onaylandı (Ali — "başla", otonom oturum; DESeq2 motoru = conda env kararı onaylı)
**İlgili:** PLAN.md v1.3 §4.1, §5, §14; m05 spec (count matrisi sözleşmesi)

## 1. Problem

m05 gen×örnek count matrisini (`counts.tsv`) üretti. m06 bunu diferansiyel ekspresyona
çevirir — **pipeline'ın gerçek biyolojik çıktısı**: hangi genler koşullar arasında anlamlı
değişiyor. İlk R/Bioconductor entegrasyonudur.

Bu, count matrisi sözleşmesinden sonra ilk **ortak** (organizma-agnostik) adımdır: prokaryot
ve ökaryot yolları burada birleşir; m06 organizma tipini bilmez.

## 2. Onaylanan kararlar

1. **Motor: R/Bioconductor DESeq2, izole conda env (`rnaforge-de`).** `Rscript` ile çağrılır.
   PLAN §4.1 kararı: DESeq2 birincil (altın standart; Katman A doğrulaması yayımlanmış sonuçlarla
   en çok onunla güven verir). **pydeseq2 çapraz-kontrolü m06'nın PER-RUN yolunda DEĞİL** — PLAN
   onu yalnız doğrulama aktivitesi olarak konumlar; ayrı iş, m06 kapsamı dışı.
2. **Girdi = m05 `counts.tsv` + metadata + design formülü.** Design `config.de.design`'dan
   (`~condition`, `~batch + condition`). coldata metadata'dan (condition, varsa batch).
3. **Kontrast: config-driven referans.** `config.de.reference` (YENİ, opsiyonel): DESeq2 factor
   referans seviyesi. Varsayılan = koşul seviyelerinin alfabetik ilki. Karşılaştırma
   "diğer seviye vs referans"; MVP iki-koşul odaklı (çok-seviyede DESeq2 `results()` varsayılanı,
   son-vs-referans). Yapılan kontrast çıktıya + özete YAZILIR (sessiz varsayım yok).
4. **Veri kapısı: `replicate_correlation` — WARN kapısı.** Koşul-içi replikaların normalize
   sayımlarının minimum ikili Pearson korelasyonu profil eşiğinin (`replicate_correlation`,
   prok 0.85) altındaysa **WARN** (sonuç ŞÜPHELİ damgalanır, ama ÜRETİLİR). Düşük replika
   korelasyonu DE'yi GEÇERSİZ kılmaz (DESeq2 dağılımı yine tahmin eder, gücü düşer) → FAIL değil
   WARN. İkili politika. Korelasyon R'da hesaplanır, kapı Python'da kurulur.
5. **Ön koşul: m05 done** (`counts.tsv` gerekir). Zincir m01→m03→m04→m05→m06.
6. **R betiği repo'da versiyonlanır** (`rnaforge/scripts/deseq2.R`) — Python'dan gömülü string
   değil; R kodu R dosyasında yaşar (okunabilir, R-linter'lanabilir).

### Reddedilen seçenekler
- **pydeseq2 birincil:** PLAN'ın "R DESeq2 birincil" kararını bozar; altın standart referansı kaybeder.
- **replicate_correlation = FAIL:** Düşük korelasyon DE'yi geçersiz kılmaz (güç meselesi) → WARN doğru.
- **R kodu Python'da gömülü string:** okunmaz/test edilmez; ayrı `.R` dosyası doğru.
- **Kontrastı sessizce tahmin:** kullanıcı hangi karşılaştırmanın yapıldığını bilmeli → rapora yazılır.

## 3. Arayüz (public sözleşme)

```python
def run_de(config: Config, metadata_path: Path, run_dir: Path, force: bool = False) -> dict: ...
```
- **Ön koşul:** m05 done değilse `ValueError`.
- **Dönüş:** özet dict: `n_genes`, `n_significant` (padj<fdr & |log2fc|≥lfc), `contrast`,
  `min_replicate_correlation`, `gate_counts`, `resumed?`.

```python
# rnaforge/deseq2.py
def run_deseq2(counts_tsv: Path, coldata_tsv: Path, design: str, out_dir: Path,
               reference: str | None = None, env: str = "rnaforge-de") -> DeseqResult: ...
def parse_deseq2_results(results_text: str) -> list[dict]: ...       # gene başına satır
def parse_de_metrics(metrics_text: str) -> dict: ...                 # key<TAB>value
```

## 4. Bileşenler

- **`rnaforge/scripts/deseq2.R`** — argümanlar: counts.tsv, coldata.tsv, design, reference,
  out_dir. `DESeqDataSetFromMatrix` → `DESeq` → `results`. Yazar:
  - `deseq2_results.tsv`: `gene baseMean log2FoldChange lfcSE stat pvalue padj`.
  - `normalized_counts.tsv`: gen×örnek normalize sayım (`counts(dds, normalized=TRUE)`).
  - `de_metrics.tsv`: `key<TAB>value` — `min_replicate_correlation`, `contrast`, `n_genes`.
  - Hata varsa stderr'e yazıp sıfırdan farklı exit (sessiz kısmi çıktı yok).
- **`rnaforge/deseq2.py`** — saf parserlar + runner (m05 deseni):
  - `parse_deseq2_results(text) -> list[dict]`: TSV → gen satırları (`NA`→None; sayısal alanlar float).
  - `parse_de_metrics(text) -> dict`: `key<TAB>value` → dict (sayısal olanlar float).
  - `run_deseq2(...)`: coldata zaten yazılmış varsayar; `conda run -n rnaforge-de Rscript
    scripts/deseq2.R ...` çağırır; üç çıktıyı parse edip `DeseqResult` döner.
  - `DeseqResult` frozen: `results: list[dict]`, `metrics: dict`, `results_path: Path`,
    `normalized_path: Path`.
  - `DeseqRunError(RuntimeError)`, `DeseqParseError(ValueError)`.
- **`rnaforge/modules/m06_de.py`** — `run_de` orkestrasyonu + `build_de_gates(min_corr: float,
  profile) -> list[GateResult]` (replicate_correlation WARN).
- **`rnaforge/config.py`** — `DE.reference: str | None = None` (opsiyonel; `de.reference`).
- **`rnaforge/cli.py`** — `de` subcommand + `_cmd_de`.

## 5. Veri akışı ve çıktı yapısı (PLAN §14)

```
runs/<ts>_<run_id>/
├── quantification/counts.tsv          # m05 girdisi (okunur)
├── differential_expression/           # DENEY-bazlı
│   ├── coldata.tsv                     # sample condition [batch] (m06 yazar, R okur)
│   ├── deseq2_results.tsv              # gen başına log2FC/padj...
│   ├── normalized_counts.tsv          # gen×örnek normalize
│   └── de_metrics.tsv                  # min_replicate_correlation, contrast, n_genes
├── statistics/de_statistics.json      # özet (n_significant, contrast, min_corr)
├── quality/gates.json                 # m06 WARN kapısı EKLENİR
└── logs/de.log
```

## 6. Hata yönetimi (Kural 7)
- **m05 done değil:** `ValueError` ("run m05 (counts) first").
- **R/DESeq2 env yok / Rscript sıfırdan farklı exit / çıktı yok:** `DeseqRunError` (cmd, exit, stderr).
- **Sonuç TSV bozuk/eksik sütun:** `DeseqParseError`.
- **Tek koşul seviyesi (karşılaştırma yok):** DESeq2 design'ı reddeder → R hata verir → `DeseqRunError`
  (aslında m01 replication/design kapıları bunu daha erken yakalar).
- **WARN yazma sırası:** m06 FAIL üretmez; yine de sıra tutarlı — results yaz → gates yaz →
  (raise gerekmez, WARN durdurmaz) → `mark_done`. Güvence kartı WARN'ı SUSPECT olarak toplar.
- **Resume:** m06 done + `de_statistics.json` varsa tekrar koşmaz; `--force`.

## 7. Test stratejisi (TDD)
**Birim (fixture/monkeypatch, env'siz):**
- `parse_deseq2_results`: fixture TSV → gen satırları; `NA`→None; sayısal alanlar float.
- `parse_de_metrics`: `key<TAB>value` → dict; sayısal float.
- `build_de_gates`: min_corr < eşik → WARN + measured; ≥ eşik → PASS; config ezimi → overridden;
  **asla FAIL değil**.
- `run_de` (run_deseq2 monkeypatch): coldata.tsv yazılır (condition[/batch]); deseq2_results/
  normalized/de_metrics çıktı yolları + de_statistics.json + gates; öncekiler korunur; m05 done
  değilse ValueError; n_significant config fdr/lfc eşikleriyle hesaplanır; resume Rscript'i tekrar
  çağırmaz; düşük korelasyon → WARN (koşu DURMAZ, exit 0).
- CLI: `de` exit 0 (validate→...→de), verdict basılır (WARN'da SUSPECT).

**Entegrasyon (gerçek R DESeq2, `rnaforge-de` yoksa skip):**
- Küçük sentetik count matrisi (2 koşul × replikalar, birkaç gen, birkaç gende açık sinyal) +
  coldata → gerçek `run_deseq2` → results.tsv parse edilebilir, sinyalli genlerin padj'ı küçük,
  metrics min_replicate_correlation makul. DESeq2'nin (rnaforge-de) gerçek çıktısıyla parser uyumu.

## 8. Kapsam dışı (sonraki spec'ler)
- **pydeseq2 çapraz kontrolü / Katman A doğrulaması:** ayrı doğrulama aktivitesi (PLAN).
- **m07 figürler (PCA/Volcano/Heatmap):** `deseq2_results.tsv` + `normalized_counts.tsv` tüketir.
- **m08 HTML rapor.**
- Çoklu-kontrast, LRT, shrinkage (apeglm/ashr), GO/KEGG (Faz 2+).
- Demo veri seti seçimi (Katman B) — kodu bloke etmez.
