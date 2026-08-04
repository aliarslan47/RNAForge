# m07 — Visualization (Figürler) — Tasarım

**Tarih:** 2026-08-03
**Durum:** Onaylandı (Ali — "başla"; figür stili = statik yüksek çözünürlük seçildi)
**İlgili:** PLAN.md v1.3 §2.1, §7 (Figure 4 PCA / 5 Heatmap / 6 Volcano / 7 MA), §14; m06 spec (DE çıktı sözleşmesi)

## 1. Problem

m06 diferansiyel ekspresyonu (`deseq2_results.tsv`, `normalized_counts.tsv`) üretti — ama
sadece tablo. m07 bunu **yayın kalitesinde figürlere** çevirir. Bu, kullanıcının açık isteği:
figürler pipeline'ın kalıcı parçası olsun, **her koşuda otomatik** çıksın, **modern + yüksek
çözünürlük**. m06 gibi **organizma-agnostik** (prok/ökaryot ayrımı bilmez) — yalnız koşunun
KENDİ verisini tüketir.

**Kapsam dışı (bilinçli):** Makale/konkordans karşılaştırma figürü m07'YE GİRMEZ. O yalnız
tek seferlik doğrulama aracıydı; müşteri koşusunda karşılaştırılacak yayımlanmış referans
olmaz. Dış referans gerektiren hiçbir görsel standart çıktı değildir.

## 2. Onaylanan kararlar

1. **Figür seti (MVP): PCA · Volcano · Heatmap · MA plot.** PLAN §7'de Fig 4/5/6/7. Hepsi
   koşunun kendi verisinden; dış referans yok. GO/KEGG/GSEA/dashboard = Faz 2+ (dahil değil).
2. **Motor: R/ggplot2, izole `rnaforge-de` env** (kurulu; m06 ile aynı env — R zaten orada).
   R betiği repo'da versiyonlanır: `rnaforge/scripts/figures.R` (Python'da gömülü string DEĞİL,
   m06/`deseq2.R` deseni). Palet renk-körü güvenli, tema modern-minimal (tek yerde tanımlı).
3. **Çıktı formatı: statik yüksek çözünürlük — her figür hem PNG (300 DPI) hem SVG (vektör).**
   Ali onayı. İnteraktif (plotly) PLAN Faz 2+; MVP değil. HTML rapor (m08) PNG'leri gömer,
   SVG opsiyonel keskin baskı için durur.
4. **Girdi = m06 çıktıları + config + GFF.** `differential_expression/{normalized_counts.tsv,
   deseq2_results.tsv,coldata.tsv}`; eşikler `config.de.{fdr_threshold,log2fc_threshold}`;
   gen adları GFF'den (`quantification.attribute`→gene) — etiketleme için, YOKSA locus_tag'e düşer.
5. **Veri kapısı YOK — m07 ASLA FAIL üretmez.** Görselleştirme biyolojiyi geçersiz kılmaz
   (m06 mantığı). m07 yeni veri-kapısı yazmaz; kalite kartını yeniden yazar (önceki kapılar
   taşınır), verdict değişmez. Bir figür render EDİLEMEZSE → yüksek sesle hata + exit 1
   (sessiz boş figür YASAK, [[feedback_gurultulu_hata]]); bu bir kapı değil, sağlamlık kuralı.
6. **Ön koşul: m06 done** (`deseq2_results.tsv` gerekir). Zincir …→m05→m06→m07. Aynı `--run-id`.
7. **Çıktı dizini:** `runs/.../figures/{01_pca,02_volcano,03_heatmap,04_ma}.{png,svg}` + figür
   başına bir `manifest.json` (hangi figür, dosya, caption) — m08 raporu bunu tüketir (sözleşme).

### Reddedilen seçenekler
- **Konkordans/makale figürü standart çıktı:** dış referans gerektirir → yalnız doğrulama, RED.
- **İnteraktif plotly (MVP):** JS bağımlılığı + env kurulumu + ağır rapor; PLAN Faz 2+ → şimdi değil.
- **matplotlib/Python:** hiçbir env'de yok; R/ggplot2 zaten `rnaforge-de`'de → gereksiz yeni bağımlılık.
- **m07'ye FAIL kapısı:** görsel biyolojiyi geçersiz kılmaz → kapı yok (m06 gibi).
- **R kodu Python'da gömülü string:** okunmaz/lint'lenmez → ayrı `.R` dosyası.

## 3. Arayüz (public sözleşme)

```python
# rnaforge/figures.py — saf yardımcılar + gerçek runner (env)
def build_figures(norm_counts: Path, de_results: Path, coldata: Path,
                  gene_map: dict[str,str], fdr: float, lfc: float,
                  out_dir: Path, env: str = "rnaforge-de") -> list[Path]: ...

# rnaforge/modules/m07_figures.py
def run_figures(config: Config, metadata_path: Path, run_dir: Path,
                force: bool = False) -> dict: ...
```
- **Ön koşul:** m06 done değilse `ValueError` (aynı `--run-id` ile önce `rnaforge de` koş).
- **Dönüş:** özet dict: `n_figures`, `figures` (ad→yol), `formats` (["png","svg"]), `resumed?`.
- **CLI:** `rnaforge figures --config … --metadata … --run-id …` → `_cmd_figures`, kalite kartı
  yazar + verdict basar (m06 CLI deseni). `resumed` ise "already completed" uyarısı.

## 4. Figür şartnamesi

| # | Figür | Girdi | İçerik |
|---|-------|-------|--------|
| 01 | **PCA** | normalized_counts (log2), en değişken ~500 gen | PC1/PC2, koşula göre renk, örnek etiketi, %varyans ekseni |
| 02 | **Volcano** | deseq2_results | log2FC vs -log10 padj; Up/Down/NS renk; eşik çizgileri; en güçlü genler etiketli |
| 03 | **Heatmap** | normalized_counts, en güçlü ~40 DEG | gen z-skoru, hiyerarşik sıralı, örnek sütunları |
| 04 | **MA** | deseq2_results | log2(baseMean) vs log2FC; anlamlı vurgulu; eşik çizgisi |

- Etiketleme: `ggrepel` (temiz, çakışmasız) — `rnaforge-de`'ye eklenecek (yoksa kurulacak).
- Ortak tema/palet `figures.R` tepesinde tek tanım (tutarlılık).

## 5. Test stratejisi (TDD)

- **Saf birim:** gene_map çıkarımı (GFF→ad), en-değişken-gen seçimi, z-skor, Up/Down/NS
  sınıflama, manifest üretimi — R'sız test edilebilir Python yardımcıları.
- **Orkestrasyon (monkeypatch):** runner'ı sahteleyip `run_figures` ön koşul/resume/heartbeat/
  kart/manifest sırasını doğrula; m06 yoksa `ValueError`.
- **Gerçek-araç entegrasyon (env yoksa skip):** küçük gerçek matriste `figures.R` 4 PNG+SVG
  üretiyor mu; dosyalar var + boş değil.
- **Canlı smoke:** var olan `runs/20260803_143036_GSE300731` üstünde `rnaforge figures` →
  4 figür + manifest; verdict değişmedi (m07 kapı eklemez).

## 6. Sıradaki (bu spec DIŞI)
- **m08 — Reporting:** kendine yeten HTML (figürleri base64 gömer + kalite kartı + kapı tablosu
  + istatistik + en güçlü DEG tablosu; dil config'ten). Ayrı spec+plan, m07 biter bitmez.
