# DURUM — RNAForge

> Bu dosya "nerede kaldık" anlık görüntüsüdür. Tüm karar detayı Claude belleğindedir
> (`rnaforge-project` memory). Claude bunu anlamlı her durakta ve `/clear` öncesi günceller.

**Konum:** `/home/ali/rnaforge-pipeline/` (git deposu)
**GitHub:** `github.com/aliarslan47/RNAForge` — **PRIVATE**, remote `origin` (SSH)
**Referans doküman:** `PLAN.md` **v1.3** (tek referans — Kural 1)
**Son güncelleme:** 2026-08-05

## Şu an nerede kaldık
- **★ m13 AMR + VIRÜLANS OVERLAY TAMAM ve `main`'de (2026-08-05, merge `eeb1cd0`, push). Dalga 2 #A.**
  abricate (yeni `rnaforge-amr` env, 1.4.0) genome'u **CARD** (AMR) + **VFDB** (virülans) tarar → koordinatla
  locus_tag → **DE durumu overlay**. Yeni `rnaforge amr`; gate YOK. **332 test.**
  - **Kod:** `abricate.py` (runner + parser + koordinat örtüşme eşleme + DE overlay); `modules/m13_amr.py`.
    Eşleme gen adıyla DEĞİL koordinatla (abricate adı sembolle eşleşmeyebilir). `config.amr` (AMRConfig:
    amr_db/virulence_db/env/min_id/min_cov). Rapor "Direnç ve Virülans" bölümü + **DB-tarih dürüst notu**
    (abricate bundled DB) + abricate/CARD 2020/VFDB 2019 kaynak.
  - **GSE300731 canlı (yayın kalitesinde konkordans):** 43 AMR (24 DE), 74 virülans (43 DE).
    **ARTAN:** marA (mar regulonu master aktivatör — antibiyotik stres imzası), acrAB/acrD/emrAB (multidrug
    efflux İNDÜKLENDİ), ugd/arn (peptid direnci), enterobaktin siderofor (fep/ent), rcsB (Rcs kapsül).
    **AZALAN:** gadW/gadX (asit direnci). **GO/KEGG/GSEA'nın tümüyle örtüşüyor.** Verdict SUSPECT değişmedi.
  - Spec/Plan: `docs/superpowers/{specs,plans}/2026-08-05-m13-amr*`.
- **★ m12 SEMANTIC REDUCTION (REVIGO) TAMAM ve `main`'de (2026-08-05, merge `7f0cb89`, push). Dalga 1 #3 (son).
  DOWNSTREAM DALGA 1 TAMAMEN BİTTİ (KEGG+GSEA+REVIGO).** Saf Python stdlib (numpy YOK). Gate YOK. **316 test.**
  - **Motor** (`semantic.py`): IC = −log(terim frekansı, arka plan `build_gene2go`'dan) + **Lin** benzerliği
    (obo `_ancestors` yeniden kullanılır, MICA) + **REVIGO-benzeri greedy indirgeme** (namespace-içi, padj
    sıralı; max Lin ≥ eşik → tek temsilci). Eşik `enrichment.revigo_similarity` (0.7).
  - **Kaynaklar:** m09 ORA GO (up/down) + m11 GSEA GO. Çıktı `semantic/reduced_{ora_up,ora_down,gsea_go}.tsv`
    (temsilci + n_collapsed + members). Yeni `rnaforge semantic`. **Figür YOK** (MDS scatter numpy/R → bilinçli sonraya).
  - **Rapor:** "Anlamsal İndirgeme (REVIGO)" bölümü (temsilci tablo + "N→M" özet); Lin 1998 + Supek 2011 kaynak.
  - **GSE300731 canlı:** ora_up 58→24, ora_down 51→24, gsea_go 81→32. Polisakkarit/kolanik asit/slime layer
    ailesi temsilcilerde toplandı; **kapsül/zarf-stres teması korundu**. Verdict SUSPECT değişmedi. Rapor 4.4 MB.
  - Spec/Plan: `docs/superpowers/{specs,plans}/2026-08-05-m12-semantic*`.
- **★ m11 GSEA TAMAM ve `main`'de (2026-08-05, merge `835b651`, push). Dalga 1 #2.**
  Motor **fgsea** (Bioconductor, altın standart — DESeq2 kararıyla tutarlı; `rnaforge-de` env'ine kuruldu,
  `envs/rnaforge-de.yml` güncel). ORA'dan farklı: **tüm genlerin ranked listesi** (DESeq2 `stat` = Wald).
  Gen-seti kurucuları (GO/KEGG) m09/m10'dan yeniden kullanıldı. Yeni `rnaforge gsea`; gate YOK. **298 test.**
  - **Girdi** (`gsea.py`): `write_rnk` (stat, NA atılır) + `invert_to_gmt` (gen→set ters çevir → GMT).
  - **Motor** (`scripts/gsea.R`): `fgsea::fgsea(minSize/maxSize)`; çıktı işaretli **NES** + öncü genler
    (locus_tag→sembol); NES dot-plot (okunur layout, 0-çizgisi). `gsea_min_size`(15)/`gsea_max_size`(500) config.
  - **Rapor:** yeni "Gen Seti Zenginleştirme (GSEA)" bölümü — GO+KEGG için ±NES tablo (term/NES/padj/size/
    **öncü genler**) + figür; Yöntem/Kaynak (Subramanian 2005 + fgsea) yalnız gsea koştuysa. Rapor 4.4 MB,
    14 gömülü figür (8 DE + 2 GO + 2 KEGG + 2 GSEA).
  - **GSE300731 canlı:** GO +34/−47, KEGG +3/−7 anlamlı NES. **Pozitif (artan):** polisakkarit/external
    encapsulating (kapsül), peptidoglikan biyosentezi, ribozom. **Negatif (azalan):** oksidatif fosforilasyon,
    glikoliz, **quorum sensing (NES −2.25)**, amino asit metabolizması. **ORA (GO+KEGG) ile birebir uyumlu**,
    enterololin zarf-stres mekanizmasıyla tutarlı. Verdict SUSPECT değişmedi.
  - Spec/Plan: `docs/superpowers/{specs,plans}/2026-08-05-m11-gsea*`.
- **★ m10 KEGG PATHWAY ORA TAMAM ve `main`'de (2026-08-05, merge `b78c62a`, push). Dalga 1 #1.**
  m09 motorunu (`enrichment.py`) **DEĞİŞTİRMEDEN** kullanır (jenerik gen→set). Yeni `rnaforge kegg`
  subcommand; gate YOK, verdict m06'dan taşınır. **282 test yeşil.**
  - **Annotation** (`kegg_annotation.py`): KEGG REST 3 dosyası (link/pathway, list/pathway, list/<org>)
    → gen→pathway; join **KEGG b-number → gen sembolü → locus_tag** (TAM+BENZERSİZ, m09 `_symbol_to_locus`
    yeniden; belirsiz ATILIR). Global/overview map'ler (01100 vb.) ORA'dan hariç. **Organizma-agnostik**
    (`enrichment.kegg_organism` config'ten; eco/hsa/mmu → **ökaryota taşınır**).
  - **Figür/rapor:** `enrichment.R` + manifest **parametrize** (title/basename prefix) — GO ve KEGG ortak;
    m09 R yolu bozulmadı. Rapordaki "Fonksiyonel Zenginleştirme" bölümü **GO + KEGG alt-bölümleri** gösterir
    (tolerant). Yöntemler'e KEGG paragrafı + Kaynaklar'a Kanehisa & Goto 2000 (yalnız kegg koştuysa).
  - **Referans (gitignore'lu):** `references/kegg/eco/{pathway_links,pathway_names,gene_list}.tsv` (KEGG REST).
  - **GSE300731 canlı:** 6 UP + 9 DOWN anlamlı yolak, 1557 gen KEGG-eşlemeli. **UP = peptidoglikan biyosentezi
    (hücre duvarı stresi) + ekzopolisakkarit + siderofor/enterobaktin** (GO "kolanik asit/slime layer/
    enterobaktin" ile birebir). **DOWN = oksidatif fosforilasyon (padj 1e-12)/glikoliz/nitrojen** (respirasyon —
    GO ile birebir). Rapor 3.6 MB, 12 gömülü figür (8 DE + 2 GO + 2 KEGG). Verdict SUSPECT değişmedi.
  - Spec/Plan: `docs/superpowers/{specs,plans}/2026-08-05-m10-kegg*`.
- **★ m09 GO FONKSİYONEL ZENGİNLEŞTİRME (ORA) TAMAM ve `main`'de (2026-08-05, merge `874068a`, push).**
  A yalın yol. Zincir: m06→m07→**m09**→m08. Yeni `rnaforge enrich` subcommand; **gate YOK**, verdict
  m06/m07'den değişmeden taşınır. **263 test yeşil** (37 yeni).
  - **Annotation birleştirme** (`go_annotation.py`): GFF **otorite** (Ontology_term + go_process/
    function/component → id/namespace/ad) → **GAF güvenli doldurma** (yalnız GFF-GO'suz genlere,
    TAM+BENZERSİZ gen sembolü; belirsiz eşleşme ATILIR; kaynak damgalı GFF|GOA) → **obo propagation**
    (`go-basic.obo` is_a+part_of ata-kapanış, döngü/obsolete korumalı).
  - **ORA** (`enrichment.py`): stdlib **hipergeometrik** (`math.comb`) + **BH FDR** namespace başına;
    UP/DOWN ayrı; arka plan = anotasyonlu test edilen genler; `min_term_size` gürültü filtresi.
    Çıktı `enrichment/enrichment_{up,down}.tsv` + `gene2go.tsv` denetim izi.
  - **Figür** (`scripts/enrichment.R`, rnaforge-de): namespace-facet dot-plot (fold×padj×gen), PNG300+SVG,
    boş-durum panelli. **m08 rapora GO bölümü** (opsiyonel/tolerant, çift dilli, çalıştırılmadıysa dürüst not).
  - **Referans (gitignore'lu):** `references/go/go-basic.obo` (32 MB, 48329 term, 2026-06 sürümü);
    `references/ecoli_bw25113/ecoli.gaf` = **EBI-GOA `18.E_coli_MG1655.goa`** (54437 anot., 3973 sembol;
    K-12 MG1655 aynı gen sembolleri). QuickGO indirmesi kapalıydı → GOA FTP proteome dosyası kullanıldı.
  - **GSE300731 canlı doğrulandı:** 58 UP + 51 DOWN anlamlı GO terimi. **UP = colanic acid/slime layer/
    polysaccharide biyosentezi + membran** (Rcs/kapsül zarf-stres imzası — makaleyle uyumlu), **DOWN =
    respirasyon/enerji metabolizması** (büyüme durması). Anotasyonlu gen **2278→3835** (GAF doldurma katkısı).
    Rapor 3.14 MB, 10 gömülü figür (8 DE + 2 enrichment). Verdict SUSPECT değişmedi (gate yok).
  - Spec: `docs/superpowers/specs/2026-08-05-m09-go-enrichment-design.md` · Plan: `.../plans/2026-08-05-m09-go-enrichment.md`.
  - **Cila (2026-08-05, `714eef6`, push):** enrichment figürü sıkışıktı → düzeltildi (uzun etiket sarma,
    geniş panel 10in, dinamik yükseklik, 0-tabanlı x-ekseni, BP/MF/CC facet tam Türkçe adlı). Rapora
    eklendi: GO tablolarının altına **sütun + kısaltma açıklaması** (çift dilli), **Yöntemler'e GO ORA
    paragrafı**, **Kaynaklar'a GO/GOA/Benjamini–Hochberg** (yalnız enrich koştuysa). 267 test.

### (önceki durak)
- **DEG tablolarına KOŞUL-BAŞI ORTALAMA EKSPRESYON sütunları eklendi ve `main`'de (2026-08-04, merge
  `bc4fba2`, push).** Up/Down tablolarında artık her koşul için ortalama normalize ekspresyon
  (`<control> ort.`, `<enterololin> ort.`) + baseMean. `normalized_counts.tsv` + `coldata.tsv`
  m08 girdisine eklendi (load_report_inputs zorunlu). Örn. gadE: kontrol 28071 → enterololin 12.6.
  226 test yeşil. **NOT: GO/KEGG enrichment kullanıcı isteğiyle SONRAYA bırakıldı** (bkz. bellek
  `reminder_rnaforge_go_enrichment`) — ökaryotla birlikte sıradaki büyük işler.
- **Methods + References güncellendi ve `main`'de (2026-08-04, merge `5b1d545`, push).** Yöntemler bölümü
  artık çift dilli düzgün BİLİMSEL anlatı (DESeq2 makalesi Love ve ark. 2014'ten: medyan-oran norm. →
  empirical-Bayes dispersiyon → negatif binom GLM → Wald → Benjamini–Hochberg; config'ten parametreli).
  References 7 kaynak + **doğrulanmış DOI linkleri** (doi.org). 220 test yeşil.
- **RAPOR ZENGİNLEŞTİRME BİTTİ ve `main`'de (2026-08-04, merge `67bfd6d`, push).** Güncel DE raporlama
  konvansiyonlarına (nf-core/differentialabundance, DEGreport, EnhancedVolcano) hizalandı. **220 test yeşil.**
  - **m06:** `deseq2.R` artık `dispersions.tsv` üretir (gene_id/baseMean/dispGeneEst/dispFit/dispFinal);
    `de_statistics.json`'a `n_up`/`n_down` eklendi (`count_up_down` helper).
  - **m07:** figür sayısı **4→8**, anlatı sırasıyla: PCA · **örnek korelasyon heatmap** · **ekspresyon
    boxplot** · **dispersiyon (plotDispEsts)** · **p-değeri histogramı** · Volcano · MA · Heatmap.
    basename'ler 01–08 yeniden numaralandı; runner'a dispersions arg; korelasyon NaN korumalı (sabit örnek).
  - **m08:** DEG tablosu **ARTAN 25 / AZALAN 25 ayrı tablo**; her figür altında **çift dilli caption**
    (ne gösterir+nasıl okunur); her bölüm başında **intro** (sabit/sayısal, uydurma yorum yok);
    DE bölümünde n_up/n_down gösteriliyor.
  - Code review: Critical yok; 1 Important (n_up/n_down raporda göster) düzeltildi. Gerçek GSE300731'de
    doğrulandı: 8 gömülü figür, n_up 807 + n_down 827 = 1634, verdict SUSPECT (değişmez).
- **★ PROKARYOT MVP TAMAM (2026-08-04) — m08 HTML RAPOR `main`'de (merge `3c96644`, push).** Zincir
  uçtan uca: FASTQ → validate → QC → trim → quant → count → DESeq2 → 4 figür → **tek self-contained
  HTML rapor**. `rnaforge report` subcommand: m06/m07 çıktı sözleşmelerinden (istatistik JSON'ları +
  güvence kartı + figür manifesti + deseq2_results.tsv) her koşuda OTOMATİK, çift dilli (`tr|en`),
  figürler base64 gömülü tek `report/report.html` (~1.4 MB). **213 test yeşil.** Saf Python + stdlib
  (R/jinja2 yok), YENİ kapı YOK, verdict m06/m07'den taşınır.
  - Kod: `rnaforge/report_html.py` (girdi okuyucular + 8 bölüm kurucu + `render_report`),
    `rnaforge/modules/m08_report.py` (`run_report`, m07 ön koşul). 7 task TDD, dal silindi.
  - 9 bölüm: Başlık(run_id) · **Güvence banner(renk-kodlu)** · Dataset/Örnek · Kalite&İşleme ·
    DE Sonuçları · 4 gömülü figür · Top-50 DEG tablo · Methods(config'ten) · References.
  - Code review: Critical/Important YOK; 3 Minor cila uygulandı (başlığa run_id, verdict class
    enum-güvenli, table None koruması). Kalan Minor (defer): DE'de UP/DOWN ayrımı + WARN damgası
    (de_statistics split taşımıyor), figür caption dil-yerelleştirme (başlıklar özel ad).
  - Gerçek GSE300731'de doğrulandı: verdict SUSPECT (m06/m07 ile birebir), 4 figür gömülü,
    gen adları eşlenmiş (pspA/gad*/hde*/ugd), DE özeti "1634/4398".
- **m07 FİGÜRLER BİTTİ ve `main`'de (2026-08-04, merge `c46baf5`, push edildi).** `feat/m07-figures`
  6 task TDD ile tamamlandı, branch silindi. `rnaforge figures` subcommand: m06 DE çıktısından her
  koşuda OTOMATİK 4 statik figür (PCA·Volcano·Heatmap·MA), **PNG 300dpi + SVG**, `runs/.../figures/`
  + `manifest.json`. YENİ veri-kapısı YOK; verdict m06'dan değişmeden taşınır. **191 test yeşil.**
  - Kod: `rnaforge/figures.py` (saf yardımcılar + `run_figures_r`), `rnaforge/scripts/figures.R`
    (ggplot2 + ggrepel, `rnaforge-de` env), `rnaforge/modules/m07_figures.py` (`run_figures`).
  - `rnaforge-de` env'e **r-ggrepel + r-svglite** eklendi (kuruldu + `envs/rnaforge-de.yml`'de).
  - **Code review 2 GERÇEK Critical yakaladı (ikisi de düzeltildi):** heatmap `hclust`, <2 anlamlı DEG
    VEYA sıfır-varyanslı DEG'de çöküp tüm m07'yi exit 1 yapıyordu → "geçerli biyolojide FAIL yok"
    kuralını ihlal. Fix: NaN ölçek satırlarını temizle, <2 satırda kümeleme yapma, 0 satırda boş-durum
    paneli. Regresyon testi eklendi (0-DEG/1-DEG/sıfır-varyans, env-gated). R stdout/stderr artık
    `logs/figures.log`'a yazılıyor.
  - **Gerçek GSE300731 koşusunda doğrulandı:** 4 figür üretildi, biyoloji birebir doğru
    (Volcano top UP=pspA/ugd/ycfJ zarf-stres+kapsül, top DOWN=gadABCE/hdeABD asit-direnç; PCA PC1 %95.8).
  - **Kalan CİLA (Minor, acil değil):** PCA'da sağdaki örnek etiketleri panel kenarından kırpılıyor
    (düz `geom_text`; ggrepel/expand ile çözülebilir); heatmap başlığı TR, diğerleri EN (dil tutarlılığı).
    Manifest'te caption YOK (bilinçli — m08 üretecek). `run_figures(metadata_path)` kullanılmıyor (simetri).
- **GERÇEK YAYIN-VERİSİ DOĞRULAMASI BİTTİ (2026-08-03) — açık konu KAPANDI.** Veri: **GSE300731**
  (Nature Microbiology 2025, "enterololin" dar-spektrum antibiyotik, Brown Lab). *E. coli* BW25113
  ΔbamBΔtolC, **5× MIC enterololin vs kontrol, 4h, 3'er replika** (6 örnek, PRJNA1281986). Referans
  BW25113 GCF_000750555.1. Girdiler gitignore'lu `raw/GSE300731/` + `references/ecoli_bw25113/`.
  Config: prok GFF → `feature_type=CDS`, `attribute=locus_tag`; `de.reference=control`.
- **Uçtan uca canlı koştu** (`runs/20260803_143036_GSE300731/`): hizalama %98.9–99.3, atama %85–86,
  replika kor. min 0.98, **DE 1634 anlamlı/4398 gen**. Tek WARN: ctrl_rep3 ham-QC GC → dürüstçe
  SUSPECT damgalandı (kapı sistemi çalışıyor). Kalan tüm kapılar PASS.
- **KONKORDANS (Katman A/B) — güçlü:** makalenin kendi kallisto abundance'ları (GEO suppl) → aynı
  DESeq2 → bizim ham-FASTQ (bowtie2+featureCounts) LFC'leriyle **Pearson r=0.972 / Spearman 0.958**
  (3592 ortak gen), **makale DEG recall %92.6**, **yön uyumu %99.9**. Biyolojik: top UP = zarf-stres
  imzası (pspA/pspC/spy + Rcs/kapsül rcsA/ugd/gmd/wzc) = makalenin doğruladığı **LolCDE** (lipoprotein
  taşıma) hedefiyle birebir; top DOWN = asit-direnç (gad*/hde*) baskılanması.
- **Analiz/figür script'leri scratchpad'de** (repo dışı): `run_all.sh`, `figures.R`, paper_de.
  Konkordans + figürler tek seferlik prototip — **ASIL İŞ: bunları pipeline'a m07/m08 olarak kur.**
- **YENİ İSTEK (kullanıcı, 2026-08-03):** figürler+rapor pipeline'ın kalıcı parçası olsun, **her koşuda
  OTOMATİK** çıksın; **güncel/modern, güzel, YÜKSEK ÇÖZÜNÜRLÜK** görseller. → m07 (figürler) + m08 (rapor)
  şimdi sıradaki iş; tasarım kalitesi çıtası yüksek. Plotlama: **rnaforge-de'de ggplot2 KURULU** (ggrepel
  YOK — ya kur ya kaçın). matplotlib hiçbir env'de yok. Konkordans figürü STANDART DEĞİL (referans DEG
  tablosu gerektirir; müşteri koşusunda olmaz) → yalnız doğrulama aracı, m07'ye girmez.
- **m07 SPEC + PLAN YAZILDI ve commit'lendi (2026-08-03), dal `feat/m07-figures`.** Kod HENÜZ YOK.
  - Spec: `docs/superpowers/specs/2026-08-03-m07-figures-design.md` · Plan: `docs/superpowers/plans/2026-08-03-m07-figures.md` (6 task, TDD).
  - Onaylı kararlar: figürler **PCA·Volcano·Heatmap·MA** (koşunun KENDİ verisinden); **statik yüksek çöz — PNG 300dpi + SVG**;
    R/ggplot2 `rnaforge-de` env (ggrepel+svglite EKLENECEK); m07 **gate yok/FAIL yok** (m06 gibi); ön koşul m06; çıktı `runs/.../figures/` + `manifest.json` (m08 tüketir).
  - **Konkordans/makale figürü m07'ye GİRMEZ** (Ali netleştirdi) — yalnız tek seferlik doğrulama aracıydı.
  - **DEVAM: Task 1'den başla** (`gene_name_map`, saf TDD) → executing-plans/inline. `conda run -n rnaforge-core --cwd <repo> python -m pytest -q`.
- (Önceki) **m06 DESeq2 `main`'de**, 181 test. PROKARYOT MVP DE zinciri tam.
- **AÇIK KONU (ARTIK KAPALI):** ~~GERÇEK yayımlanmış veri seti ile doğrulama henüz YOK.~~
  Testler + canlı smoke SENTETİK (Kural 8: gerçek/müşteri verisi repo'da yok — doğru). Ama
  Katman A/B doğruluğu için gerçek bakteri veri seti seçilmeli (bkz. "Açık konu" bölümü,
  demo veri kriterleri). SIRADAKİ GERÇEK İŞ olabilir.
- Önceki: m01, m02, m03, m04(prok), m05(prok) `main`'de.
- **Test komutu:** `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q`
  (repo dışından çağırırsan `tests.conftest` importu kırılır — yanlış alarm verir.)
- **m05 detay:** featureCounts TÜM BAM'lere tek çağrı → native matris; sütun→sample_id KONUMLA
  (BAM adıyla değil). Veri kapısı `assignment_rate` (ÜÇÜNCÜ FAIL kapısı). featureCounts params
  config-driven: `quantification.feature_type`(exon)/`attribute`(gene_id) — prok GFF3 için ez
  (CDS/locus_tag). `quantification` KNOWN_TOP_LEVEL_KEYS'e eklendi. Yanlış feature_type → yüksek
  sesle hata, exit 1 (sessiz boş matris YOK). Ön koşul m04; zincir m01→m03→m04→m05.
- **Çalıştırma NOTU:** `python -m rnaforge.cli` ÇALIŞMAZ (main-guard yok); entry point
  `rnaforge` kullan. Referans yolları config'te göreliyse CWD'ye bağlı — smoke'ta cd gerekti.

## Tamamlanan (kalite kapıları planı)
| Task | İş | Durum |
|---|---|---|
| 1 | `gates.py` sözleşmesi (PASS/WARN/FAIL, gates.json) | ✅ re-review temiz |
| 2 | Profiller veri olarak (`profiles/*.yml`) + eşik ezme | ✅ re-review temiz |
| 3 | `subject` sütunu + `looks_paired()` detektörü | ✅ review temiz |
| 4 | `validate_design` → `GateResult` döndürüyor (B kararı) | ✅ re-review temiz |
| 5 | m01 kapıları yazıp zorluyor | ✅ review temiz (mutasyon testli) |
| 6 | Güvence kartı + CLI verdict | ✅ controller + final review (m01 gates.json'ı enforce'tan ÖNCE yazıyor → FAIL'de kart INVALID, doğrulandı) |
| 7 | Dokümantasyon (README/PLAN v1.4) | ⏸ **ERTELENDİ** (m02-m03 sonrası) |

**2026-07-30 tamamlanan (bu oturum):**
- **config sertleştirme:** bilinmeyen üst-seviye anahtar (`design:` → `de.design`, yazım
  hataları) artık `ConfigError` ile reddediliyor (`KNOWN_TOP_LEVEL_KEYS`). Sessiz yutma bitti.
- **test temizliği:** `test_m01_validate.py`'de ikinci `_illumina` gölgeleme kaldırıldı,
  tek kanonik helper (200-read), gates importu tepeye taşındı.
- **Task 6 final review:** temiz. m02'ye geçiş serbest.

Reviewer'ların yakaladığı 3 gerçek bug (hepsi düzeltildi):
1. `gates.json` bozulunca önceki modüllerin kaydı **sessizce siliniyordu** (atomik yazma yok).
2. `~subject + condition` doygun tasarımda **üç kapı da yeşil yanıyordu** → DESeq2 kriptik
   çökme + "TRUSTWORTHY" damgası. (Spec'in deliğiydi, implementer hatası değil.)
3. Güvence kartı **FAIL'de hiç yazılmıyordu** — en çok gerektiği anda yok.

## SIRADAKİ İŞ
1. ~~config.py sessiz hata fix'i~~ ✅ BİTTİ (2026-07-30)
2. ~~Final whole-branch review + `_illumina` gölgeleme~~ ✅ BİTTİ (2026-07-30)
3. ~~`feat/kalite-kapilari` → `main` merge~~ ✅ BİTTİ (2026-07-30)
4. ~~m02 = FastQC~~ ✅ BİTTİ (2026-07-30) — `main`'de (merge `db4b501`).
5. ~~m03 = fastp (nazik trimming)~~ ✅ BİTTİ — `main`'de (`25b3ca2`).
6. ~~m04 = prokaryot hizalama (bowtie2)~~ ✅ BİTTİ — `main`'de (`7595d66`).
7. ~~m05 = prokaryot count matrisi (featureCounts)~~ ✅ BİTTİ — `main`'de (`caba117`).
8. ~~m06 = DESeq2~~ ✅ BİTTİ (2026-07-30) — spec+plan `.../2026-07-30-m06-deseq2*`; env `rnaforge-de`
   (bioconductor-deseq2 1.50.2), betik `rnaforge/scripts/deseq2.R` (dispersiyon fallback'li).
9. ~~`feat/m06-deseq2` → `main` merge + push~~ ✅ BİTTİ.
10. ~~GERÇEK VERİ doğrulaması (Katman A/B, GSE300731)~~ ✅ BİTTİ (2026-08-03) — konkordans r=0.972.
11. ~~m07 figürler (PCA/Volcano/Heatmap/MA, otomatik, PNG300+SVG)~~ ✅ BİTTİ (2026-08-04) —
    `main`'de (merge `c46baf5`). Manifest m08'in tüketeceği sözleşme.
12. ~~m08 = HTML rapor~~ ✅ BİTTİ (2026-08-04) — `main`'de (merge `3c96644`). **PROKARYOT MVP TAMAM.**
13. ~~m09 = GO ENRICHMENT + rapora ekleme~~ ✅ BİTTİ (2026-08-05) — `main`'de (merge `874068a`).
    A yalın yol; GFF+GAF+obo; `rnaforge enrich`; GSE300731 canlı doğrulandı. **PROKARYOT MVP + GO TAMAM.**
14. ~~m10 = KEGG pathway ORA (Dalga 1 #1)~~ ✅ BİTTİ (2026-08-05) — `main`'de (merge `b78c62a`).
    Motor yeniden kullanım; `rnaforge kegg`; GSE300731 canlı (peptidoglikan/respirasyon, GO ile uyumlu).
14b. ~~m11 = GSEA (Dalga 1 #2)~~ ✅ BİTTİ (2026-08-05) — `main`'de (merge `835b651`). fgsea (rnaforge-de),
    işaretli NES; `rnaforge gsea`; GSE300731 canlı (ORA GO+KEGG ile birebir uyumlu).
14c. ~~m12 = Semantic/REVIGO (Dalga 1 #3, son)~~ ✅ BİTTİ (2026-08-05) — `main`'de (merge `7f0cb89`). Lin,
    saf Python; `rnaforge semantic`; 58→24 vb. indirgeme. **DOWNSTREAM DALGA 1 TAMAM.**

### Downstream analiz kuyruğu (Ali seçti 2026-08-05; ökaryot yolundan ÖNCE)
Karar: prokaryot odaklı ama ökaryota taşınabilenler agnostik tasarlanır. WGCNA ELENDİ (6 örnek zayıf).
- **Dalga 1 ✅ TAMAM:** ~~KEGG ORA~~ (m10) → ~~GSEA~~ (m11, fgsea) → ~~Semantic/REVIGO~~ (m12, Lin). Hepsi agnostik.
- **Dalga 2 (bakteri overlay):** ~~#A AMR+Virülans~~ ✅ (m13, abricate CARD/VFDB) → **★ SIRADAKİ #B: Operon
  analizi** (operon tahmini + DEG koordinasyonu) → #C: PPI (STRING) + community detection.
  **DİĞER SEÇENEKLER:** Ökaryot yolu (m04-euk salmon + m05-euk tximport); küçük cila (REVIGO MDS, m07 PCA etiket).
- **Dalga 2 (bakteri overlay, antibiyotik verisine birebir):** AMR (CARD/AMRFinder) + Virulence (VFDB) →
  Operon → PPI (STRING) + community detection.
- **Agnostik-KIRAN (en sona, opsiyonel/E.coli):** Regulon / Sigma factor (RegulonDB/EcoCyc — yalnız E.coli).
- **Düşük değer (bakteride):** Reactome, WikiPathways — uyarlanır ama içi boş.
- **Küçük:** Batch correction (ComBat/removeBatchEffect; DE'de `~batch+condition` zaten mümkün).
15. **SONRAKİ (downstream sonrası):** Ökaryot yolu (m04-euk salmon 2.3.4 + m05-euk tximport) — MVP'nin
    ikinci kolu. salmon 2.3.4 CLI/index ÖNCE doğrulanmalı. Bellek `reminder_rnaforge_eukaryote`.
    m06–m10 zaten agnostik; ayrım YALNIZ m04/m05'te. Ayrıca GO agnostik (eggNOG B yolu).
    - Cila: m07 PCA etiket kırpması; m09 build_gene2go GFF'i iki kez parse (küçük verimlilik).

## Kalite kapıları — Ali ile onaylanan kararlar (2026-07-20)
Gerekçe: *"Yalancı sonuç asla istemem, müşteri güvenceli alsın"* + *"Sorun varsa sorun var
densin, patladıysa bilelim."* Pipeline doğru olsa bile kötü girdi MAKUL GÖRÜNEN sahte sonuç
üretir (%8 hizalama da count matrisi + p-değeri + rapor üretir).
- **Hedef kitle: D (karma)** — bakteriyel derin, altyapı genele açık.
- **İkili politika:** FAIL = sonuç GEÇERSİZ (durur, biyolojik çıktı YOK) · WARN = ŞÜPHELİ
  (üretilir + damgalanır).
- **Eşikler veri:** `profiles/{prokaryote,eukaryote}.yml`, `organism_type` seçer.
  Ökaryot BİLİNÇLİ gevşek + "geniş toleranslı" damgası (ökaryot doğrulaması YOK; uydurma
  eşik kapı sistemini itibarsızlaştırır). Ezilen eşik rapora YAZILIR — sessiz gevşetme yok.
- **FAIL çıktısı:** teşhis raporu (hangi kapı, ölçüm, eşik, sorumlu örnek, ne yapılmalı).
  "Damgalı ama üretilmiş sonuç" REDDEDİLDİ (grafik kopyalanır, damga kaybolur).
- **Eşleşmiş tasarım (öncesi/sonrası):** metadata'ya `subject` sütunu. Tespit VAR, karar YOK —
  eşleşmiş görünüp design kullanmıyorsa DUR ve sor; `paired: false` ile bilerek geçilir.
- **Kapsam:** çerçeve + tasarım kapıları (m01, metadata'dan çıkar) ŞİMDİ; veri kapıları
  kendi modülüyle (m04 yazılırken `alignment_rate` da yazılır).

## Ortam
- 54 GB RAM · 16 çekirdek · R sistemde kurulu (`/usr/bin/Rscript`)
- `rnaforge-core` (python 3.11 + pyyaml + pytest, `pip install -e .`)
- **Araç env'leri KURULDU (2026-07-20):** `rnaforge-qc` (FastQC 0.12.1, fastp 1.3.6) ·
  `rnaforge-quant-prok` (bowtie2 2.5.5, samtools, featureCounts 2.1.1) ·
  `rnaforge-quant-euk` (salmon 2.3.4)
- **DİKKAT:** salmon **2.3.4** geldi, PLAN 1.x varsayıyordu. m04 yazılmadan ÖNCE
  index/CLI davranışı doğrulanmalı — körlemesine güvenme.

## Onaylanan kararlar
1. **Orkestrasyon:** Python paketi, `ali-wgs-pipeline` deseni (config.yaml, runs/, conda env'ler).
2. **Yönlendirme:** `organism_type` (prokaryote|eukaryote) ZORUNLU, varsayılanı yok.
   - prokaryote → bowtie2 + featureCounts · eukaryote → Salmon + tximport/tx2gene
   - Ayrım YALNIZCA `m04`/`m05`'te; ikisi de **gen × örnek count matrisi** sözleşmesinde
     buluşur. `m06`/`m07`/`m08` organizma tipini bilmez.
3. **Platform:** FASTQ'dan tespit (Illumina/ONT/PacBio). MVP **Illumina-only**;
   ONT/PacBio tespit edilir + net hatayla REDDEDİLİR. Kütüphane kimyası
   (stranded/rRNA-polyA) FASTQ'dan tespit EDİLEMEZ → config'ten.
4. **Trimming:** fastp NAZİK — adapter + zorunlu min-uzunluk filtresi;
   agresif kalite trimming KAPALI. Gerekçe: Williams et al. 2016 (BMC Bioinformatics)
   agresif trimming genlerin >%10'unun ekspresyon tahminini bozuyor.
5. **DE motoru:** R/Bioconductor DESeq2 birincil + **pydeseq2 çapraz kontrol**
   (yalnız doğrulamada, her run'da DEĞİL).
6. **Doğrulama iki katmanlı:** Katman A (yayımlanmış count matrisi → DE → birebir yakın
   beklenir) · Katman B (ham FASTQ → tüm pipeline → yayınla KONKORDANS; birebir değil).
7. **Dosyalama:** deney-bazlı `runs/<ts>_<id>/{raw_qc,quantification,...}` (PLAN §14).
   WGS'deki örnek-bazlı `<id>/analiz/` düzeninden BİLİNÇLİ ayrılış — DE/PCA/volcano
   tek örneğe değil deneye ait.
8. **Dil:** kod/log EN · `README.md` (EN) + `README.tr.md` (TR) · PLAN/DURUM/spec TR ·
   HTML rapor dili config'ten (`tr`|`en`).

## MVP kapsamı (PLAN §2.1)
VAR: doğrulama+platform tespiti → FastQC → fastp → quant (2 yol) → count matrisi →
DESeq2 → PCA/Volcano/Heatmap → HTML rapor.
YOK (Faz 2+): ONT, MultiQC, GO/KEGG, GSEA, workflow diyagramı, dashboard, PDF.

## Açık konu
- **Demo veri seti SEÇİLMEDİ.** Kriterler PLAN §3.2'de yazılı: bakteriyel, public
  Illumina FASTQ, üst segment dergi (Frontiers/PLOS One/Sci Rep HARİÇ), 2024-2026
  (tercih 2025-26), GÜÇLÜ sinyal, temiz tasarım, yayında açık DEG tablosu.
- Literatür taraması yapıldı, henüz tüm kriterleri sağlayan net kazanan yok.
  En güçlü aday sinyal açısından: S. coelicolor Scr1 overekspresyonu (1308 DEG) ama
  dergisi (Microbial Cell Factories) üst segment değil.
- MVP kodunu BLOKE ETMEZ (pipeline organizma-agnostik) — Katman B doğrulamasından
  önce seçilmeli.

## Ortam
- 54 GB RAM · 16 çekirdek · 654 GB boş · R sistemde kurulu (`/usr/bin/Rscript`)
- Conda env **`rnaforge-core` KURULU** (python 3.11 + pyyaml + pytest; `pip install -e .` yapıldı).
  Test: `conda run -n rnaforge-core python -m pytest -q`
- Plan B araç env'leri YOK: **fastqc, fastp, salmon, bowtie2 KURULU DEĞİL**
  (featureCounts ve samtools sistemde var).

## Kritik kurallar
- Müşteri verisi/PII asla commit edilmez (`runs/`, `raw/`, `references/` .gitignore'da).
- PLAN.md yeniden yazılmaz, yalnız revize + sürüm yükseltilir (Kural 3).
- Tespit etmek ≠ desteklemek: desteklenmeyen girdi sessizce işlenmez (Kural 7).
