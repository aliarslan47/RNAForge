# DURUM — RNAForge

> Bu dosya "nerede kaldık" anlık görüntüsüdür. Tüm karar detayı Claude belleğindedir
> (`rnaforge-project` memory). Claude bunu anlamlı her durakta ve `/clear` öncesi günceller.

**Konum:** `/home/ali/rnaforge-pipeline/` (git deposu)
**GitHub:** `github.com/aliarslan47/RNAForge` — **PRIVATE**, remote `origin` (SSH)
**Referans doküman:** `PLAN.md` **v1.3** (tek referans — Kural 1)
**Son güncelleme:** 2026-08-18

## Şu an nerede kaldık
- **★★★ ÖKARYOT UZUN-OKUMA YOLU BİYOLOJİK UÇTAN-UCA DOĞRULANDI — `main` `00bcad3`, 515 test (2026-08-18). → DÖRT KOL DA
  BİYOLOJİK DOĞRULANDI.** Veri (2 ajan buluştu, bellek lead'iyle örtüştü): **PRJNA1231053** *S. cerevisiae* Glukoz vs
  Galaktoz (Microorganisms 2025, ONT **PCR-cDNA SQK-PCB114.24**, 4 replika). Referans R64-1-1 cDNA (önden hazırdı) +
  tx2gene başlıklardan. Koşu `runs/20260818_144112_yeast_carbon/`, raw `raw/yeast_carbon/` (4.7GB).
  - **UÇTAN UCA:** validate→NanoPlot→**Pychopper**(survival 0.76-0.87 — SQK-PCB'de SSP/VNP var → Pychopper yolu gerçek
    ONT-kit cDNA'da doğrulandı; Typhi dscDNA'nın aksine)→minimap2→**transkriptom** (mapping %86-96.5)→primer-sayım→tx2gene
    →DESeq2. 6333 gen, **1425 DEG (673↑galaktoz/752↓)**.
  - **★★ GAL REGULONU KONKORDANSI KUSURSUZ:** GAL7 +10.2, GAL10 +8.9, GAL2 +8.8, GAL1 +7.8 (hepsi padj<1e-20, galaktozda
    UP), GCY1 +5.2, GAL80 +3.4; **GAL4 doğru şekilde DEĞİŞMİYOR** (post-translasyonel, aktivatör). Mayanın en ikonik
    galaktoz indüksiyonu ham ONT cDNA'dan üretildi.
  - **★ GERÇEK KOŞUNUN AÇTIĞI BUG + DÜZELTME (`00bcad3`):** `profile_name_for(eukaryote,long)`→`eukaryote_long` profili yoktu
    (prokaryote_long vardı, eukaryote_long yoktu). Eklendi (permissive+damgalı, ONT eşikleri) + regresyon testi.
  - **SIRADA:** cila / temizlik / yeni istek — dört kol tamam. [[reminder_rnaforge_eukaryote]]
- **★★ ÖKARYOT UZUN-OKUMA YOLU (gen düzeyi, transkriptom-hizalama) İMPLEMENTE EDİLDİ — `main` `ec006ca`, 514 test (2026-08-18).**
  Spec `docs/superpowers/specs/2026-08-18-eukaryote-longread-design.md`, plan `.../plans/2026-08-18-eukaryote-longread.md`.
  Kararlar (onaylı): **gen düzeyi** (m06+ reuse), **transkriptoma minimap2** hizalama. m04/m05'te ökaryot dalı içine
  **read_type alt-dallanması** eklendi; m06+ DEĞİŞMEZ; referans mevcut (transcriptome_fasta+tx2gene), genom/GTF yok.
  - **m04** `_quant_euk_long` (YENİ): minimap2 `-ax map-ont`/`map-hifi` (platform'dan) → **transkriptom** (transkript=hedef,
    splice yok); `run_minimap2` reuse; mapping_rate **diagnostik, FAIL kapısı YOK** (tüm long yolları gibi).
  - **m05** `_counts_euk_long` (YENİ): `minimap2.count_primary_alignments` (`samtools view -F 2308` = unmapped/secondary/
    supplementary hariç, okuma başına tek primer) → `tximport.parse_tx2gene` ile **gen'e topla** (gen-içi izoform çoklu-eşleşme
    doğal toplanır) → counts.tsv. Salmon-EM yok (ONT'a uygun). Diagnostik. Gen evreni = gözlenen genler birleşimi (eksik=0).
  - TDD 4 görev; **gerçek samtools sayım testi geçti** (primer-only doğru sayıldı). Regresyon yok (ökaryot-short + prokaryot
    yolları değişmedi).
  - **SIRADA (bilimsel adım):** gerçek ökaryot ONT/PacBio cDNA seti ile uçtan-uca + konkordans (Typhi/airway disiplini,
    kör indirme yok). [[reminder_rnaforge_eukaryote]]
- **★★★ ÖKARYOT YOLU BİYOLOJİK UÇTAN-UCA DOĞRULANDI (gerçek insan verisi) — `main` `9c398b1`, 505 test (2026-08-18).**
  Veri seçimi = bilimsel karar (2 paralel ajan, ENA+literatür buluştu): **PRJNA229998 · insan hava yolu düz kas,
  deksametazon vs kontrol** (Himes 2014, PLoS ONE — **Bioconductor/DESeq2 kanonik seti**). 4 hücre hattı × untreated/Dex,
  paired, `~batch + condition` (hücre hattı=batch, Himes tasarımı). Referans: Ensembl release-110 GRCh38 cDNA (76MB) +
  tx2gene başlıklardan türetildi (207k tx, versiyon eşleşmeli). Decoy YOK (transkriptom-only, bilinen-DEG için kabul).
  Koşu `runs/20260818_110541_airway_dex/` (rapor+8 figür). raw `raw/airway_dex/` (21GB, gitignore).
  - **UÇTAN UCA:** validate→FastQC→fastp→**salmon**(mapping %93.5–94.8)→**tximport**(207k tx→**38366 gen**, toplama çalıştı)
    →DESeq2. **842 DEG (449↑/393↓)**, ref=untreated. profile=eukaryote (permissive, SUSPECT damgalı — 1 WARN dedup, normal).
  - **★★ KONKORDANS KUSURSUZ:** 8/8 glukokortikoid-yanıt işaret geni Dex'te anlamlı **UP** — CRISPLD2 +2.59 (5e-40, Himes
    başlık geni), DUSP1 +2.97, KLF15 +4.49, PER1 +2.98, FKBP5 +4.04, TSC22D3/GILZ +3.20, ZBTB16 +6.30 (1e-179), SPARCL1 +4.52.
  - **★ GERÇEK VERİNİN AÇTIĞI BUG + DÜZELTME (`9c398b1`):** m07 `gene_name_map` `annotation_gff`'e güveniyordu; ökaryotta
    None → `os.fspath(None)` çöküşü. None-güvenli yapıldı (boş gen-map, figürler ENSG ID'yle etiketler) + 2 regresyon testi.
  - **★ GO/KEGG/GSEA/REVIGO DOWNSTREAM ÖKARYOTTA KOŞTU + YOLAK KONKORDANSI (`ac263f0`, spec
    `2026-08-18-eukaryote-enrichment-annotation-design.md`):** anotasyon katmanı GFF-tabanlıydı; insanda GFF-GO yok →
    **transkriptom-sembol köprüsü** eklendi (`parse_transcriptome_symbols`+`parse_annotation_symbols`; gff=None→ENSG→sembol
    FASTA başlıklarından, gene2go boş başlar, GAF/KEGG-by-symbol doldurur; m09/m10/m11/m12 transcriptome_fasta geçer;
    imzalar geriye uyumlu). Referans: `goa_human.gaf` (EBI GOA 828k) + KEGG `hsa`. Config `enrichment: {obo, gaf, kegg_organism: hsa}`.
    **Canlı (airway):** GO 82↑/477↓ (annotated 17568), KEGG 3↑/14↓ (hsa, 7824), GSEA go+136/-12, REVIGO indirgeme.
    **Dex glukokortikoid imzası:** KEGG DOWN=Cytokine-cytokine receptor(5.6e-5)/TNF/chemokine/rheumatoid arthritis (anti-inflamatuar);
    GO UP=response to hormone. 510 test.
  - **SIRADA:** ökaryot UZUN-okuma (minimap2 splice-aware + izoform, yeni alt-sistem) VEYA temizlik. [[reminder_rnaforge_eukaryote]]
- **★★★ ÖKARYOT KISA-OKUMA YOLU (Salmon + tximport) İMPLEMENTE EDİLDİ — `main` `b8f554a`, 503 test (2026-08-18).**
  MVP'nin ikinci ana kolu. Brainstorm→spec (`docs/superpowers/specs/2026-08-18-eukaryote-shortread-design.md`)→
  plan (`docs/superpowers/plans/2026-08-18-eukaryote-shortread.md`)→TDD 4 görev→main→push. Kararlar: kısa-okuma
  önce (uzun sonra), Salmon **decoy-aware** (genom opsiyonel), tximport **`countsFromAbundance="lengthScaledTPM"`**
  (uzunluk-düzeltilmiş sayım → m06 DESeq2 offset gerektirmez, **m06+ agnostik korunur, DEĞİŞMEDİ**).
  - **m04-euk** (`rnaforge/salmon.py` + m04 router dalı): `organism_type==eukaryote` → Salmon; eski
    `NotImplementedError` KALDIRILDI. Index: genome_fasta varsa gentrome+decoys.txt decoy-aware, yoksa
    transkriptom-only + yüksek sesle not. `salmon quant -l A` (paired destekli). mapping_rate → alignment_rate
    FAIL kapısı (eukaryote.yml permissive; `_MappingAdapter` ile `build_alignment_gates` yeniden kullanılır).
    Env `rnaforge-quant-euk` (salmon 2.3.4 kuruluydu). **Gerçek salmon entegrasyon testi GEÇTİ** (index+quant).
  - **m05-euk** (`rnaforge/tximport.py` + `scripts/tximport.R` + m05 router dalı): quant.sf×örnek + tx2gene →
    gen matrisi; `counts.tsv` (`gene\t<sid>`) ortak sözleşme. `bioconductor-tximport` `rnaforge-de` env'ine
    kuruldu (`envs/rnaforge-de.yml`; ikinci R env yok). assignment FAIL kapısı yok (salmon zaten atadı, diagnostik).
  - **config:** `Reference.genome_fasta` ökaryotta opsiyonel (alanlar zaten vardı); `REQUIRED_REFERENCE[eukaryote]`
    = (transcriptome_fasta, tx2gene) DEĞİŞMEDİ.
  - **SIRADA (bilimsel adım, kör indirme yok):** Typhi disiplini — yayınlanmış DEG'li gerçek ökaryot Illumina seti
    seç (insan/fare/maya/Arabidopsis; ≥3 replika; GENCODE/Ensembl transkriptom + tx2gene) → uçtan uca + konkordans.
    ([[reminder_rnaforge_eukaryote]])

- **★★★ UZUN-OKUMA KOLU BİYOLOJİK UÇTAN-UCA DOĞRULANDI (gerçek bakteri ONT cDNA) — `main` `44c2e06`, 492 test (2026-08-18).**
  Uzun-okuma kolu şimdiye dek yalnız teknik doğrulanmıştı; ilk kez **yayınlanmış bir çalışmayla biyolojik konkordans**.
  **Veri seçimi = bilimsel karar** (kör indirme yok): 2 paralel araştırma ajanı ENA+literatür taradı, ikisi de aynı sette
  buluştu → **PRJNA1254696 · *S.* Typhi + rifampin** (Lee & Song, **eBioMedicine 2025**; 2 koşul × 3 gerçek replika, ONT
  cDNA, GEO GSE295448). Elenen adaylar: NAR2025 ısı-stresi (replikalar birleştirilmiş n=1), microbepore (untreated-vs-TEX =
  teknik), Nano3P (havuzlu/protokol). Referans **S. Typhi CT18** (GCF_000195995.1; Vi kapsül lokusu tvi/vex GO-anotasyonlu).
  - **★ GERÇEK VERİNİN AÇTIĞI BOŞLUK + DÜZELTME:** bu kütüphane **rastgele-primer'lı dscDNA + native barcoding** (ONT-kit cDNA
    DEĞİL) → SSP/VNP strand-switch primeri yok → **Pychopper 20k okumadan 429'unu tuttu (%2)**. Yeni **`library.full_length_cdna`**
    bayrağı (varsayılan `true`, geriye uyumlu): `false` → cdna yolu **Pychopper'ı atlar, chopper-only** (yüksek sesle loglar).
    TDD: config `_as_bool` + Library alanı + m03 dispatch + 4 test. **Gerçek veride survival 0.02 → 0.9974.** Commit `44c2e06`.
  - **UÇTAN UCA (izole env'ler):** validate→NanoPlot→chopper→minimap2(map-ont)→featureCounts-L→DESeq2, hepsi TRUSTWORTHY,
    profile=prokaryote_long. 5032 gen × 6 örnek; **174 anlamlı DEG (39 UP / 135 DOWN)**, ref=dmso.
  - **★★ KONKORDANS ~KUSURSUZ (makale ana bulgusu):** makale "rifampin Vi biyosentez genlerini **ve** SPI-1 T3SS'i azalttı" +
    RT-qPCR tviA/tviB. Bizde: **Vi kapsül tvi/vex 10/10 DOWN** (tviA −4.58 padj 9e-16; vexA −4.09 5e-11; 8/10 padj<0.05),
    **SPI-1 T3SS 16/16 anlamlı DOWN** (invA −4.49 **2e-45**; hilA −4.87; sopE −5.66). Yön+anlamlılık birebir. RNAP-inhibitörü
    baskın DOWN (yüksek-AT ada baskılanması, makaleyle tutarlı). Koşu: `runs/20260818_085055_typhi_rif/` (figürler+rapor üretildi).
  - **★ GO+KEGG DOWNSTREAM DE TYPHI'DE KOŞTU + YOLAK-DÜZEYİ KONKORDANS (`2e9db51`):** KEGG `sty` referansı çekildi
    (`references/kegg/sty/`, gitignore); config'e `enrichment.obo`+`kegg_organism: sty` eklendi. GO 35↑/62↓ (GFF Ontology_term
    otorite, GAF yok), KEGG 1↑/3↓, GSEA (go +48/-14), REVIGO indirgeme. **GO DOWN**: `capsular polysaccharide transport`
    (Vi kapsül). **KEGG DOWN**: `Bacterial invasion of epithelial cells` (padj 2e-6, sip/sop), `Salmonella infection` (2e-5),
    `Flagellar assembly`, + nucleotide-sugars→tviB/tviC. Makalenin iki ana gen seti (Vi + SPI-1) yolak düzeyinde de DOWN.
    Rapor `--force` yenilendi. AMR/operon/PPI KOŞULMADI (STRING taxid/abricate prep gerekir; opsiyonel, E.coli'de doğrulanmış).
  - **SIRADA:** ökaryot yolu VEYA PacBio HiFi canlı VEYA temizlik (Nano3P 26GB gereksiz). 
- **★★ m00 BASECALL (ham sinyal FAST5/POD5 → FASTQ) TAMAM ve `main`'de (2026-08-17, merge `e60d981`, 485 test).**
  Ali: "FAST5/POD5 gelirse pipeline görsün, çözümlesin, başlatsın." **Fizibilite: GPU VAR (RTX 4050, 6GB, driver
  566.14) → dorado GPU basecalling uygulanabilir** (CPU olmazdı). Plan `docs/superpowers/plans/2026-08-17-m00-basecall.md`.
  Kuruldu: dorado 2.1.1 (`/home/ali/tools/dorado-2.1.1-linux-x64/`) + `rnaforge-basecall` env (pod5 0.3.44,
  `envs/rnaforge-basecall.yml`). **Canlı spike GEÇTİ:** gerçek r10.4.1 POD5 → dorado hac (model otomatik indi,
  cuda:0 RTX 4050) → FASTQ. Kod: `rnaforge/basecall.py` (`is_signal_input`/`convert_fast5_to_pod5`/`run_dorado`/
  `basecalled_metadata_path`) + `modules/m00_basecall.py` (per-sample: FAST5→POD5→dorado→FASTQ, FASTQ passthrough,
  **resolved metadata** yazar, diagnostik/kapı yok) + `config.basecall` (dorado_bin/model=hac/device=cuda:all/env/
  models_dir) + CLI `basecall` (ilk aşama) + **m01 resolved metadata'yı otomatik devralır** (handoff sözleşmesi).
  **Canlı e2e:** POD5 4 örnek → `rnaforge basecall` (GPU) → FASTQ → `rnaforge validate` resolved metadata'yı otomatik
  kullandı → TRUSTWORTHY. **Bug (canlı e2e yakaladı+düzeltildi):** resolved metadata göreli yol yazınca load_metadata
  ikiliyordu → mutlak yol yazılır. Böylece ham-sinyal ONT setleri (FASTQ olmasa bile) artık kullanılabilir.
- **★★ UZUN-OKUMA (ONT/PacBio) YOLU — ADIM 1 (tespit→yönlendirme + `library.chemistry`) TAMAM ve `main`'de
  (2026-08-06, merge `a2acbca`, push).** Plan: `docs/superpowers/plans/2026-08-06-longread-step1-routing.md`.
  5 görev TDD, **423 test yeşil** (412→423). Yapılanlar: `platform.py` `read_type_for()` (illumina→short,
  ont/pacbio_hifi→long); long okumalar artık **REDDEDİLMİYOR, YÖNLENDİRİLİYOR** (`SUPPORTED_PLATFORMS`
  genişledi, yalnız `unknown` reddedilir); `library.chemistry` (cdna|direct_rna) config; **m01** read_type+
  chemistry'yi `raw_statistics.json`'a yazar + ONT-long için chemistry ZORUNLU; yeni `rnaforge/routing.py`
  (`resolve_read_type` + `require_short_read`) → **m02–m05'te muhafız**: long okuma gelince YÜKSEK SESLE durur
  (NotImplementedError, yanlış araçla sessiz koşmaz). Canlı smoke: GSE300731 short etkilenmedi; sentetik ONT
  yönlendi→m02 dürüstçe durdu.
- **★ UZUN-OKUMA ADIM 2 (m02-long NanoPlot QC) TAMAM ve `main`'de (2026-08-06, merge `02e30db`, push).**
  Plan `docs/superpowers/plans/2026-08-06-longread-step2-m02-nanoplot.md`, 4 görev TDD, **427 test**.
  `rnaforge/nanoplot.py`: `parse_nanostats` (NanoPlot `--tsv_stats` NanoStats.txt) + `run_nanoplot`
  (`rnaforge-longread` env, `--tsv_stats --no_static`). **m02 read_type'a göre dispatch** (m04 router deseni):
  short→FastQC (`_qc_short`, aynen), long→NanoPlot (`_qc_long`); ikisi de DIAGNOSTIK (kapı yok, durmaz). CLI qc
  dala göre araç adı yazar. `envs/rnaforge-longread.yml` eklendi.
  **★ HAZIRLIK (2026-08-06):** `rnaforge-longread` env KURULDU (minimap2 2.31 · NanoPlot 1.47.1 · Pychopper ·
  chopper 0.13 · samtools 1.24). **microbepore verisi SEÇİLDİ + İNDİRİLDİ:** `raw/microbepore_mg1655/` 10 run
  (PRJNA731531, ENA'dan FASTQ, 5.1 GB, ONT/long doğrulandı). **MG1655 referansı** `references/ecoli_mg1655/`
  (genome+GFF, GCF_000005845.2 — GO/KEGG(eco)/GAF ile birebir tutarlı). **Canlı smoke** (`runs/*_mbp_smoke`,
  4 PCR-cDNA örnek subset): validate→read_type=long, qc→NanoPlot koştu (N50 472-633, meanQ 8.1-8.6), m03 hâlâ
  dürüstçe durdu. microbepore çoğunlukla tek-koşul → DE-sinyal smoke için değil (aday B: glucose-vs-pyruvate sonra).
- **★ UZUN-OKUMA ADIM 3 (m03-long Pychopper+chopper) TAMAM ve `main`'de (2026-08-06, merge `a5aed57`, push).**
  Plan `docs/superpowers/plans/2026-08-06-longread-step3-m03-pychopper.md`, 4 görev TDD, **435 test**.
  `rnaforge/chopper.py` (`run_chopper`, ONT uzunluk/kalite filtresi, stdin→stdout) + `rnaforge/pychopper.py`
  (`run_pychopper` tam-boy cDNA yönlendir/kes + `parse_pychopper_stats`). **m03 read_type dispatch:** short→fastp
  (`_trim_short`, survival FAIL kapısı korunur), long→`_trim_long` (kimyaya göre: **cdna→Pychopper+chopper**,
  **direct_rna→yalnız chopper**); long DIAGNOSTIK (FAIL kapısı yok — long profil Step 6). Trimlenmiş çıktı m04
  sözleşme yolunda. **ÖNEMLİ:** pychopper 2.7.10 pandas-3'te sonda PDF rapor çizerken çöküyor (`_plot_stats`
  `float(Series)`) AMA çekirdek çıktı+stats tam → `run_pychopper` bu **bilinen çöküşü** (çıktı var + stderr'de
  `_plot_stats`) tolere eder (yüksek sesle uyarır), başka her hatada fırlatır. Canlı smoke (mbp_smoke cdna):
  survival 0.54-0.63, m04 hâlâ durdu.
- **★ UZUN-OKUMA ADIM 4 (m04-long minimap2 hizalama) TAMAM ve `main`'de (2026-08-17, merge `931c539`).**
  Plan `~/.claude/plans/rosy-sleeping-lamport.md`, 4 görev TDD, **449 test** (435→449). Durable:
  (1) `rnaforge/minimap2.py` = bowtie2.py deseni: `run_minimap2` (`minimap2 -ax <preset> -t <n> genome fastq
  -o sam` → samtools sort/index; SAM dosyaya, belleğe değil) + `parse_flagstat_mapped` (hizalama oranı =
  `samtools flagstat` **primary-mapped/primary**, minimap2 "overall alignment rate" yazmaz; yazdırılan %'e
  değil SAYIMLARA bakar) + `minimap2_preset` (ont→map-ont, pacbio_hifi→map-hifi). (2) `routing.resolve_platform`
  (raw_statistics.json'dan platform; resolve_read_type deseni) — preset seçimi için. (3) **m04 read_type
  dispatch** (m02/m03 deseni): short→`_quant_short` (bowtie2 gövdesi aynen, alignment_rate FAIL kapısı korunur),
  long→`_quant_long` (minimap2, preset platformdan). Step-1'in m04'teki `require_short_read` muhafızı dispatch'le
  DEĞİŞTİ (m05'te kaldı). (4) **Long DIAGNOSTIK — FAIL kapısı YOK** (m03-long deseni): alignment_rate yalnız
  istatistik; Illumina 0.70 eşiği ONT'yi yanlış FAIL'lerdi (long profil Step 6). Çıktı BAM `quantification/<sid>/
  aligned.sorted.bam` (m05 sözleşmesi). **Canlı smoke** (mbp_smoke, 4 PCR-cDNA microbepore): minimap2 MG1655'e
  hizaladı, oranlar **0.71–0.81** (ctrl1=0.7122 → 0.70 eşiğine sınırda = diagnostik kararı doğruladı); BAM'ler
  üretildi; `rnaforge counts` hâlâ dürüstçe durdu (m05 muhafızı). **NOT:** microbepore trimmed adları `_1`'siz
  (m03 zamanı); smoke metadata `_1`'siz symlink'lerle eşlendi.
- **★ UZUN-OKUMA ADIM 5 (m05-long featureCounts `-L`) TAMAM ve `main`'de (2026-08-17, 453 test).** Plan
  `docs/superpowers/plans/2026-08-17-longread-step5-m05-featurecounts.md`, 3 görev TDD. Durable: (1)
  `run_featurecounts(long_read=True)` → `-L` (ONT/PacBio tek-molekül; `-L`/`-p` bağdaşmaz, paired yok sayılır);
  aynı binary (subread 2.1.1, `rnaforge-quant-prok`). (2) **m05 read_type dispatch** (m04 deseni): short→
  `_counts_short` (featureCounts, assignment_rate FAIL kapısı korunur), long→`_counts_long` (`-L`, paired=False,
  **DIAGNOSTIK — FAIL kapısı YOK**, long profil Step 6). counts.tsv + tpm/fpkm yazımı ortak `_write_count_outputs`
  yardımcısında. (3) Step-1'in m05'teki `require_short_read` muhafızı KALDIRILDI — bu **uzun-okuma kolunun SON
  Step-1 muhafızıydı**; long run artık uçtan uca count matrisine ulaşıyor. **Canlı smoke (mbp_smoke):**
  featureCounts -L → **4308 gen × 4 örnek**, atama %6–16 (düşük = diagnostik kararını doğruladı; 0.50 eşiği
  yanlış FAIL'lerdi), m05 gate YOK; **`rnaforge de` (m06 DESeq2) uzun-okuma matrisi üzerinde ÇALIŞTI** →
  "0 significant / 4308 genes" (microbepore tek-koşul, ctrl/trt yapay → 0 beklenen). **Uzun-okuma kolu artık
  m06+ ile birebir buluşuyor (kod değişmeden).**
- **★ UZUN-OKUMA ADIM 6 (long-read kalite profili + kapıları) TAMAM ve `main`'de (2026-08-17, 462 test).** Plan
  `docs/superpowers/plans/2026-08-17-longread-step6-profile-gates.md`, 4 görev TDD. Durable: (1) Yeni
  `rnaforge/profiles/prokaryote_long.yml` (`permissive: true`, DAMGALI) — ONT bilinçli permissive eşikler:
  `alignment_rate=0.50` (katastrofik/yanlış-referans → **FAIL**), `survival_rate=0.20` + `assignment_rate=0.05`
  (ONT'de doğal düşük → **WARN**, geçersiz kılmaz). Gerekçe belleğe: *"uydurma eşik kapı sistemini itibarsızlaştırır"*
  → temsili ONT veri seti gelene dek permissive+damgalı (ökaryot profili deseni). (2) `quality.profile_name_for(
  organism_type, read_type)` (long→`<organism>_long`) — read_type→profil eşlemesinin tek kaynağı. (3) `build_trim_gates`
  float-tabanlı yapıldı + `build_trim_gates`/`build_count_gates`'e `warn_only` (eşik altı WARN). (4) Long dalları
  kapı yazar: **m04-long alignment FAIL** (`build_alignment_gates`+raise_if_failed), **m03-long survival WARN**,
  **m05-long assignment WARN**; m02-long diagnostik kalır (short m02 gibi). (5) CLI `_load_run_profile(config,run_dir)`
  — güven kartı long run'da `prokaryote_long` (permissive) damgalar (m01 öncesi kısa'ya düşer); TÜM stage kart
  yazımlarına uygulandı (m06+ kart damgasını short'a geri çevirmesin). **Canlı smoke (mbp_smoke, quant+counts
  --force):** kart **profile=prokaryote_long permissive=True**; **alignment PASS** (0.7122>0.50), **assignment PASS**
  (0.0635>0.05, kıl payı = floor doğru), verdict SUSPECT (m06 replicate_correlation WARN, yapay ctrl/trt → dürüst).
  Yalancı FAIL yok — Illumina eşikleri reddederdi, long profil geçirdi + damgaladı.
- **★ UZUN-OKUMA ADIM 7 (rapor read_type rozeti + long-read notları) TAMAM ve `main`'de (2026-08-17, 469 test).
  → UZUN-OKUMA (ONT/PacBio) KOLU KOMPLE (Adım 1-7).** Plan `docs/superpowers/plans/2026-08-17-longread-step7-report.md`.
  Hepsi `report_html.py`: (1) `section_dataset`'e **read_type rozeti** (Okuma tipi: uzun/long, platform yanında;
  kısa→kısa/short, uzun→uzun/long glossu). (2) `_SOFTWARE` platform araçlarına short/long `cond` + yeni long
  girdileri (NanoPlot 1.47.1/Pychopper 2.7.10/chopper 0.13.0/minimap2 2.31); featureCounts "-L" notlu, paylaşımlı.
  (3) `_METHODS_TEXT_LONG` (tr/en): NanoPlot→Pychopper+chopper→minimap2→featureCounts -L anlatısı + permissive
  prokaryote_long notu; `section_methods(read_type=)`. (4) Referanslar read_type'a bölündü: base + `_REFERENCES_SHORT`
  (FastQC/fastp/Bowtie2/Williams) + `_REFERENCES_LONG` (NanoPack2 btad311, minimap2 bty191, Pychopper repo) —
  **kullanılmayan aracı atıflamaz (dürüstlük)**. render_report read_type'ı software flags + methods + refs'e geçirir.
  **Canlı uçtan-uca (mbp_smoke ONT): figures→report OK**, rapor (582 KB) profile=prokaryote_long damgalı, read_type
  rozeti + minimap2/NanoPlot/Pychopper + featureCounts -L + NanoPack2/minimap2 atıfları; **FastQC/Bowtie2 HİÇ geçmiyor**.
- **★ CİLA TURU TAMAM ve `main`'de (2026-08-17, 472 test).** Kalan küçük maddelerden 3'ü otonom kapatıldı:
  (#8) **CLI `--help` çöküşü** — `cli.py:183` seqqc help'inde bare `%` (`rRNA%`) argparse'ı patlatıyordu →
  `%%` escape (`20a8c86`); regresyon testi `tests/test_cli_help.py`. (#10) **m09 çift GFF-parse** —
  `build_gene2go` artık `gene_symbol` de döndürüyor, m09'daki ikinci `parse_gff_go` kaldırıldı (`9d9b9d1`).
  (#9) **README (EN+TR) + PLAN v1.4** — uzun-okuma kolu + downstream (m09-15) + QC (m16-18) + read_type
  yönlendirme + `library.chemistry` + `prokaryote_long` profili belgelendi; PLAN "ONT reddedilir" →
  "yönlendirilir" güncellendi (`8c7b272`). **KALAN (Ali kararı gerekir, kör-otonom değil):** #1 DE-sinyalli
  bakteri ONT veri seti SEÇİMİ (kayıtlı "aday B" YANLIŞTI = maya/ökaryot; bakteri seti hâlâ seçilmedi),
  #3-4 ökaryot kolu (kendi spec'i gerekir), #5-7 downstream ekleri (regulon/batch/pydeseq2).
- **★★ UZUN-OKUMA KOLU TAMAM ÖZET:** validate→route (chemistry) · m02 NanoPlot · m03 Pychopper+chopper · m04 minimap2 ·
  m05 featureCounts -L · long profil+kapılar (align FAIL / survival+assign WARN, permissive+damgalı) · rapor read_type
  farkında. Ortak count matrisinde m06+ ile **kod değişmeden** buluşuyor. **SIRADAKİ SEÇENEKLER (Ali seçecek):**
  (a) DE-sinyalli gerçek ONT veri seti (aday B *E. coli* glucose-vs-pyruvate) ile biyolojik-anlamlı uçtan-uca doğrulama.
  (b) Ökaryot yolu (m04-euk salmon + m05-euk tximport). (c) PacBio HiFi canlı doğrulama (`map-hifi` bağlı, veri yok).
  (d) README/PLAN v1.4 dokümantasyon. Tasarım `docs/superpowers/specs/2026-08-05-longread-arm-design.md`.
- **★ QC TAMAMLAMA (5 düşük-öncelik eksik) TAMAM ve `main`'de (2026-08-05, merge `4a9bd88`, push).** Ali
  "hepsini sırayla ekle, otonom" dedi. 5'i de eklendi, hepsi **diagnostik** (verdict'i asla FAIL ile bozmaz):
  - **F1 per-base baz kompozisyonu + duplikasyon** → m02: `fastqc.py` `parse_per_base_content` +
    `parse_deduplication` (mevcut FastQC zip'inden, yeni araç yok); **dedup_fraction WARN kapısı** (profilde
    prok 0.20/euk 0.15, asla FAIL); per-base kompozisyon figürü. Canlı: dedup %55-62 → PASS.
  - **F2 insert-size** + **F3 read-distribution** + **F4 coverage** → yeni **`m17_alignqc`** (ortak BAM döngüsü,
    kapı yok). `alignqc.py`: samtools stats (insert-size+IS histogram), samtools coverage (kontig derinliği),
    RSeQC read_distribution (CDS/UTR/Intron/Intergenic %). insert-size figürü (paired-only, SE'de atla) +
    coverage figürü. Yeni `rnaforge alignqc`. Canlı: insert **269.6bp**, genom derinlik **324×**, read-dist
    **CDS %91.5 / Intergenic %8.5** (bakteriyel mRNA doğru).
  - **F5 MultiQC toplu görünüm** → yeni **`m18_multiqc`** (kapstone, en son koşar). `multiqc.py` + run dizinini
    tarar → toplu HTML; rapora göreli link. **multiqc 1.35 rnaforge-seqqc'e pip** (env yml güncel). Canlı:
    **6 modül / 20 örnek** topladı.
  - **Figür altyapısı:** bağımsız `qc_plot.py` (matplotlib, rnaforge-seqqc) + `qcplots.py` sarmalayıcı;
    lines/bars tipleri. Figürler best-effort ama **sessiz değil** (hata log+stats'e yazılır).
  - **Rapor:** kalite bölümü genişledi — benzersiz-% sütunu, read-distribution tablosu, 3 QC figürü,
    MultiQC linki (çift dilli). Spec `docs/.../2026-08-05-qc-completion-design.md`. **412 test.**
  - **GSE300731_final canlı:** rapor 6.35 MB, F1-F5 hepsi gömülü; verdict SUSPECT değişmedi (PASS 13/WARN 1;
    dedup PASS eklendi, hiçbir tanısal verdict'i bozmadı). insert 269.6bp · derinlik 324× · CDS %91.5.
- **★ EKSİK KAPATMA (TPM/FPKM + Software/DB tabloları + genel kaynaklar) TAMAM ve `main`'de (2026-08-05,
  merge `16576a3`, push).** Ali 5 referans PDF verdi (r1 ticari şablon + 4 metodoloji makalesi); RNAForge'u
  bunlarla karşılaştırdım (çekirdeği tam karşılıyor, downstream'de ötesinde; "fazla" yok). Kapatılan eksikler:
  - **TPM/FPKM ekspresyon değerleri:** `featurecounts.py` `tpm_fpkm(gene_ids,counts,lengths)` (Length raw
    featureCounts'tan; lengths artık `FeatureCountsResult`'ta). m05 → `quantification/{tpm,fpkm}.tsv`. Rapor
    DE bölümüne "Ekspresyon Düzeyi" notu. (TPM sütun toplamı 1e6 doğrulandı.)
  - **Software + Database sürüm tabloları:** yeni rapor bölümü "Yazılım ve Veritabanları" (koşan araçlar+
    sürümler: FastQC/fastp/Bowtie2/featureCounts/DESeq2/fgsea/abricate/AMRFinderPlus/SortMeRNA/RSeQC/networkx
    + kullanılan DB'ler CARD/VFDB/STRING/KEGG/GO). N_SECTIONS=17, 21 tablo.
  - **Genel kaynaklar:** Dawadi 2025 (Front Genet), Deshpande 2023 (Front Genet), Pola-Sánchez 2024 (Curr
    Protoc), Claussen/EICC 2023 — koşulsuz References'a eklendi. **391 test.** GSE300731_final güncel.
  - **Kalan düşük-öncelik eksik:** insert-size dağılımı, MultiQC toplu, per-base baz komp./duplikasyon,
    hizalama görselleştirme/coverage, RSeQC read-distribution. Spec `docs/.../2026-08-05-m13b-*` (ilgisiz).
- **★ m16 SEKANS QC (rRNA% + STRANDEDNESS) TAMAM ve `main`'de (2026-08-05, merge `05b04ea`, push).** Ali'nin
  "QC ekle" isteği — eksik girdi-QC kapıları. İki bağımsız araç (Ali seçti): **SortMeRNA** (rRNA%) + **RSeQC**
  (strandedness), yeni `rnaforge-seqqc` env. **İki yeni WARN kapısı** (asla FAIL): `rrna_fraction` (>eşik,
  profilde 0.20) + `strandedness_match` (çıkarım≠beyan). Güvence kartına akar. Yeni `rnaforge seqqc`. **385 test.**
  - **Kod:** `seqqc.py` (rRNA-ref **genomdan çıkarılır** — indirme yok, agnostik; subsample gzip/düz oto;
    sortmerna/infer_experiment runner+parser; GFF→BED12). `m16_seqqc.py` ön koşul m04. Rapor "Kalite ve İşleme"ye
    rRNA%/strandedness/uyum satırları.
  - **GSE300731_final canlı:** rRNA **%0.6** (depletion mükemmel → PASS), strandedness **unstranded = beyan**
    (PASS). İki kapı da PASS (PASS 10→12), verdict SUSPECT değişmedi. Yalancı alarm yok, gerçek sinyal doğru.
  - **DÜZELTME:** m03 trimlenmiş çıktı gzipli DEĞİL → `subsample_fastq` gzip/düz oto-algılama (`05b04ea`).
  - Spec `docs/.../2026-08-05-m16-seqqc-design.md`.
- **★ RAPOR YENİDEN YAPILANDIRMA TAMAM ve `main`'de (2026-08-05, merge `9b3bb27`, push).** Ali'nin standart
  RNA-seq bileşen listesine göre: (1) **Bölüm sırası** = liste sırası (Kalite→DESeq2→Figürler→Top DEG→**GO**→
  **KEGG**→GSEA→REVIGO→[AMR/Operon/PPI ek]→Yöntem/Kaynak). (2) **GO ve KEGG artık AYRI üst-bölümler**
  (`section_go`+`section_kegg`; eski birleşik `section_enrichment` bölündü). (3) **Figürler numaralı**
  (Şekil N / Figure N, belge sırasına göre post-pass regex). (4) **Her analiz bir sayfada** — bölüm-başı
  `page-break-before` (print/PDF). Operon/PPI/Güvence Kartı **KORUNDU** (Ali "dur" dedi, silme iptal). N_SECTIONS=16.
  - **★ SIFIRDAN TAM KOŞU doğrulandı** (`runs/20260805_160103_GSE300731_final`, ham FASTQ→rapor, 15 adım):
    DE **1634 anlamlı (öncekiyle birebir → tekrarlanabilir)**, rapor 6.2 MB, 15 bölüm doğru sırada, 19 numaralı
    figür, verdict SUSPECT (ctrl_rep3 GC WARN) değişmedi. **370 test yeşil.**
- **★ m13b AMRFinderPlus (İKİNCİ AMR ARACI, YAN-YANA) TAMAM ve `main`'de (2026-08-05, merge `b04480b`, push).**
  Reviewer geri bildirimi (bağımsız DB-eşleşmeli AMR modülü) karşılandı — aslında m13 zaten CARD/VFDB'yi bağımsız
  modül olarak veriyordu; şimdi **AMRFinderPlus** ikinci araç olarak eklendi ve AMR tablosunda CARD ile **yan yana**
  (konkordans). `abricate.py`: `run_amrfinder`+`parse_amrfinder` (Type AMR/STRESS, abricate ile aynı dict şekli →
  map/overlay yeniden kullanılır). `m13_amr.py`: CARD ∪ AMRFinderPlus locus_tag'te birleşir; tablo sütunları
  `gene·CARD·AMRFinderPlus·%id·DE`. `config.amr.amrfinder_organism`(verilmezse eski davranış)/`amrfinder_env`
  (**ali-amrfinder** mevcut, DB 2026-05-15). Rapor Feldgarden 2021 kaynak. **370 test.** GSE300731 canlı:
  **CARD 43, AMRFinderPlus 5, her ikisi 4** (efflux mdtM/acrF/emrE + ampC β-laktamaz; ariR yalnız AMRFinderPlus;
  ~39 yalnız CARD = küratörlü vs geniş). Verdict değişmedi. Spec `docs/.../2026-08-05-m13b-amrfinder-design.md`.
- **★ README (EN+TR) TAM GÜNCELLENDİ (2026-08-05, `5446a72`).** m01-m15, tüm subcommand'lar, referans-prep
  komutları, kalite kapıları, envler. (Eski "sadece m01" metni kaldırıldı.)
- **★ GÖRSEL CİLA TAMAM ve `main`'de (2026-08-05, merge `604bf7e`, push).** 4 figür işi:
  (1) **m07 PCA** etiketleri ggrepel + eksen expand (kırpılma giderildi). (2) **m12 REVIGO** semantik-uzay
  **MDS scatter** (`semantic.R`, base R `cmdscale` + Lin uzaklık) kaynak başına. (3) **m14 operon** koordineli
  operon **bar figürü** (`operon.R`). (4) **m15 PPI** en büyük modüller **ağ figürü** (`ppi.R`, networkx spring
  layout — **numpy+scipy rnaforge-core'a pip** eklendi; beyaz zemin). Tüm yeni figürler **best-effort** (çekirdek
  tabloları bozmaz, hata loglanır) + rapora gömüldü. **366 test.** Rapor **6.2 MB, 19 gömülü figür.** Canlı
  doğrulandı (PPI modülleri: kapsül wcaJ/ugd · hücre-duvarı fts · ribozom rpl · glikoliz pgi/pykF; operon bar
  wca/wz UP; MDS semantik harita). Verdict SUSPECT değişmedi.
- **★ m15 PPI + COMMUNITY DETECTION TAMAM ve `main`'de (2026-08-05, merge `2163e86`, push). Dalga 2 #C (son).
  DOWNSTREAM DALGA 2 TAMAMEN BİTTİ (AMR/Operon/PPI).** STRING PPI alt-ağı + Louvain (networkx, rnaforge-core'a
  pip). Yeni `rnaforge ppi`; gate YOK. **366 test.**
  - **Kod:** `ppi.py` (STRING parser + sembol-join + `build_deg_network` + `louvain_communities` seed=42);
    `modules/m15_ppi.py`. `config.ppi` (PPIConfig: taxid/string_dir/min_score=700/min_community_size=3).
    Referans gitignore'lu `references/string/<taxid>/` (STRING v12, prep'te indirilir). Rapor "Protein
    Etkileşim Modülleri" bölümü + kanıt-skoru dürüst notu + STRING2023/Louvain2008 kaynak.
  - **GSE300731 canlı (kapstone):** 1226/1634 DEG ağda, 4072 kenar, 37 modül. Modüller tüm hikâyeyi
    protein-kompleks düzeyinde topladı: zarf-stres/kapsül (arn/cps/cpx/basR), efflux (acrABDE/emrABK),
    hücre-duvarı (fts/dac/lpoB), ABC transporter UP; karbon/yağ asidi/respirasyon DOWN. GO/KEGG/GSEA/AMR/
    operon ile birebir. Verdict SUSPECT değişmedi. Rapor 4.4 MB, **14 bölüm**.
  - Spec/Plan: `docs/superpowers/{specs,plans}/2026-08-05-m15-ppi*`.
- **★ m14 OPERON ANALİZİ TAMAM ve `main`'de (2026-08-05, merge `90112ef`, push). Dalga 2 #B.**
  İntergenik-mesafe sezgiseli (saf Python, GFF'ten; Ali harici araçtan yalın yerele döndü — Operon-mapper
  web/veri-dışarı, Rockhopper ağır). Yeni `rnaforge operon`; gate YOK; operonlar TAHMİN (dürüst not). **347 test.**
  - **Kod:** `operon.py` (`predict_operons` aynı-strand+gap≤max_gap; `aggregate_operon_de` koordineli);
    `modules/m14_operon.py`. `config.operon.max_gap`(50). Çıktı `operon/operons.tsv`. Figür yok (tablo-öncelikli).
    Rapor "Operon Analizi" bölümü + Moreno-Hagelsieb 2002 kaynak.
  - **GSE300731 canlı (yayın kalitesi):** 2781 operon (804 çok-genli, 217 koordineli DE). thrABC doğru.
    Koordineli operonlar tüm hikâyeyi topladı: **wca kolanik asit/kapsül (11 gen UP, log2FC +7.62)**, dcw
    peptidoglikan/hücre-duvarı (mraZ/fts/mur), ribozom UP, wec ECA; **his/arn/nrf DOWN**. GO/KEGG/GSEA/AMR
    ile birebir. Verdict SUSPECT değişmedi.
  - Spec/Plan: `docs/superpowers/{specs,plans}/2026-08-05-m14-operon*`.
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
- **Dalga 2 (bakteri overlay) ✅ TAMAM:** ~~#A AMR~~ (m13) · ~~#B Operon~~ (m14) · ~~#C PPI+community~~ (m15).
  ~~Cila (PPI/operon figür, REVIGO MDS, PCA etiket)~~ ✅ BİTTİ (604bf7e).
  **★ SIRADAKİ SEÇENEKLER (Ali seçecek):** (a) **Ökaryot yolu** (m04-euk salmon 2.3.4 + m05-euk tximport) —
  MVP'nin ikinci kolu, `reminder_rnaforge_eukaryote`. (b) GO agnostik (eggNOG). (c) Regulon/Sigma (E.coli-özel,
  agnostik-kırar — en sona). (d) README/PLAN v1.4 dokümantasyon güncellemesi (uzun süredir ertelendi).
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
