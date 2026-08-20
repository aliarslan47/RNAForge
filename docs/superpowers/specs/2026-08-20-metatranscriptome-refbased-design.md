# Metatranskriptom Kolu — Referans-Tabanlı, Kısa-Okuma (M1) Tasarımı

**Tarih:** 2026-08-20
**Durum:** Onaylandı (Ali, brainstorm) — implementasyon planı yazılacak
**Kapsam:** Milestone 1 — referans-tabanlı, kısa-okuma (Illumina) metatranskriptom.
De novo assembly (M2), uzun-okuma (ONT), taxon-stratified fonksiyon ve eşli-metagenom
(DNA) normalizasyonu bu spec'in DIŞINDA (sonraki milestone'lar).

## 1. Amaç ve gerekçe

RNAForge'a bir **topluluk (community) RNA-seq** kolu ekliyoruz: tek organizma yerine
mikrobiyal bir topluluğun metatranskriptomu. Alan iki paradigmaya ayrılır — referans-tabanlı
(assembly-free) ve de novo assembly. **M1 = referans-tabanlı, kısa-okuma**, çünkü:

- RNAForge'un mevcut altyapısını en çok yeniden kullanır (SortMeRNA, fastp, DESeq2, KEGG/GO,
  rapor, kapılar), en düşük risk.
- 2025 benchmark'ı (MT-Enviro, ISME Communications) iyi-karakterize topluluklarda
  referans-tabanlı + TPM'i önerir.
- Hazır bir nf-core ürünü M1'i karşılamıyor: `nf-core/metatdenovo` yalnız de-novo (Nextflow,
  bizim M2 şablonumuz); `nf-core/taxprofiler` yalnız metagenom-taksonomi (fonksiyon/DE yok);
  HUMAnN nf-core değil ve monolitik. Native modüler inşa, aile standardıyla (RNA/Virus/Bac
  Forge tek rapor, güven kartı, DAG, TR/EN) entegre tek çözüm.

## 2. Temel mimari karar: üçüncü `organism_type`

RNAForge zaten `organism_type` (`prokaryote`/`eukaryote`) ve `read_type` (`short`/`long`)
üzerinden dallanır (`m04_quant.py`, `m05_counts.py` router deseni). Metatranskriptomu
**üçüncü bir `organism_type = "metatranscriptome"`** olarak ekliyoruz.

- `ORGANISM_TYPES` → `("prokaryote", "eukaryote", "metatranscriptome")` (`config.py:9`).
- `REQUIRED_REFERENCE["metatranscriptome"]` = `("gene_catalog_fasta", "catalog_annotation")`
  (aşağıda §6). Kraken2 DB ve rRNA DB config'te ayrı bölümlerde (zorunlu-referans değil,
  taksonomi/temizlik yardımcıları).
- `m04`/`m05` router'larına, mevcut `eukaryote` dalının yanına `metatranscriptome` dalı eklenir.
- Downstream (`m06`+) **kod değişmeden** ortak sayım matrisinde buluşur — ökaryot/uzun-okuma
  kollarındaki aynı sözleşme.

## 3. Uçtan uca akış

```
m01 validate         → organism_type=metatranscriptome tanınır (PAYLAŞILAN, ufak ekleme)
m02 FastQC           → PAYLAŞILAN (değişmez)
m03 fastp            → PAYLAŞILAN (değişmez)
m_rrna rRNA DEPLETION → YENİ: SortMeRNA ile rRNA'yı ÇIKAR (sadece ölçme değil) → rRNA'sız FASTQ
m_tax  TAKSONOMİ     → YENİ: Kraken2 + Bracken (kim var + göreli bolluk); diagnostik, kapı yok
m04 KATALOG-KUANT    → YENİ DAL: Bowtie2/BBMap ile gen kataloğuna hizala; mapping DİAGNOSTİK
m05 KATALOG-SAYIM    → YENİ DAL: featureCounts → KO/gen sayım matrisi (ortak counts.tsv sözleşmesi)
m06 DESeq2           → PAYLAŞILAN: topluluk-düzeyi diferansiyel ekspresyon (özellik=gen/KO)
m09-m12 KEGG/GO/GSEA → PAYLAŞILAN: topluluk fonksiyonları (KO→KEGG pathway)
m08 rapor            → PAYLAŞILAN + taksonomi-kompozisyon bölümü + metatranskriptom yöntem/atıf
```

**Paylaşılan (değişmez):** m02, m03, m06, m07, m09–m12.
**Paylaşılan + ufak ekleme:** m01 (organism_type), m08 (rapor bölümü).
**Yeni:** rRNA-depletion işlem adımı, taksonomi modülü, m04/m05 metatranscriptome dalı, profil,
referans hazırlığı, config bölümleri.

## 4. Bileşenler (izole birimler)

### 4.1 rRNA depletion — İŞLEM adımı (kritik fark)
Mevcut `seqqc.py`'deki `run_sortmerna` yalnız **ölçer** (`--num_alignments 1`, rRNA% döndürür).
Metatranskriptomda rRNA'yı **çıkarıp** rRNA'sız FASTQ üretmemiz gerekir.

- **Yeni fonksiyon** `rnaforge/rrna_deplete.py::run_sortmerna_deplete(reads, rrna_db, workdir, ...)`
  → SortMeRNA `--fastx --aligned <rrna> --other <non_rRNA>` ile paired-aware çalışır; **`--other`
  (rRNA'sız) çıktıyı** downstream sözleşme yoluna yazar; depletion fraksiyonunu döndürür.
- Mevcut `run_sortmerna` (ölçüm) DEĞİŞMEZ — yeni fonksiyon ayrı, tek-sorumluluk.
- **rRNA referansı:** community'de tek genom yok → GFF'ten çıkarma çalışmaz. SortMeRNA'nın
  paketlenmiş rRNA DB'leri (smr_v4 default) veya SILVA/Rfam kullanılır; config'te `rrna.db_fasta`.
- **Modül** `modules/m_rrna_deplete.py`: per-sample döngü, atomic state, heartbeat, resume
  (aile deseni). Çıktı: `rrna_depleted/<sid>/*.fastq.gz` + `depletion_stats.json`.
- **Kapı:** depletion verimi WARN (profil `rrna_depletion_rate`), asla FAIL.

### 4.2 Taksonomik profilleme (yeni modül)
- **Yeni** `rnaforge/kraken2.py`: `run_kraken2(reads, db, ...)` + `run_bracken(...)` +
  parser'lar (`parse_kraken2_report`, `parse_bracken`) → tür/cins düzeyi göreli bolluk tablosu.
- **Modül** `modules/m_taxonomy.py`: rRNA'sız okumalar girdisi; **tamamen diagnostik, kapı yok**.
  Çıktı: `taxonomy/<sid>.bracken` + birleşik `taxonomy/abundance_matrix.tsv`.
- Env: yeni `rnaforge-meta` (kraken2, bracken).

### 4.3 Katalog-tabanlı kuantifikasyon (m04/m05 yeni dal)
- **m04** `_quant_meta` (yeni): rRNA'sız okumalar → Bowtie2 (mevcut `bowtie2.py` yeniden
  kullanılır) veya BBMap ile **gen kataloğuna** hizala. `alignment_rate` **DİAGNOSTİK — FAIL
  kapısı YOK** (kataloglar eksiktir, düşük oran normaldir; uzun-okuma kolundaki aynı felsefe).
  Çıktı BAM `quantification/<sid>/aligned.sorted.bam` (m05 sözleşmesi).
- **m05** `_counts_meta` (yeni): `featureCounts` katalog anotasyonu (GFF/SAF) ile → gen/KO sayımı.
  Katalog anotasyonu KO/fonksiyon taşır → sayımlar KO'ya toplanabilir. `assignment_rate`
  DİAGNOSTİK. Ortak `_write_count_outputs` (counts/tpm/fpkm) yeniden kullanılır.
- Router: `m04_quant.py`/`m05_counts.py` içine `if config.organism_type == "metatranscriptome":`
  dalı (mevcut `eukaryote` dalının yanına), read_type=short.

### 4.4 Downstream (paylaşılan, değişmez)
- **m06 DESeq2:** özellik = gen/KO; topluluk-düzeyi DE. Değişmez.
- **m09–m12 KEGG/GO/GSEA/REVIGO:** KO→KEGG pathway zenginleştirme. Katalog KO anotasyonu
  `enrichment` altyapısına köprülenir (ökaryot transkriptom-sembol köprüsü deseni). Değişmez kod.
- **m08 rapor:** yeni **taksonomi-kompozisyon bölümü** (Bracken top-N tür + figür) +
  `organism_type=metatranscriptome` rozeti + metatranskriptom yöntem metni (TR+EN) + atıflar
  (SortMeRNA, Kraken2/Bracken, ilgili benchmark). Kullanılmayan aracı atıflamaz (dürüstlük).

## 5. Kalite profili
Yeni `rnaforge/profiles/metatranscriptome.yml` (`permissive: true`, DAMGALI —
`prokaryote_long.yml` deseni):

- `alignment_rate`: yok / çok düşük eşik → **FAIL kapısı YOK** (katalog eksikliği doğal). Yalnız
  katastrofik (~0) durum için opsiyonel çok-düşük eşik, WARN.
- `rrna_depletion_rate`: WARN (verim düşükse şüpheli, geçersiz kılmaz).
- `assignment_rate`: WARN.
- `replicate_correlation`: 0.75 (korunur — replika tutarlılığı topluluğda da anlamlı).
- Gerekçe (belleğe): *uydurma eşik kapı sistemini itibarsızlaştırır* → temsili bir metatranskriptom
  seti gelene dek permissive + damgalı (ökaryot/uzun-okuma profili felsefesi).

`quality.profile_name_for` metatranscriptome + short → `metatranscriptome` döndürür (mevcut
mantık zaten short'ta organism_type'ı döndürür; ek kod gerekmez).

## 6. Referans hazırlığı ve config

### 6.1 config bölümleri (yeni)
- `REQUIRED_REFERENCE["metatranscriptome"] = ("gene_catalog_fasta", "catalog_annotation")`
- `Reference`'a alanlar: `gene_catalog_fasta`, `catalog_annotation` (GFF/SAF; KO/fonksiyon taşır).
- Yeni bölüm `taxonomy`: `kraken2_db` (dizin), `bracken_read_len`, `bracken_level` (S), `env`.
- Yeni bölüm `rrna`: `db_fasta` (SortMeRNA rRNA DB), `env`.
- `KNOWN_TOP_LEVEL_KEYS`'e `taxonomy`, `rrna` eklenir.

### 6.2 prepare_references.sh genişletme
- Kraken2 DB indir/derle (standart veya kullanıcı-DB; parametreli `--kraken2-db`).
- Bracken kmer dağıtımı üret.
- Gen kataloğu (ör. insan bağırsak için UHGG; config ile değiştirilebilir) + KO anotasyonu.
- rRNA DB (SortMeRNA smr default veya SILVA). Hepsi gitignore + `.sha256`.

## 7. CLI ve orkestrasyon
- Yeni subcommand'lar: `rnaforge rrna-deplete`, `rnaforge taxonomy` (aile deseni; her modül tek
  komut). Mevcut `quant`/`counts` metatranscriptome dalını otomatik seçer (organism_type'tan).
- `rnaforge run` orkestratörü metatranscriptome için sıralamaya rRNA-depletion + taxonomy'yi
  organism_type'a göre ekler.

## 8. Test stratejisi (TDD)
Her birim ayrı test edilir (aile disiplini):
- `rrna_deplete`: SortMeRNA `--other` çıktısını üretir; depletion fraksiyonu doğru parse (gerçek
  küçük SortMeRNA entegrasyon testi, conda yoksa skip).
- `kraken2`/`bracken`: rapor/bracken parser birim testleri; küçük gerçek Kraken2 çalıştırma (skip'li).
- `m04/m05 meta dal`: küçük katalog + sentetik okumalarla Bowtie2→featureCounts, KO'ya toplama doğru.
- Router regresyonu: prokaryot/ökaryot yolları DEĞİŞMEDİ.
- Profil: metatranscriptome yüklenir, permissive damgası, kapı davranışı (align FAIL yok, WARN'lar).
- E2E smoke: sentetik küçük topluluk → `rnaforge run` (rrna→tax→quant→counts→de→report), conda yoksa skip.

## 9. Biyolojik doğrulama (aile disiplini — kod sonrası)
Typhi/airway disipliniyle: 2 paralel ajan (ENA + literatür) yayınlanmış, bilinen-DE'li bir
kısa-okuma metatranskriptom seti seçer (aday: iHMP/HMP2 IBD bağırsak metatranskriptomu; ya da
tanımlı mock-topluluk). Uçtan uca koşulur, konkordans kontrol edilir. Opsiyonel: aynı sette
nf-core/metatdenovo benchmark olarak koşulur. Veri seçimi kör indirme değil — o aşamada birlikte.

## 10. Kapsam dışı (YAGNI — sonraki milestone'lar)
- **M2:** de novo assembly kolu (rnaSPAdes/MEGAHIT → Prodigal/TransDecoder → eggNOG/KofamScan →
  featureCounts → DE); nf-core/metatdenovo araç seti şablon.
- Uzun-okuma (ONT) metatranskriptom.
- Taxon-stratified fonksiyon (HUMAnN-tarzı "hangi tür bu değişimi sürüyor").
- Eşli metagenom (DNA) ile RNA/DNA normalizasyonu.

## 11. Riskler
- **DB boyutu:** Kraken2 standart DB (~50-100GB) ve gen kataloğu büyük — prepare-references
  parametreli, küçük DB (ör. MiniKraken/16GB) opsiyonu belgelenir.
- **Katalog eksikliği → düşük mapping:** tasarımla diagnostik (FAIL yok); rapor bunu dürüstçe
  damgalar. Taksonomi (Kraken2) tam-eşleşme gerektirmez, tamamlayıcı sinyal verir.
- **KO köprüsü:** katalog anotasyon formatı değişken; parser esnek + yüksek-sesle hata (sessiz yok).
