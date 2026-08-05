# m16 — Sequencing QC (rRNA% + strandedness) · Tasarım Spec'i

**Tarih:** 2026-08-05 · **Dal:** `feat/m16-seqqc` · **Referans:** kalite kapıları çerçevesi (m01)
**Amaç:** İki bakteri-kritik QC metriği + **iki WARN kapısı**: rRNA depletion verimi (SortMeRNA) ve
strandedness doğrulama (RSeQC). Kötü girdiyi yakalar → "makul görünen sahte sonuç" korumasını güçlendirir
([[feedback_dogruluk_kontrol]]). Ali kararı: ayrı araçlar (SortMeRNA+RSeQC), uyuşmazlık → WARN (FAIL değil).

## 1. rRNA% (SortMeRNA)

- Trimlenmiş okumalar (m03) rRNA referans DB'sine hizalanır → **% rRNA** = rRNA / toplam. Depletion verimi.
- **Subsample**: örnek başına ilk ~200k okuma (hız; QC tahmini yeterli). `qc.rrna_subsample`(200000).
- `run_sortmerna(r1, r2, ref_fasta, workdir, n_subsample, env) -> str` + `parse_sortmerna_log(log) -> float`
  (aligned yüzdesi). Ortalama örnek rRNA fraksiyonu koşu metriği.
- **Referans (indirme YOK, agnostik):** `rrna_fasta_from_reference(genome_fa, gff, out)` — GFF `rRNA`
  feature koordinatlarından genom rRNA dizilerini çıkarır (− strand ters-tümler). SortMeRNA `--ref` bu.
  Referans organizmanın rRNA'sına hizalama = depletion için isabetli; büyük SILVA indirmesi gerekmez.

## 2. Strandedness (RSeQC)

- `gff_to_bed(gff, out_bed)`: GFF `gene` feature'ları → BED12 (prokaryotta gen = tek blok). Saf Python.
- `run_infer_experiment(bam, bed, env) -> str` + `parse_infer_experiment(out) -> (strand, fwd, rev)`:
  sense/antisense oranları → **çıkarılan strand** (`unstranded`/`stranded`/`reverse`; eşik ~0.8 baskınlık).
- Çıkarılan strand `config.library.strandedness` ile karşılaştırılır.

## 3. Kapılar (`gates.py` sözleşmesi)

- `rrna_fraction`: WARN if ortalama rRNA% > profil eşiği (`max_rrna_fraction`, varsayılan 0.20).
- `strandedness_match`: WARN if çıkarılan ≠ beyan edilen (mesaj: çıkarım vs config). Eşleşiyorsa PASS.
- Her ikisi de m16 tarafından gates.json'a yazılır → **güvence kartına akar** (WARN → verdict SUSPECT).
  FAIL üretmez.

## 4. Modül (`rnaforge/modules/m16_seqqc.py`) + CLI

- Ön koşul **m04** (BAM `quantification/<sample>/aligned.sorted.bam`) + m03 (trimlenmiş okuma). Resume/heartbeat.
- Çıktı `seqqc/{rrna.json, strand.json, genes.bed}` + gates + `statistics/seqqc_statistics.json`
  (mean_rrna_fraction, per-sample rRNA%, inferred_strandedness, declared_strandedness, match).
- Yeni `rnaforge seqqc` subcommand.

## 5. Rapor (`report_html.py`)

"Kalite ve İşleme" bölümüne satırlar: **rRNA %** (ortalama + örnek-başı), **çıkarılan strandedness**,
**beyan ile uyum** (✓/uyuşmuyor). WARN durumları güvence kartında zaten görünür.

## 6. Config & env

Yeni **rnaforge-seqqc** env (sortmerna 7.0 + rseqc 5.0.5). Config: `qc.rrna_subsample`(200000).
Profil eşiği `max_rrna_fraction`(0.20) → `profiles/{prokaryote,eukaryote}.yml`. `qc` KNOWN_TOP_LEVEL_KEYS'e
(veya mevcut `quality`/`library`'e). SortMeRNA DB `references/sortmerna/` (prep).

## 7. Doğrulama

- **Birim** (`rnaforge-core`): `gff_to_bed` (BED12, tek blok, 0-tabanlı), `parse_sortmerna_log` (%),
  `parse_infer_experiment` (unstranded/forward/reverse eşiği), kapı mantığı (WARN eşikleri), orkestrasyon
  (m04 ön koşul, resume). Araç çağrıları monkeypatch.
- **GSE300731 canlı smoke:** rRNA% (rRNA-depletion kütüphanesi → düşük beklenir, WARN yok); strandedness
  çıkarımı (config "unstranded" beyanıyla karşılaştır). Rapor satırları; güvence kartı.

## 8. İş akışı

spec → writing-plans → rnaforge-seqqc env + SortMeRNA DB → TDD (~8 task) → GSE300731 smoke → merge → DURUM/bellek.
