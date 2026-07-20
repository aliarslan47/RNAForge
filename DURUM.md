# DURUM — RNAForge

> Bu dosya "nerede kaldık" anlık görüntüsüdür. Tüm karar detayı Claude belleğindedir
> (`rnaforge-project` memory). Claude bunu anlamlı her durakta ve `/clear` öncesi günceller.

**Konum:** `/home/ali/rnaforge-pipeline/` (git deposu)
**GitHub:** `github.com/aliarslan47/RNAForge` — **PRIVATE**, remote `origin` (SSH), yalnız `main`
**Referans doküman:** `PLAN.md` **v1.2** (tek referans — Kural 1)
**Son güncelleme:** 2026-07-17

## Şu an nerede kaldık
- **Plan A 7/7 görev YAZILDI ve her görev review'dan geçti (2026-07-17). 43/43 test geçiyor.**
  Branch `feat/plan-a-temel-m01` (GitHub'a push edildi). Son commit `b1d03db`.
- **ANCAK final whole-branch review (Opus) "merge ÖNCESİ düzeltme gerekli" dedi.**
  Critical yok, **6 Important** var. Merge EDİLMEDİ. Düzeltmeler YAPILMADI — burada durduk.
- `rnaforge validate --config X --metadata Y` çalışıyor: config + metadata + design formülü
  doğrulama, FASTQ'dan platform tespiti, ONT/PacBio reddi. Saf Python (harici araç yok).
- **DİKKAT — kanıtlanmış yanlış iddia:** "resume/heartbeat altyapısı çalışıyor" DEĞİL.
  `new_run_dir` her çağrıda yeni zaman damgalı dizin açıyor → `state.json` hiç geri okunmuyor;
  `is_done`/`completed_modules` üretimde çağrılmıyor. Yani PLAN §15'in arkasında **kayıt var,
  dayanıklılık yok**. (Reviewer aynı `--run-id` ile iki kez koşturup kanıtladı.)
- Plan: `docs/superpowers/plans/2026-07-16-rnaforge-plan-a-temel-m01.md`
- **İlerleme ledger'ı: `.superpowers/sdd/progress.md`** (git-ignored; commit SHA'ları,
  Minor bulgular ve final review triyajı orada — devam ederken ÖNCE onu oku).

## Devam edince İLK İŞ: final review'un 6 Important'ı
1. **Resume unwired** (`state.py:18-22`, `cli.py`): `--run-id` verilse bile her koşu yeni dizin.
   Fix: mevcut `<stamp>_<run_id>` varsa onu kullan + `run_validation` `is_done` ise erken dön.
   Ayrıca `HEARTBEAT_INTERVAL_SECONDS=10` hiçbir yerde kullanılmıyor → ya kullan ya sil+PLAN düzelt.
   **m02 gelmeden ÖNCE çözülmeli** — 8 modül yanlış varsayıma yaslanmasın.
2. **`validation.log` sadece başarıda yazılıyor** (`m01_validate.py:54-107`): hata anında bellekte
   kalıp kayboluyor, `logs/` boş çıkıyor. Fix: satır satır aç+flush; `mark_done`'a log'u da ekle.
3. **Bozuk config bölümü ham traceback veriyor** (`config.py:123-127,140`): `library: "foo"` →
   `AttributeError`, `ConfigError` değil. Ticari üründe kabul edilemez. Fix: `_section()` +
   `_as_int(value, field)` helper'ları.
4. **`test_trimming_defaults_are_gentle` BOŞ test** (`tests/test_config.py:85`) — **KULLANICI KARARI
   BEKLİYOR** (aşağıda).
5. **`REQUIRED_REFERENCE` iki yerde** (`config.py:16-19` ve `m01_validate.py:21-24`): yönlendirme
   sözleşmesi çift kaynaklı → sürüklenme riski. Fix: m01 config'ten import etsin.
6. **Rank-deficient design'lar geçiyor** (`metadata.py:102-106`): batch condition'la tam confounded
   ya da tek seviyeli olsa bile kabul ediliyor → DESeq2'de kriptik hata. m06 gelmeden ucuz.

Minor'lar (fsync eksik, exit-2 stdout'a gidiyor, `modules/__init__.py`'de `from __future__` yok,
duplike FASTQ yolu kabul ediliyor) ledger'da.

## KULLANICI KARARI BEKLEYEN (plan-dayatmalı çelişki)
`tests/test_config.py:85` planın dayattığı şekilde `assert cfg.trimming.min_length >= 1` diyor.
Reviewer: bu **vacuous** — `config.py:133-134` zaten `min_length < 1`'de `ConfigError` atıyor,
yani assertion'ın başarısız olması İMKÂNSIZ. Williams 2016 gerekçesini koruduğu sanılan test
hiçbir şey korumuyor (tek işe yarar satır `aggressive_quality is False`).
**Soru: test `== 36` yapılsın mı, yoksa plan metni mi geçerli?** Ali cevaplamadı.

## Yöntem (devam ederken)
Subagent-driven development: her görev için `scripts/task-brief` ile brief çıkar →
implementer subagent (haiku) → `scripts/review-package` → task reviewer (sonnet) →
ledger'a kaydet. Script'ler:
`~/.claude/plugins/cache/claude-plugins-official/superpowers/6.1.1/skills/subagent-driven-development/scripts/`

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
