# m03 — Read Trimming (fastp, nazik) — Tasarım

**Tarih:** 2026-07-30
**Durum:** Onaylandı (Ali — "ne gerekiyorsa yap / devam edelim", otonom oturum)
**İlgili:** PLAN.md v1.3 §4.2, §5, §14; `docs/superpowers/specs/2026-07-30-m02-fastqc-design.md`

## 1. Problem

m02 ham okuma kalitesini raporladı ama düzeltmedi. m03 fastp ile okumaları hazırlar:
gerçek adapter kontaminasyonunu temizler ve çok kısa okumaları eler — ama **nazikçe**.

**Neden nazik (kritik, sezgiye aykırı):** Williams et al. 2016 (BMC Bioinformatics,
doi:10.1186/s12859-016-0956-2) agresif kalite trimming'in genlerin **>%10'unun**
ekspresyon tahminini bozduğunu gösterir (kısalan okumaların yanlış hizalanması). Sapmanın
çoğu **minimum uzunluk filtresiyle** kalkar. Agresif trimming ile koşulsaydı sonuçlar
referans yayından sapar ve bu sapma pipeline hatası sanılırdı (Katman B doğrulamasını bozar).

## 2. Onaylanan kararlar

1. **fastp nazik varsayılanlarla (PLAN §4.2):** adapter tespiti/temizliği AÇIK (fastp otomatik,
   `-A` verilmez); agresif kalite trimming (sliding-window `-r/--cut_right`) KAPALI; minimum
   uzunluk filtresi AÇIK ve zorunlu (`-l <min_length>`); kalite filtresi fastp varsayılanı
   (nazik, `-Q` verilmez).
2. **Veri kapısı: `survival_rate` — FAIL kapısı.** `after_filtering.total_reads /
   before_filtering.total_reads` profil eşiğinin (`survival_rate`, prok/euk = 0.50) altındaysa
   **FAIL → koşu durur** (sonuç GEÇERSİZ: kütüphanenin yarısından çoğu gittiyse kantifikasyon
   çöp üretir). Bu, m03'ü m02'den ayırır: **m03 durabilir** (ilk gerçek veri-FAIL kapısı),
   m02 (ham QC diagnostik) asla durmaz. Eşik zaten `profiles/{prokaryote,eukaryote}.yml`'de var.
3. **Paired-end destekli.** metadata `fastq_2` taşır, m01 tespit eder; fastp native paired
   (`-I/-O`). Tek-uçlu ve çift-uçlu ikisi de MVP'de.
4. **Config bağlama:** `min_length ← config.trimming.min_length` (varsayılan 36).
   `config.trimming.aggressive_quality: true` → fastp'ye `-r` eklenir (BİLİNÇLİ agresif; nazik
   varsayılandan sapma). Nazik varsayılan daima güvenli taraftadır.
5. **Ön koşul: m01 tamamlanmış olmalı** (doğrulanmış girdi; m02 ön koşul DEĞİL — QC diagnostik,
   trimming için gerekli değil). Aynı `run_dir`/`--run-id` deseni (m02 ile aynı).

### Reddedilen seçenekler
- **survival_rate = WARN:** Yarıdan fazla okuma kaybı sonucu GEÇERSİZ kılar (kalan %40'ı
  saymak makul görünen sahte sonuç üretir). Dual politika: geçersiz = FAIL.
- **Agresif kalite varsayılanı:** Williams 2016 — doğruluğu bozar; reddedildi (PLAN §4.2).
- **m02 ön koşulu zorlamak:** QC trimming için gerekli değil; gereksiz bağ.

## 3. Arayüz (public sözleşme)

```python
def run_trim(config: Config, metadata_path: Path, run_dir: Path,
             force: bool = False) -> dict: ...
```

- **Ön koşul:** m01 bu `run_dir`'de done değilse `ValueError` (m02 deseni).
- **Dönüş:** özet dict: `n_samples`, örnek başına `survival_rate`/`reads_before`/`reads_after`,
  `gate_counts`, `resumed?`.

## 4. Bileşenler

- **`rnaforge/fastp.py`** — saf parser + runner ayrımı (m02 `fastqc.py` deseni):
  - `run_fastp(fastq_1, out_dir, min_length, fastq_2=None, aggressive_quality=False,
    env="rnaforge-qc") -> FastpResult` — fastp'yi shell'ler; trimlenmiş FASTQ + `fastp.json`
    + `fastp.html` üretir; parse edilmiş `FastpResult` döner.
  - `parse_fastp_json(json_text: str) -> FastpResult` — *saf* string parser (I/O yok, hızlı test).
  - `FastpResult` frozen dataclass: `reads_before: int`, `reads_after: int`,
    `survival_rate: float` (`reads_after/reads_before`, before=0 ise 0.0), `out1: Path|None`,
    `out2: Path|None`.
  - `FastpRunError(RuntimeError)`.
- **`rnaforge/modules/m03_trim.py`** — `run_trim` orkestrasyonu (m02 `run_qc` deseni:
  `RunState` resume+heartbeat, `trim.log`, çıktı yazımı, kapı üretimi) + `build_trim_gates(
  results: dict[str, FastpResult], profile: Profile) -> list[GateResult]`.
- **`rnaforge/cli.py`** — `trim` subcommand + `_cmd_trim`.

## 5. Veri akışı ve çıktı yapısı

```
runs/<ts>_<run_id>/
├── trimmed/<sample_id>/         # trimlenmiş FASTQ (out1[, out2]) + fastp.json + fastp.html
├── statistics/trimming_statistics.json  # örnek başına survival/before/after
├── quality/gates.json           # m03 kapısı EKLENİR (m01/m02'ye dokunma)
└── logs/trim.log
```

fastp komutu (nazik): `fastp -i <in1> -o <out1> [-I <in2> -O <out2>] -l <min_length>
-j fastp.json -h fastp.html [-r]`. `-r` yalnız `aggressive_quality: true` ise.

## 6. Kapı eşlemesi (`build_trim_gates`)

Örnek başına tek `survival` kapısı DEĞİL — örnekler boyunca tek toplu `survival_rate` kapısı:
- Herhangi bir örnek `survival_rate < profile.threshold("survival_rate")` ise **FAIL**;
  `samples` = eşiğin altındaki örnekler; `measured` = en düşük survival; `threshold` = profil eşiği.
- Hepsi eşikte/üstündeyse **PASS**.
- `remedy`: "Bu örneklerde okumaların çoğu min_length filtresini geçemedi — okuma uzunluğunu/
  platformu doğrulayın; yanlış veri veya çok kısa okumalar olabilir."
- Profil eşiği config'ten ezildiyse `overridden=True` (mevcut `Profile.overrides()` üzerinden).

## 7. Hata yönetimi (Kural 7)
- **m01 done değil:** `ValueError` ("run m01 (validate) first").
- **fastp binary/env yok, sıfırdan farklı exit, çıktı yok:** `FastpRunError` (env, cmd, exit, stderr).
- **fastp.json bozuk/eksik alan:** `parse_fastp_json` net hata.
- **FAIL kapısı — yazma sırası (m01 deseni):** tüm örnekler için fastp koşulur ve
  `trimming_statistics.json` yazılır → kapılar `write_gate_results` ile gates.json'a yazılır →
  EN SON `raise_if_failed`. Böylece FAIL'de bile hem trimlenmiş FASTQ hem istatistik hem
  gates.json diskte kalır; teşhis raporu bunları okur. `state.mark_done` yalnız FAIL yoksa
  çağrılır (FAIL'li koşu "tamamlandı" işaretlenmez).
- **Resume:** m03 done + `trimming_statistics.json` varsa tekrar koşmaz (`resumed: True`); `--force`.

## 8. Test stratejisi (TDD)
**Birim (fixture'lara karşı, env'siz):**
- `parse_fastp_json`: fixture json → doğru before/after/survival; before=0 → survival 0.0.
- `build_trim_gates`: survival < eşik → FAIL + doğru `samples`/`measured`; hepsi ≥ eşik → PASS;
  eşik config'ten ezilmiş → `overridden=True`.
- `run_trim` (fastp monkeypatch): trimmed/<sample>/ + trimming_statistics.json + gates yazılır;
  m01 kapıları korunur; resume FastQC gibi fastp'yi tekrar çağırmaz; m01 done değilse ValueError;
  düşük survival → GateFailure + gates.json'da FAIL kaydı.
- CLI: `trim` subcommand exit 0 (validate→trim), verdict basılır; düşük survival'da exit 1.

**Entegrasyon (gerçek fastp, `rnaforge-qc` yoksa skip):**
- Sentetik FASTQ → gerçek `run_fastp` → trimlenmiş FASTQ + parse edilebilir fastp.json;
  survival ~1.0 (nazik trimming kısa-olmayan okumaları elemez). Parser'ın fastp 1.3.6 uyumu.

## 9. Kapsam dışı (Faz 2+ / sonraki modüller)
- MultiQC (fastp+FastQC+Salmon toplayıcı, quant sonrası, Faz 2+).
- Trim-sonrası FastQC (MVP'de yalnız ham QC; m02 tek QC adımı).
- UMI/dedup, poly-G/poly-X aşırı işleme (nazik politika dışı).
- Kantifikasyon (m04) — m03 yalnız trimler; hizalama/sayım sonraki modül.
