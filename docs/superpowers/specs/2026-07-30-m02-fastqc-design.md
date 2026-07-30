# m02 — Quality Control (FastQC) — Tasarım

**Tarih:** 2026-07-30
**Durum:** Onaylandı (Ali, brainstorm oturumu)
**İlgili:** PLAN.md v1.3 §5, §14; `docs/superpowers/specs/2026-07-20-kalite-kapilari-design.md`

## 1. Problem

Pipeline'ın ilk gerçek biyoinformatik adımı. Her örneğin **ham** FASTQ'suna FastQC
koşup okuma kalitesini görünür kılar: müşteri ve sonraki modüller, verinin trimming
öncesi durumunu (kalite profili, adapter içeriği, aşırı-temsil edilen diziler) burada
görür. m01 girdinin *biçimini* doğruladı; m02 girdinin *kalitesini* raporlar.

**Kritik ayrım:** Ham okumada kötü kalite BEKLENEN bir durumdur ve m03 (fastp trimming)
tam da bunu düzeltmek içindir. Bu yüzden m02 **koşuyu asla durdurmaz** — ham QC
diagnostiktir (PLAN §4.3, "ham QC diagnostiktir"). FAIL-durdurma, sonucu gerçekten
geçersiz kılan downstream veri kapılarına (m04 hizalama oranı) aittir.

## 2. Onaylanan kararlar (brainstorm 2026-07-30)

1. **Kapı politikası: yalnız PASS/WARN.** FastQC modül bayrakları bizim kapı
   sistemimize eşlenir: `PASS→PASS`, `WARN→WARN`, **`FAIL→WARN`**. m02 **hiçbir zaman
   `FAIL` üretmez** → koşu m02'de durmaz. Kötü ham kalite → ŞÜPHELİ damga, GEÇERSİZ değil.
   İkili politika (kalite kapıları spec §2.2) korunur: durdurma yalnız sonuç geçersizse.
2. **Orkestrasyon: ayrı `qc` subcommand.** `rnaforge qc` komutu, m01'in `validate`
   komutunu birebir yansıtır; aynı `run_dir`'de resume-aware. Tam pipeline `run` zinciri
   sonraki modüller geldikçe eklenecek (şimdilik YAGNI).
3. **Test: parser birim + 1 entegrasyon.** Parser (`summary.txt`/`fastqc_data.txt` →
   dataclass) commit'li küçük fixture'lara karşı TDD ile test edilir (hızlı, deterministik,
   env'siz). Ayrıca gerçek FastQC'yi shell'leyen 1 entegrasyon testi — `rnaforge-qc` env
   / `fastqc` binary yoksa `skip` (yanlış alarm vermesin).

### Reddedilen seçenekler
- **m02'de FAIL-durdurma:** Ham okumada FAIL beklenen; trimming'e şans vermeden durdurmak
  diagnostik amaçla çelişir ve çok sayıda meşru koşuyu boşuna reddeder.
- **Kapı yok, sadece topla:** "Her modül kendi veri kapısıyla" ilkesinden sapar; ham QC
  durumunun güvence kartında görünmesi gerekir.
- **Tek `run` zinciri şimdi:** m03+ henüz yok; erken orkestrasyon.

## 3. Arayüz (public sözleşme)

```python
def run_qc(config: Config, metadata_path: Path, run_dir: Path,
           force: bool = False) -> dict: ...
```

- **Girdi:** m01'in doğruladığı config + metadata; her örneğin `fastq_1` (varsa `fastq_2`).
- **Ön koşul:** m02, m01'in bu `run_dir`'de tamamlanmış olmasını bekler. m01 yapılmamışsa
  net hata (`ValueError`) ver — sessizce ham FASTQ'ya koşma (girdi doğrulanmamış olabilir).
- **Çıktı (dönüş):** özet dict (örnek sayısı, örnek başına FastQC modül bayrakları özeti,
  toplu kapı durumu). m01'in summary sözleşmesiyle uyumlu; `resumed` bayrağı resume'de.

## 4. Bileşenler ve dosya düzeni

Tek sorumluluk sınırlarıyla iki yeni birim:

- **`rnaforge/fastqc.py`** — *saf* parser + runner ayrımı:
  - `run_fastqc(fastq: Path, out_dir: Path, env: str) -> Path`: FastQC'yi `rnaforge-qc`
    env'de shell'ler, üretilen `<name>_fastqc.zip` yolunu döner. Yan etki burada izole.
  - `parse_fastqc_zip(zip_path: Path) -> FastQCReport`: zip'ten `summary.txt` +
    `fastqc_data.txt` okur → `FastQCReport` (saf, I/O yalnız okuma; shell yok → hızlı test).
  - `FastQCReport` dataclass (frozen): `modules: dict[str,str]` (modül→PASS/WARN/FAIL),
    `basic_stats: dict` (total sequences, %GC, sequence length (FastQC'nin verdiği
    aralık string'i, ör. "150"), encoding). FastQC ortalama uzunluk VERMEZ — aralık verir;
    ortalama okuma uzunluğu m01'in işidir, burada tekrarlanmaz.
- **`rnaforge/modules/m02_qc.py`** — `run_qc` orkestrasyonu (m01 desenini izler:
  `RunState` resume+heartbeat, `qc.log` satır-satır flush, çıktı yazımı, kapı üretimi).

**Yeni birim yaratmayan değişiklikler:** `rnaforge/cli.py`'ye `qc` subcommand + `_cmd_qc`.

## 5. Veri akışı ve çıktı yapısı (PLAN §14)

```
runs/<ts>_<run_id>/
├── raw_qc/<sample_id>/          # her örneğin ham FastQC çıktısı (HTML + zip KORUNUR)
├── statistics/qc_statistics.json # parse edilmiş: örnek başına modül bayrakları + basic stats
├── quality/gates.json          # m02 kapıları EKLENİR (m01'inkine dokunulmaz)
└── logs/qc.log
```

Adım-bazlı bölme (PLAN §14): örnek kırılımı `raw_qc/` *içinde* korunur. Ham FastQC HTML'i
silinmez — müşteri/rapor onu doğrudan gösterebilir.

## 6. Kapı eşlemesi

FastQC'nin ilgilendiğimiz modülleri **isimli kapılara** eşlenir; her kapı örnekler
boyunca toplanır (`GateResult.samples` = bayrak veren örnekler), durum = örnekler
arasındaki en kötü bayrak (`FAIL→WARN` eşlemesiyle). MVP curated set:

| Kapı adı | FastQC modülü | remedy (WARN'da) |
|---|---|---|
| `per_base_quality` | Per base sequence quality | m03 trimming kalite profilini iyileştirir |
| `adapter_content` | Adapter Content | m03 fastp adapter'ı otomatik temizler |
| `overrepresented` | Overrepresented sequences | rRNA/adapter kontaminasyonu; kütüphane kimyasını doğrula |
| `gc_content` | Per sequence GC content | beklenmedik GC → kontaminasyon/tür karışımı olabilir |

- Tüm FastQC modüllerinin bayrakları `qc_statistics.json`'a yazılır (curated set kapıya
  dönüşür; gerisi bilgi olarak durur). Kapı `measured/threshold` alanları FastQC'de
  kategorik olduğu için boş; `message` hangi örneklerin ne bayrağı verdiğini söyler.
- Kapılar `write_gate_results` ile yazılır (atomik, m01'inkine dokunmaz). Güvence kartı
  bunları otomatik toplar (m02 WARN → kart `SUSPECT`).

## 7. Hata yönetimi (Kural 7: sessiz devam yok)

- **m01 yapılmamış:** `ValueError` — "run m01 (validate) first".
- **FastQC binary yok / env eksik:** `run_fastqc` net hata verir (hangi env, ne çalıştırıldı,
  exit kodu, stderr). Sessizce boş rapor üretme.
- **FastQC sıfırdan farklı exit / zip yok:** yüksek sesle hata; o örnek için sahte "PASS"
  uydurma.
- **Bozuk/eksik `fastqc_data.txt`:** `parse_fastqc_zip` hangi alanın eksik olduğunu söyleyen
  net hata verir.
- **Resume:** m02 bu `run_dir`'de `is_done` ve `qc_statistics.json` varsa işi tekrarlamaz,
  önceki özeti döner (`resumed: True`); `--force` ile ezilir (m01 deseni).

## 8. Test stratejisi (TDD)

**Birim (fixture'lara karşı, env'siz):**
- `parse_fastqc_zip`: commit'li küçük `summary.txt` + `fastqc_data.txt` fixture → doğru
  `modules` + `basic_stats`.
- Kapı eşlemesi: `FAIL→WARN`, en-kötü-bayrak toplama, `samples` doğru dolar.
- Curated set: yalnız 4 kapı üretilir; diğer modüller `qc_statistics.json`'da kalır.
- Hata yolları: eksik alan → net hata; m01 yapılmamış → `ValueError`.
- Resume: ikinci koşu FastQC'yi tekrar çağırmaz (`resumed: True`).
- CLI: `qc` subcommand exit 0; verdict satırı basılır.

**Entegrasyon (gerçek FastQC, `rnaforge-qc` yoksa skip):**
- Sentetik FASTQ → gerçek `run_fastqc` → `raw_qc/<sample>/` zip üretilir → parse edilir →
  beklenen modüller mevcut. Parser'ın gerçek FastQC 0.12.1 çıktısıyla uyumunu kanıtlar.

Fixture'lar sentetik (Kural 8: gerçek/müşteri verisi yok); FastQC çıktısı da sentetik
FASTQ'dan üretilmiş küçük dosyalar.

## 9. Kapsam dışı (Faz 2+ / sonraki modüller)
- MultiQC toplama (PLAN: quant sonrası, Faz 2+).
- Trimming (m03) — m02 yalnız raporlar, düzeltmez.
- Trim sonrası ("post-QC") FastQC — MVP'de yalnız ham QC.
- FastQC HTML'inin rapora gömülmesi (m08'in işi).
