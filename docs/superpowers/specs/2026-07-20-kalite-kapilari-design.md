# Kalite Kapıları ve Güvence Raporlaması — Tasarım

**Tarih:** 2026-07-20
**Durum:** Onaylandı (Ali, brainstorm oturumu)
**İlgili:** PLAN.md v1.3, `docs/superpowers/specs/2026-07-16-rnaforge-mvp-design.md`

## 1. Problem

Pipeline'ın kendisi doğru olsa bile, kötü girdi **makul görünen sahte sonuç** üretir.
En sinsi örnek hizalama oranıdır: %8 hizalama da bir count matrisi üretir, DESeq2 ona
bir p-değeri verir, rapor basılır. Hiçbir adımda hata yoktur — sadece sonuç yanlıştır.

Bu, PLAN §3'teki iki katmanlı doğrulamadan **farklı** bir problemdir:

| | Soru | Ne zaman |
|---|---|---|
| PLAN §3 (Katman A/B) | Pipeline doğru mu? | Bir kez, geliştirmede |
| **Bu spec** | *Bu koşunun* sonucuna güvenilir mi? | Her müşteri koşusunda |

Ali'nin gereksinimi: *"Yalancı sonuç asla istemem. End-to-end çalışacak, müşteri
sorunsuzca datalarının sonucunu güvenceli alabilmelidir."* ve *"Sorun varsa sorun var
densin. Patladıysa pipeline bilelim, oturup bal beklemeyeyim."*

## 2. Onaylanan kararlar

1. **Hedef kitle: karma (D).** Bakteriyel derinlik korunur, altyapı genele açık kalır.
2. **Kapı politikası: ikili (C).** Sonucu GEÇERSİZ kılan → `FAIL`, pipeline durur.
   Sonucu ŞÜPHELİ kılan → `WARN`, sonuç üretilir ve damgalanır.
3. **Eşikler veri, kod değil (C).** `organism_type` bir profil seçer. Config'ten
   ezilebilir; ezilen eşik rapora **açıkça** yazılır — sessiz gevşetme yok.
4. **FAIL çıktısı: teşhis raporu (B).** Biyolojik sonuç ÜRETİLMEZ; hangi kapının
   düştüğünü, ölçümü, eşiği, sorumlu örnekleri ve ne yapılacağını anlatan HTML üretilir.
5. **Eşleşmiş tasarım: tespit et, karar verme (C).** Veri eşleşmiş görünüp design
   kullanmıyorsa DUR ve sor. `paired: false` ile bilerek geçilebilir.

### Reddedilen seçenekler ve gerekçeleri

- **Tam otomatik design seçimi:** PLAN'ın kendi ilkesine ters (`organism_type`
  "varsayılanı yok, tahmin sessiz hataya yol açar"). Design tahmini daha tehlikeli,
  çünkü çıktı makul görünür.
- **Hepsi WARN:** Uyarı okunmayan şeydir. Müşteri tabloyu alır, uyarıyı atlar.
- **Hepsi FAIL:** Sınır durumlarda ürün sürekli "hayır" der, müşteri kaçar.
- **Damgalı ama üretilmiş sonuç:** Damgalı tablo yine tablodur. Birisi grafiği kopyalayıp
  sunuma koyar, damga kaybolur — korkulan senaryo tam bu yoldan gerçekleşir.
- **Sabit tek eşik seti:** "Herkes gelsin" kararıyla çelişir. Bakteriye göre ayarlanırsa
  insan örneği yalancı FAIL alır; insana göre gevşetilirse bakteride gerçek bozukluk kaçar.

## 3. Mimari

### 3.1 `rnaforge/gates.py` — çekirdek

```python
@dataclass(frozen=True)
class GateResult:
    name: str            # "alignment_rate"
    module: str          # "m04"
    status: str          # PASS | WARN | FAIL
    measured: float | None
    threshold: float | None
    message: str         # neden önemli
    remedy: str          # ne yapılmalı
    overridden: bool     # eşik config'ten gevşetildi mi
```

`GateFailure(Exception)` — bir veya daha fazla `FAIL` taşır; pipeline'ı durdurur.

### 3.2 Profiller — eşikler veri olarak

`rnaforge/profiles/prokaryote.yml`, `rnaforge/profiles/eukaryote.yml`.
`organism_type` hangisinin yükleneceğini seçer. Yeni profil eklemek kod değil, dosya
eklemektir.

Config'teki `quality:` bölümü eşikleri ezer. Ezilen her eşik `overridden=True` alır ve
güvence kartında **ayrı bir bölümde** listelenir.

**Ökaryot profili bilinçli olarak gevşektir** ve raporda "geniş toleranslı profil"
damgası taşır. Gerekçe: elde gerçek ökaryot doğrulaması yok. Uydurma bir sayıyı "eşik"
diye koymak kapı sisteminin tüm anlamını bozar. Ökaryot doğrulaması geldiğinde sıkılır.

### 3.3 Akış

Her modül işini bitirince kapılarını koşar → sonuçları `runs/<id>/quality/gates.json`'a
**ekler** (üzerine yazmaz; resume ile uyumlu). `FAIL` varsa `GateFailure` fırlatır,
sonraki modül **çalışmaz**.

### 3.4 Güvence kartı

`runs/<id>/quality/confidence_card.json` + HTML raporun en üstünde bir bölüm:
her kapı, ölçüm, eşik, karar, kullanılan profil, ezilen eşikler.

**Normal (PASS) raporda da bulunur** — müşteri neyin kontrol edildiğini görür.
Görünmeyen güvence, güvence değildir.

### 3.5 Rapor modülü iki kipli

- **Normal kip:** biyolojik sonuçlar + güvence kartı.
- **Teşhis kipi:** `FAIL` durumunda. DEG tablosu, volcano, PCA — hiçbiri üretilmez.
  Üretilen: düşen kapılar, ölçümler, sorumlu örnekler, düzeltme önerileri.

### 3.6 Metadata değişikliği

`subject` sütunu (opsiyonel) eklenir — eşleşmiş tasarımın dayanağı.

## 4. Kapı kataloğu

### 4.1 Tasarım kapıları — m01, veri okunmadan

Metadata'dan çıkar, saniyeler sürer, tek FASTQ okunmaz. Müşteri 6 saatlik koşunun
sonunda değil, **başlamadan** öğrenir. Hepsi `FAIL` — hepsi sonucu geçersiz kılar.

| Kapı | Ne yakalar | Durum |
|---|---|---|
| `design_rank` | confounded batch, tek seviyeli faktör | ✅ yazıldı |
| `replication` | replikasız koşul | ✅ yazıldı |
| `paired_declared` | subject tekrarlıyor ama design kullanmıyor | ⬜ bu spec |
| `reference_match` | referans dosyaları var mı, tutarlı mı | ✅ kısmen |

### 4.2 Veri kapıları — kendi modülleriyle gelir

| Modül | Kapı | Prokaryot | Ökaryot | Tip |
|---|---|---|---|---|
| m02 | `read_depth` | 1M | 10M | FAIL altında / WARN yakınında |
| m02 | `base_quality` | Q30 | Q30 | WARN |
| m03 | `survival_rate` | %50 | %50 | FAIL |
| m04 | `alignment_rate` | %70 | %50 | FAIL |
| m04 | `rrna_fraction` | %20 | %30 | WARN |
| m05 | `genes_detected` | profil | profil | WARN |
| m06 | `replicate_correlation` | r 0.85 | r 0.80 | WARN; çok düşükse FAIL |
| m06 | `sample_swap` | replika yanlış kümede | | FAIL |
| m06 | `dispersion_fit` | model uyumu | | WARN |

`sample_swap` kataloğun en değerli kapısıdır: örnek etiketlerinin karışması laboratuvarda
gerçekten olur ve başka hiçbir kontrol yakalamaz.

## 5. Kapsam

**Bu spec'in uygulama kapsamı:**
- `gates.py` çerçevesi (GateResult, GateFailure, koşum, gates.json)
- Profil yükleme + config ile ezme + override kaydı
- Tasarım kapıları (m01) — mevcut kontroller kapı sözleşmesine taşınır + `paired_declared`
- `subject` sütunu desteği ve eşleşmiş design doğrulaması
- Güvence kartı üretimi (JSON; HTML rapor m08'de)

**Kapsam dışı — kendi modülleriyle gelecek:**
Veri kapıları (§4.2). m02-m06 henüz yazılmadı; kapılar o modüllerle birlikte yazılır.
Katalog burada tanımlı durur, boş vaat olmaz.

## 6. Test stratejisi

- Her kapı için: geçen durum, WARN durumu, FAIL durumu.
- **Yanlış pozitif testleri zorunlu:** dengeli tasarım geçmeli, sağlıklı veri PASS almalı.
  Kapı sistemi yalancı FAIL üretirse müşteri ürüne güvenmez.
- Override akışı: gevşetilen eşik `overridden=True` ile kartta görünmeli.
- FAIL akışı uçtan uca: biyolojik çıktı üretilmediği **dosya sisteminde** doğrulanmalı.
