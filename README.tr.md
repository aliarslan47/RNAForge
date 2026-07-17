# RNAForge

Yeniden üretilebilir, modüler Bulk RNA-seq analiz pipeline'ı.

İngilizce sürüm: [README.md](README.md) · Referans doküman: [PLAN.md](PLAN.md) (v1.2)

## Durum

Erken geliştirme. Şu an hazır olan: girdi doğrulama ve platform tespiti (`m01`).

## Kurulum

```bash
conda env create -f envs/rnaforge-core.yml
conda activate rnaforge-core
pip install -e .
```

## Kullanım

```bash
rnaforge validate --config config/config.yaml --metadata samples.tsv --run-id demo
```

### Metadata formatı (TSV)

| Sütun | Zorunlu | Açıklama |
|---|---|---|
| `sample_id` | evet | Benzersiz örnek kimliği |
| `condition` | evet | Deney grubu; ≥2 seviye ve her seviyede ≥2 replika gerekir |
| `fastq_1` | evet | R1 yolu (veya single-end okumalar) |
| `fastq_2` | hayır | Paired-end için R2 yolu |
| `batch` | hayır | Batch/kovaryat; design formülü `batch` kullanıyorsa zorunlu |

## Temel tasarım kararları

- **`organism_type` zorunludur, varsayılanı yoktur** (`prokaryote` | `eukaryote`).
  Kantifikasyonu yönlendirir: prokaryotta genom hizalama + featureCounts, ökaryotta
  Salmon + tximport. İkisi de aynı gen × örnek count matrisinde buluşur.
- **Yalnızca Illumina (MVP).** ONT/PacBio girdileri tespit edilir ve sessizce yanlış
  yoldan işlenmek yerine net bir hatayla reddedilir.
- **Trimming bilinçli olarak naziktir.** Agresif kalite trimming ekspresyon tahminlerini
  bozar ([Williams et al. 2016](https://doi.org/10.1186/s12859-016-0956-2)); sapmayı
  engelleyen şey minimum uzunluk filtresidir.

## Geliştirme

```bash
conda run -n rnaforge-core python -m pytest -v
```

## Gizlilik

Müşteri verisi asla commit edilmez. `runs/`, `raw/` ve `references/` git tarafından yok sayılır.
