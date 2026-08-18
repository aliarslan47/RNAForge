# Ökaryot zenginleştirme anotasyonu (GFF'siz GO/KEGG) — Tasarım

**Tarih:** 2026-08-18
**Durum:** Onaylı (brainstorming), implementasyona hazır.
**Kapsam:** m09 (GO ORA) + m10 (KEGG) anotasyon katmanını ökaryotta (GFF yok, transkriptom+tx2gene)
çalışacak şekilde adapte et. m11 (GSEA) m09/m10 kurucularını yeniden kullandığından otomatik gelir.

## 1. Problem

m09/m10 anotasyonu **GFF-tabanlı**: `parse_gff_go(gff)` → gen→GO (Ontology_term) + gen→sembol;
GAF ve KEGG bu sembol köprüsüne bağlanıyor. Ökaryot (insan) transkriptom-only koşuda GFF yok, ve
insan için `Ontology_term`'lü GFF pratikte yok → GO/KEGG koşamıyor (m10 `gff = config.reference.annotation_gff`).

## 2. Fikir

Ökaryotta **sembol köprüsünü transkriptom FASTA başlıklarından** kur. Ensembl cDNA başlıkları
zaten `gene:<ID> ... gene_symbol:<SYM>` taşıyor (bizim gen ID'miz = tx2gene'in gen sütunu =
`gene:` alanı, versiyonlu, DE tablosuyla birebir). Böylece:
- `gene2go` **boş** başlar (GFF Ontology_term yok),
- mevcut **`fill_from_gaf`** (sembol→gen) GO'yu doldurur (GAF = tek kaynak),
- mevcut **KEGG-by-sembol** makinesi pathway'i doldurur.
GFF gerekmez; mevcut GAF/KEGG/propagation/ORA kodu DEĞİŞMEZ (yalnız sembol kaynağı değişir).

## 3. Değişiklikler

### 3.1 `rnaforge/go_annotation.py`
- **Yeni** `parse_transcriptome_symbols(fasta_path: Path) -> dict[str, str]`:
  FASTA başlıklarından `gene:<ID>` → `gene_symbol:<SYM>` (gz-aware; ilk görülen kazanır).
  Anahtar = gen ID (versiyonlu), DE gen sütunuyla eşleşir.
- **Yeni** dispatcher `parse_annotation_symbols(gff_path, transcriptome_fasta) -> (gene2go, meta, gene_symbol)`:
  - `gff_path` verilmişse → `parse_gff_go(gff_path)` (prokaryot, aynen)
  - değilse `transcriptome_fasta` verilmişse → `({}, {}, parse_transcriptome_symbols(...))` (ökaryot)
  - ikisi de yoksa → yüksek sesle `ValueError` (sessiz boş anotasyon yok).
- `build_gene2go(gff_path, obo, gaf_path=None, transcriptome_fasta=None, log=None)`: ilk satır
  `parse_gff_go(gff_path)` → `parse_annotation_symbols(gff_path, transcriptome_fasta)`. Gerisi aynı.
  (Ökaryotta gene2go boş → n_gff=0; GAF her şeyi doldurur; propagation aynı.)

### 3.2 `rnaforge/kegg_annotation.py`
- `build_gene2pathway(gff_path, links, names, genelist, transcriptome_fasta=None)`: içteki
  `parse_gff_go(gff_path)` → `parse_annotation_symbols(gff_path, transcriptome_fasta)`.

### 3.3 `rnaforge/modules/m09_enrichment.py` + `m10_kegg.py`
- Artık GFF'i zorunlu tutma: `gff = config.reference.annotation_gff` (None olabilir) +
  `transcriptome_fasta = config.reference.transcriptome_fasta` ikisini de kurucuya geçir.
- İkisi de yoksa kurucu zaten yüksek sesle hata verir.

## 4. Referanslar (indirilecek)

- **GO GAF (insan):** EBI GOA `goa_human.gaf` (`references/go/goa_human.gaf`, gz açılmış — `fill_from_gaf`
  düz metin okur). config `enrichment.gaf`. FTP: ftp.ebi.ac.uk/pub/databases/GO/goa/HUMAN/goa_human.gaf.gz
- **KEGG hsa:** `references/kegg/hsa/{gene_list,pathway_links,pathway_names}.tsv` (KEGG REST list/hsa +
  link/pathway/hsa + list/pathway/hsa; sty deseni). config `enrichment.kegg_organism: hsa`.
- `go-basic.obo` zaten var.

## 5. Doğrulama (airway_dex koşusu)

Mevcut `runs/20260818_110541_airway_dex` üzerinde `rnaforge enrich/kegg/gsea/semantic/report`.
Beklenti (glukokortikoid biyolojisi): GO/KEGG'de steroid/glukokortikoid yanıtı, inflamasyon-baskılanması,
metabolik yeniden programlama temaları. Konkordans yön düzeyinde (yayın DEG'iyle tutarlı gen setleri).

## 6. Test (modül deseni)

- `parse_transcriptome_symbols`: fixture FASTA (gene:+gene_symbol: başlık) → doğru map; gz de.
- `parse_annotation_symbols`: gff verilince parse_gff_go yolu; yalnız transcriptome verilince boş
  gene2go + semboller; ikisi de yoksa ValueError.
- `build_gene2go(transcriptome_fasta=..., gaf=...)`: gene2go GAF'tan dolar (gff=None).
- `build_gene2pathway(transcriptome_fasta=...)`: KEGG sembol→ENSG join.
- m09/m10: GFF'siz (ökaryot config) çağrı çökmez, ORA çıktısı üretir (küçük fixture, GAF/KEGG monkeypatch/mini).
- Regresyon: prokaryot (GFF) yolu değişmedi (mevcut testler yeşil).

## 7. Karar özeti

- **Sembol köprüsü (GAF/KEGG-by-symbol reuse), doğrudan gene2go-TSV DEĞİL:** mevcut makineyi aynen
  kullanır, ikinci bir anotasyon yolu bakımı gerektirmez. Belirsiz sembol ATILIR (mevcut `_symbol_to_locus`,
  yalancı yok). Transkriptom zaten indirilmiş → ek sembol kaynağı bedava.
- **GAF ökaryotta tek GO kaynağı** (GFF Ontology_term yok) — Typhi'de GFF otorite + GAF doldurma; ökaryotta
  otorite yok, GAF birincil. Damga: eukaryote profili zaten permissive/SUSPECT.
