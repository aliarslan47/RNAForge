# m09 GO Enrichment (ORA) Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** DESeq2 sonuçlarından artan/azalan DEG'ler için GO over-representation analizi (ORA) yürüten, tablo+figür üreten ve m08 rapora GO bölümü ekleyen m09 modülünü kurmak.

**Architecture:** Saf Python annotation birleştirme (GFF otorite + GAF doldurma + obo propagation) → stdlib hipergeometrik ORA → R/ggplot2 dot-plot (m07 deseni) → m09 orkestratör (gate yok, verdict m06'dan taşınır) → CLI `enrich` → m08 rapor bölümü. Zincir: m06→m07→m09→m08.

**Tech Stack:** Python 3.11 stdlib (`math.comb`), R/ggplot2 (`rnaforge-de` env), pytest.

## Global Constraints

- YENİ veri-kapısı YOK; verdict m06/m07'den değişmeden taşınır (m09 dokunmaz).
- Belirsiz annotation eşleşmesi TAHMİN EDİLMEZ, atılır; her GO kaydı `source ∈ {GFF, GOA}` damgalı.
- Sessiz hata yasak: eksik obo/gaf → yüksek sesle hata + talimat; GAF yoksa log'a "yalnız GFF" yaz.
- Uydurma biyolojik yorum yok; rapor yalnız gerçek sayılar.
- Saf Python (ORA) + R (figür); stdlib dışı ağır bağımlılık yok.
- DEG eşikleri TEK kaynak: `config.de.fdr_threshold` / `config.de.log2fc_threshold`.
- Boş set / anlamlı-terim-yok → başlıklı boş çıktı + boş-durum paneli; ÇÖKME YOK.
- Test komutu (repo kökünden): `conda run -n rnaforge-core --cwd /home/ali/rnaforge-pipeline python -m pytest -q`
- Modül deseni m07 ile birebir: saf yardımcı + R runner + orkestratör + CLI + TDD.

---

### Task 1: Config — `EnrichmentConfig`

**Files:**
- Modify: `rnaforge/config.py`
- Test: `tests/test_config.py`

**Interfaces:**
- Produces: `Config.enrichment: EnrichmentConfig` alanı; `EnrichmentConfig(min_term_size:int=3, top_n:int=15, obo:Path|None, gaf:Path|None)`.

- [ ] **Step 1: Failing test** — `tests/test_config.py`'e ekle:
```python
def test_enrichment_config_defaults(tmp_path):
    cfg = load_config_from_dict({..minimal valid config..})
    assert cfg.enrichment.min_term_size == 3
    assert cfg.enrichment.top_n == 15

def test_enrichment_config_parsed(...):
    # enrichment: {min_term_size: 5, top_n: 10, obo: ref/go.obo, gaf: ref/e.gaf}
    assert cfg.enrichment.min_term_size == 5
    assert cfg.enrichment.obo == Path("ref/go.obo")

def test_unknown_enrichment_key_rejected(...):
    # enrichment: {min_term_sizeX: 5} -> ConfigError
```
- [ ] **Step 2:** testleri koştur, FAIL (AttributeError: enrichment).
- [ ] **Step 3:** `config.py`: `EnrichmentConfig` frozen dataclass ekle; `Config`'e `enrichment` alanı; `"enrichment"` → `KNOWN_TOP_LEVEL_KEYS`; parse bloğu (`enrichment_raw = raw.get("enrichment", {})`, bilinmeyen alt-anahtar `_reject_unknown` deseniyle reddedilir — diğer alt-config'lerin yaptığı gibi kontrol et ve aynısını uygula).
- [ ] **Step 4:** testler PASS.
- [ ] **Step 5:** commit `feat(m09): config EnrichmentConfig + KNOWN_TOP_LEVEL_KEYS`.

---

### Task 2: GFF GO parser

**Files:**
- Create: `rnaforge/go_annotation.py`
- Test: `tests/test_go_annotation.py`

**Interfaces:**
- Produces: `parse_gff_go(gff_path:Path) -> tuple[dict[str,set[str]], dict[str,tuple[str,str]], dict[str,str]]`
  → `(gene2go_direct, go_meta{go_id:(namespace,name)}, gene_symbol{locus_tag:symbol})`.
  namespace ∈ {"BP","MF","CC"}.

- [ ] **Step 1: Failing test** — küçük GFF fixture (2 CDS satırı, biri Ontology_term+go_process'li, biri GO'suz):
```python
def test_parse_gff_go_extracts_ids_namespace_name(tmp_path):
    gff = tmp_path/"g.gff"; gff.write_text(GFF_2CDS)
    g2go, meta, sym = parse_gff_go(gff)
    assert "GO:0008652" in g2go["BW25113_RS00010"]
    assert meta["GO:0008652"] == ("BP", "amino acid biosynthetic process")
    assert sym["BW25113_RS00010"] == "thrA"
    assert "BW25113_RS00005" not in g2go  # GO'suz gen map'te yok
```
- [ ] **Step 2:** FAIL (modül yok).
- [ ] **Step 3:** `parse_gff_go`: CDS satırlarını ayrıştır (mevcut `figures.py:_parse_attrs` deseni gibi `;`/`=` attr parse). `Ontology_term` → id seti. `go_process|go_function|go_component` → her `ad|GOID_no_prefix||kanıt` parçasından `GO:<id>` + namespace (process→BP, function→MF, component→CC) + ad. `gene=` → symbol.
- [ ] **Step 4:** PASS.
- [ ] **Step 5:** commit `feat(m09): GFF GO annotation parser`.

---

### Task 3: obo parser + propagation

**Files:**
- Modify: `rnaforge/go_annotation.py`
- Test: `tests/test_go_annotation.py`

**Interfaces:**
- Produces: `parse_obo(obo_path:Path) -> dict[str,dict]` (`{go_id:{"name","namespace","parents":set,"obsolete":bool}}`);
  `propagate(gene2go:dict[str,set], obo:dict) -> dict[str,set]` (ata-terimlerle genişletilmiş).

- [ ] **Step 1: Failing test** — küçük obo fixture (3 term zinciri C→B→A `is_a`, biri part_of):
```python
def test_parse_obo_parents_and_obsolete(tmp_path): ...
    assert obo["GO:0000003"]["parents"] == {"GO:0000002"}
    assert obo["GO:0000009"]["obsolete"] is True
def test_propagate_adds_ancestors():
    g2go = {"geneX": {"GO:0000003"}}
    out = propagate(g2go, obo)
    assert out["geneX"] == {"GO:0000003","GO:0000002","GO:0000001"}
def test_propagate_cycle_safe(): ...  # yapay döngüde takılmaz
```
- [ ] **Step 2:** FAIL.
- [ ] **Step 3:** `parse_obo`: `[Term]` bloklarını ayrıştır (`id:`, `name:`, `namespace:` → BP/MF/CC eşle, `is_a: GO:x ! ...`, `relationship: part_of GO:x`, `is_obsolete: true`). `propagate`: her gen için BFS/DFS ile parents geçişli kapanışı (visited set → döngü koruması), obsolete terimleri atla.
- [ ] **Step 4:** PASS.
- [ ] **Step 5:** commit `feat(m09): obo parser + ancestor propagation`.

---

### Task 4: GAF doldurma + `build_gene2go`

**Files:**
- Modify: `rnaforge/go_annotation.py`
- Test: `tests/test_go_annotation.py`

**Interfaces:**
- Produces: `fill_from_gaf(gene2go, gene_symbol, gaf_path, obo) -> tuple[dict[str,set], dict[str,str]]`
  (yalnız GFF'te GO'su OLMAYAN genlere, tam+benzersiz sembol eşleşmesiyle ekler; `source` haritası döndürür);
  `build_gene2go(gff, obo, gaf=None, log=None) -> tuple[gene2go(propagated), go_meta, sources, stats]`.

- [ ] **Step 1: Failing test:**
```python
def test_gaf_fills_only_ungapped_unique_symbol():
    # gene A: GFF'te GO var -> GAF eklenmez; gene B: GFF'te yok, GAF'ta benzersiz sembol -> eklenir
    # gene C: sembol GFF'te 2 locus'a karşılık -> atılır (tahmin yok)
def test_build_gene2go_source_stamped(): assert sources[("geneB","GO:x")] == "GOA"
```
- [ ] **Step 2:** FAIL.
- [ ] **Step 3:** `fill_from_gaf`: GAF satırlarını oku (tab; sütun 3 `db_object_symbol`, 5 `go_id`, 9 aspect P/F/C). Sembol→locus_tag ters haritasını `gene_symbol`'den kur; **çoklu locus'lu sembolleri diskalifiye et**. GAF'ta aynı sembol ≥2 farklı UniProt (sütun 2) → atla. Yalnız GFF-GO'suz genlere ekle. `build_gene2go`: parse_gff_go → (gaf varsa) fill_from_gaf → propagate; go_meta'yı obo adı/namespace ile tamamla; stats (n_annotated, source sayıları).
- [ ] **Step 4:** PASS.
- [ ] **Step 5:** commit `feat(m09): GAF fill (safe unique-symbol) + build_gene2go`.

---

### Task 5: Hipergeometrik + BH

**Files:**
- Create: `rnaforge/enrichment.py`
- Test: `tests/test_enrichment.py`

**Interfaces:**
- Produces: `hypergeometric_pvalue(k,n,K,N) -> float` (over-representation, üst kuyruk);
  `bh_fdr(pvalues:list[float]) -> list[float]`.

- [ ] **Step 1: Failing test** — elle hesaplanmış küçük değerler:
```python
def test_hypergeometric_known_value():
    # N=10,K=5,n=4,k=4 -> p = C(5,4)C(5,0)/C(10,4) = 5/210
    assert hypergeometric_pvalue(4,4,5,10) == pytest.approx(5/210)
def test_hypergeometric_upper_tail_sums():  # k=3: i=3,4
    assert hypergeometric_pvalue(3,4,5,10) == pytest.approx((C(5,3)*C(5,1)+C(5,4)*C(5,0))/C(10,4))
def test_bh_monotone_and_bounded():
    adj = bh_fdr([0.01,0.02,0.03,0.04]); assert all(0<=a<=1 for a in adj) and adj==sorted(adj)
```
- [ ] **Step 2:** FAIL.
- [ ] **Step 3:** `hypergeometric_pvalue`: `sum(comb(K,i)*comb(N-K,n-i) for i in range(k, min(K,n)+1))/comb(N,n)`. `bh_fdr`: sırala, `p*m/rank`, ardışık min (cumulative min from largest), 1'e kırp, orijinal sıraya döndür.
- [ ] **Step 4:** PASS.
- [ ] **Step 5:** commit `feat(m09): hypergeometric ORA + BH FDR (stdlib)`.

---

### Task 6: ORA runner (TSV)

**Files:**
- Modify: `rnaforge/enrichment.py`
- Test: `tests/test_enrichment.py`

**Interfaces:**
- Produces: `deg_sets(deseq_tsv, fdr, lfc) -> tuple[list[str], list[str]]` (up, down locus_tag);
  `run_ora(gene_set, background, gene2go, go_meta, gene_symbol, min_term_size) -> list[dict]` (terim satırları);
  `write_ora_tsv(rows, path)`.

- [ ] **Step 1: Failing test:**
```python
def test_deg_sets_split_by_direction(): up,down = deg_sets(tsv,0.05,1.0); assert "geneUp" in up
def test_run_ora_enriched_term_significant():
    # kurgulanmış: bir terim set'te aşırı temsil -> düşük p, fold>1, genes doğru
def test_run_ora_respects_min_term_size(): # bg_count<min -> terim yok
def test_run_ora_bh_per_namespace(): # padj her namespace içinde
def test_write_ora_tsv_empty_set_header_only(): # boş set -> yalnız başlık, çökme yok
```
- [ ] **Step 2:** FAIL.
- [ ] **Step 3:** `deg_sets`: deseq2_results.tsv oku, padj/lfc filtre. `run_ora`: background = gene2go'su olan genler; her namespace için terim→bg genleri say; min_term_size filtre; her terim için k (set∩term), hypergeometric; namespace başına bh_fdr; satır dict'leri (go_id, namespace, term, study_count, study_n, bg_count, bg_n, expected, fold_enrichment, p_value, p_adj, genes=set∩term sembolleri). padj artan sırala. `write_ora_tsv` sabit başlık.
- [ ] **Step 4:** PASS.
- [ ] **Step 5:** commit `feat(m09): ORA runner + deg_sets + TSV writer`.

---

### Task 7: enrichment.R figürü + runner

**Files:**
- Create: `rnaforge/scripts/enrichment.R`
- Modify: `rnaforge/enrichment.py`
- Test: `tests/test_enrichment.py` (env-gated integration)

**Interfaces:**
- Produces: `run_enrichment_r(up_tsv, down_tsv, out_dir, top_n) -> str` (rnaforge-de env, PNG+SVG);
  `build_enrichment_manifest(out_dir) -> dict` (m07 build_manifest deseni; id: enrichment_up/down).

- [ ] **Step 1: Failing test** (skip if env yok): küçük up/down TSV → `run_enrichment_r` → `enrichment_up.png` var.
- [ ] **Step 2:** FAIL/skip.
- [ ] **Step 3:** `enrichment.R`: argümanları oku (up_tsv, down_tsv, out_dir, top_n); her set için `p_adj<0.05` süz, top_n al, `ggplot(aes(fold_enrichment, reorder(term,fold_enrichment), size=study_count, color=p_adj)) + geom_point() + facet_grid(namespace~., scales="free_y", space="free")`; boş → boş-durum paneli (`annotate("text",...)`). `ggsave` PNG 300dpi + SVG. `run_enrichment_r`: `conda run -n rnaforge-de Rscript ...`, stdout döndür. `build_enrichment_manifest`: figures.py `build_manifest` mantığını yeniden kullan (id başına png/svg).
- [ ] **Step 4:** PASS/skip; elle bir kez env'de doğrula.
- [ ] **Step 5:** commit `feat(m09): enrichment dot-plot (R/ggplot2) + manifest`.

---

### Task 8: m09 orkestratör + CLI

**Files:**
- Create: `rnaforge/modules/m09_enrichment.py`
- Modify: `rnaforge/cli.py`
- Test: `tests/test_m09_enrichment.py`

**Interfaces:**
- Consumes: `build_gene2go`, `deg_sets`, `run_ora`, `write_ora_tsv`, `run_enrichment_r`, `build_enrichment_manifest`, `RunState`.
- Produces: `run_enrichment(config, metadata_path, run_dir, force=False) -> dict`; CLI `rnaforge enrich`.

- [ ] **Step 1: Failing test** (R monkeypatch'li):
```python
def test_run_enrichment_requires_m06(tmp_path): # m06 done değil -> ValueError
def test_run_enrichment_no_gate_verdict_untouched(): # gates.json'a dokunmaz
def test_run_enrichment_writes_tsvs_and_stats(monkeypatch): # up/down tsv + enrichment_statistics.json
def test_run_enrichment_missing_obo_loud(): # obo yok -> FileNotFoundError net mesaj
def test_run_enrichment_resume(): # ikinci çağrı resumed=True
```
- [ ] **Step 2:** FAIL.
- [ ] **Step 3:** `m09_enrichment.py` m07 deseninden: dizinler, RunState, resume, m06 ön koşul ValueError, obo yoksa gürültülü hata; `build_gene2go` → `deg_sets` → `run_ora`(up/down) → `write_ora_tsv` → `run_enrichment_r` → `build_enrichment_manifest` → stats JSON (`n_terms_up/down`, `n_sig_up/down`, `background_size`, `n_annotated`, source sayıları); gate YOK. `cli.py`: `enrich` subparser (config/metadata/run-id/force) + `_cmd_enrich` (run_enrichment + güvence kartı yaz, verdict taşınır) + `main` dispatch.
- [ ] **Step 4:** PASS.
- [ ] **Step 5:** commit `feat(m09): orchestrator + rnaforge enrich CLI`.

---

### Task 9: m08 rapor GO bölümü

**Files:**
- Modify: `rnaforge/report_html.py`, `rnaforge/modules/m08_report.py`
- Test: `tests/test_report_html.py`

**Interfaces:**
- Consumes: `enrichment_up.tsv`, `enrichment_down.tsv`, `enrichment_statistics.json`, enrichment manifest.
- Produces: rapor HTML'inde "Fonksiyonel Zenginleştirme (GO)" bölümü.

- [ ] **Step 1: Failing test:**
```python
def test_report_includes_go_section_when_present(tmp_path): # enrichment/ dolu -> başlık + tablo + <img> gömülü
def test_report_go_section_absent_note_when_missing(): # enrichment yok -> "çalıştırılmadı" notu, kırılmaz
```
- [ ] **Step 2:** FAIL.
- [ ] **Step 3:** `load_report_inputs`'a enrichment TSV/stats/manifest okuma ekle (yoksa None). Yeni `section_enrichment(...)` (DE'den sonra): UP/DOWN namespace başına top-N anlamlı terim `<table>` + gömülü dot-plot (`_img_data_uri` mevcut yardımcı) + çift dilli caption/intro. `render_report`'a bölümü ekle (m08 organizma-agnostik kalır). `m08_report.py`: ön koşulu m07→**m07+m09**; m09 yoksa dürüst not (kırılmaz).
- [ ] **Step 4:** PASS. Tüm suit yeşil.
- [ ] **Step 5:** commit `feat(m09): m08 raporda GO enrichment bölümü`.

---

### Task 10: Referans hazırlık + GSE300731 canlı smoke

**Files:**
- Modify: `.gitignore` (gerekliyse `references/go/`), `docs/` prep notu
- Create: `references/go/go-basic.obo`, `references/ecoli_bw25113/ecoli.gaf` (gitignore'lu, indirilir)

- [ ] **Step 1:** obo indir (`purl.obolibrary.org/obo/go/go-basic.obo`) + E. coli GAF (EBI-GOA); **formatı doğrula** (sütun/başlık kontrolü — körlemesine güvenme). `.gitignore`'da olduklarını teyit et.
- [ ] **Step 2:** GSE300731 config'ine `enrichment` bloğu (obo/gaf yolları) ekle (runs dışı; scratchpad config veya repo örnek config).
- [ ] **Step 3:** `rnaforge enrich --run-id 20260803_143036_GSE300731` (aynı run_dir) koştur.
- [ ] **Step 4:** Biyolojik akıl kontrolü: DOWN'da gad*/hde* asit-direnç GO'ları; UP'ta zarf-stres/kapsül GO'ları anlamlı mı? Figürler + `enrichment_*.tsv` incele.
- [ ] **Step 5:** `rnaforge report` yeniden koştur → GO bölümü + figürler gömülü; verdict SUSPECT (değişmez). Doğrulama notu.

---

## Self-Review

**Spec coverage:** §2 annotation→Task 2-4; §3 ORA→Task 5-6; §4 figür→Task 7; §5 orkestrasyon→Task 8; §6 CLI→Task 8; §7 config→Task 1; §8 rapor→Task 9; §10 doğrulama→her task testleri + Task 10 smoke. Tümü kapsanıyor.

**Placeholder scan:** Kod adımları somut imza + algoritma içeriyor; "uygun hata ekle" tarzı yer tutucu yok.

**Type consistency:** `build_gene2go`→`run_ora`/`run_enrichment`; `parse_gff_go` üçlü dönüşü Task 4 tüketir; `deg_sets`/`run_ora`/`write_ora_tsv` Task 8 tüketir; manifest id'leri (`enrichment_up/down`) Task 7↔9 tutarlı.
