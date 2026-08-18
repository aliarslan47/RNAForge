"""m09 — GO annotation building (pure Python). Gene(locus_tag) -> GO set.

Kaynak önceliği (spec §2): GFF OTORİTE → GAF yalnız GFF-GO'suz genlere (tam+benzersiz
sembol eşleşmesi, belirsizlik ATILIR) → go-basic.obo ile ata-terimlere propagation.
Her GO kaydı kaynak-damgalı (GFF|GOA). Yalancı anotasyon üretmez (belirsizliği tahmin etmez).
"""
from __future__ import annotations

import gzip
import re
from pathlib import Path

# GFF go_process/go_function/go_component -> GO namespace kısaltması.
_GFF_ASPECT = {"go_process": "BP", "go_function": "MF", "go_component": "CC"}
# GAF sütun 9 (aspect) -> namespace; obo namespace tam adı -> kısaltma.
_GAF_ASPECT = {"P": "BP", "F": "MF", "C": "CC"}
_OBO_NS = {"biological_process": "BP", "molecular_function": "MF", "cellular_component": "CC"}


def _attrs(line: str) -> dict[str, str]:
    """GFF 9. sütun (col 8, 0-index) attribute'larını dict'e ayrıştır."""
    field = line.rstrip("\n").split("\t")[8]
    return dict(kv.split("=", 1) for kv in field.split(";") if "=" in kv)


def parse_gff_go(gff_path: Path):
    """CDS satırlarından GO anotasyonu.

    Returns:
        gene2go: dict[locus_tag, set[go_id]]  (yalnız Ontology_term'den; GO'suz gen yok)
        go_meta: dict[go_id, (namespace, name)]  (go_process/function/component'ten)
        gene_symbol: dict[locus_tag, symbol]  (gene= alanı)
    """
    gene2go: dict[str, set[str]] = {}
    go_meta: dict[str, tuple[str, str]] = {}
    gene_symbol: dict[str, str] = {}
    for line in Path(gff_path).read_text().splitlines():
        if not line or line.startswith("#") or "\tCDS\t" not in line:
            continue
        d = _attrs(line)
        lt = d.get("locus_tag")
        if not lt:
            continue
        if lt not in gene_symbol and d.get("gene"):
            gene_symbol[lt] = d["gene"]
        # go_process/function/component: "ad|<id>||kanıt" parçaları -> namespace + ad
        for key, ns in _GFF_ASPECT.items():
            if key not in d:
                continue
            for entry in d[key].split(","):
                parts = entry.split("|")
                if len(parts) < 2 or not parts[1]:
                    continue
                go_id = f"GO:{parts[1]}"
                go_meta.setdefault(go_id, (ns, parts[0]))
        # Ontology_term: genin GO id kümesi (otorite id listesi)
        ids = d.get("Ontology_term")
        if ids:
            gene2go.setdefault(lt, set()).update(
                gid for gid in ids.split(",") if gid.startswith("GO:")
            )
    return gene2go, go_meta, gene_symbol


def _open_text(path: Path):
    """gz veya düz metin aç (transkriptom FASTA .fa/.fa.gz olabilir)."""
    path = Path(path)
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def parse_transcriptome_symbols(fasta_path: Path) -> dict[str, str]:
    """Ökaryot: transkriptom FASTA başlıklarından gen ID -> sembol.

    Ensembl cDNA başlığı: '>ENST... gene:ENSG00000103196.12 ... gene_symbol:CRISPLD2 ...'.
    Anahtar = gen ID (versiyonlu — tx2gene'in gen sütunu ve DE tablosuyla eşleşir).
    gene_symbol yoksa gen atlanır (map'e girmez). İlk görülen kazanır."""
    out: dict[str, str] = {}
    gene_re = re.compile(r"\bgene:(\S+)")
    sym_re = re.compile(r"\bgene_symbol:(\S+)")
    with _open_text(fasta_path) as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            gm = gene_re.search(line)
            sm = sym_re.search(line)
            if gm and sm and gm.group(1) not in out:
                out[gm.group(1)] = sm.group(1)
    return out


def parse_annotation_symbols(gff_path: Path | None,
                             transcriptome_fasta: Path | None):
    """Sembol/GO kaynağı dispatcher'ı (prokaryot GFF vs ökaryot transkriptom).

    Returns (gene2go, meta, gene_symbol) — parse_gff_go ile aynı şekil.
    - gff_path verilmişse: parse_gff_go (Ontology_term GO + sembol).
    - yoksa transcriptome_fasta: boş gene2go/meta + transkriptom sembolleri (ökaryot;
      GO'yu GAF dolduracak).
    - ikisi de yoksa: yüksek sesle hata (sessiz boş anotasyon yok)."""
    if gff_path is not None:
        return parse_gff_go(gff_path)
    if transcriptome_fasta is not None:
        return {}, {}, parse_transcriptome_symbols(transcriptome_fasta)
    raise ValueError(
        "anotasyon için annotation_gff (prokaryot) veya transcriptome_fasta (ökaryot) "
        "gerekli; ikisi de yok.")


def parse_obo(obo_path: Path) -> dict[str, dict]:
    """go-basic.obo -> {go_id: {name, namespace(BP/MF/CC), parents:set, obsolete:bool}}.

    parents = is_a + relationship: part_of hedefleri (propagation kenarları)."""
    terms: dict[str, dict] = {}
    cur: dict | None = None
    in_term = False
    for line in Path(obo_path).read_text().splitlines():
        line = line.rstrip("\n")
        if line == "[Term]":
            in_term = True
            cur = {"name": "", "namespace": "", "parents": set(), "obsolete": False}
            continue
        if line.startswith("["):  # başka bir stanza (Typedef vs) -> term dışı
            in_term = False
            cur = None
            continue
        if not in_term or cur is None or ":" not in line:
            continue
        key, _, val = line.partition(":")
        val = val.strip()
        if key == "id" and val.startswith("GO:"):
            terms[val] = cur
        elif key == "name":
            cur["name"] = val
        elif key == "namespace":
            cur["namespace"] = _OBO_NS.get(val, val)
        elif key == "is_a":
            cur["parents"].add(val.split("!")[0].strip())
        elif key == "relationship" and val.startswith("part_of "):
            cur["parents"].add(val.split()[1])
        elif key == "is_obsolete" and val == "true":
            cur["obsolete"] = True
    return terms


def _ancestors(go_id: str, obo: dict, cache: dict[str, set[str]]) -> set[str]:
    """go_id'nin tüm ataları (kendisi hariç), geçişli + döngü-korumalı."""
    if go_id in cache:
        return cache[go_id]
    cache[go_id] = set()  # döngü koruması: hesaplanırken boş dur
    result: set[str] = set()
    for parent in obo.get(go_id, {}).get("parents", ()):
        if parent not in obo:
            continue
        result.add(parent)
        result |= _ancestors(parent, obo, cache)
    cache[go_id] = result
    return result


def propagate(gene2go: dict[str, set[str]], obo: dict) -> dict[str, set[str]]:
    """Her genin GO setini ata-terimlerle genişlet. Obsolete terimler dışlanır."""
    cache: dict[str, set[str]] = {}
    out: dict[str, set[str]] = {}
    for gene, gos in gene2go.items():
        expanded: set[str] = set()
        for go_id in gos:
            expanded.add(go_id)
            expanded |= _ancestors(go_id, obo, cache)
        out[gene] = {g for g in expanded if not obo.get(g, {}).get("obsolete", False)}
    return out


def _symbol_to_locus(gene_symbol: dict[str, str]) -> dict[str, str]:
    """symbol -> locus_tag ters harita. Aynı sembol ≥2 locus'a giderse ATILIR (belirsiz)."""
    counts: dict[str, list[str]] = {}
    for lt, sym in gene_symbol.items():
        counts.setdefault(sym, []).append(lt)
    return {sym: lts[0] for sym, lts in counts.items() if len(lts) == 1}


def fill_from_gaf(gene2go: dict[str, set[str]], gene_symbol: dict[str, str],
                  gaf_path: Path, go_meta: dict[str, tuple[str, str]]):
    """GAF'ı yalnız GFF-GO'suz genlere ekle. Belirsiz eşleşme atılır.

    Returns:
        additions: dict[locus_tag, set[go_id]]  (GAF'tan eklenenler)
        sources: dict[(locus_tag, go_id), "GOA"]
    """
    sym2lt = _symbol_to_locus(gene_symbol)
    # GAF'ı sembole göre topla: {symbol: {"gos": {(go,aspect)}, "uniprots": set}}
    by_symbol: dict[str, dict] = {}
    for line in Path(gaf_path).read_text().splitlines():
        if not line or line.startswith("!"):
            continue
        cols = line.split("\t")
        if len(cols) < 9:
            continue
        uniprot, symbol, go_id, aspect = cols[1], cols[2], cols[4], cols[8]
        if not go_id.startswith("GO:") or symbol not in sym2lt:
            continue
        rec = by_symbol.setdefault(symbol, {"gos": set(), "uniprots": set()})
        rec["gos"].add((go_id, aspect))
        rec["uniprots"].add(uniprot)

    additions: dict[str, set[str]] = {}
    sources: dict[tuple[str, str], str] = {}
    for symbol, rec in by_symbol.items():
        if len(rec["uniprots"]) != 1:  # aynı sembol ≥2 farklı proteine -> belirsiz, atla
            continue
        lt = sym2lt[symbol]
        if lt in gene2go and gene2go[lt]:  # GFF otorite: zaten GO'su varsa dokunma
            continue
        for go_id, aspect in rec["gos"]:
            additions.setdefault(lt, set()).add(go_id)
            sources[(lt, go_id)] = "GOA"
            go_meta.setdefault(go_id, (_GAF_ASPECT.get(aspect, aspect), ""))
    return additions, sources


def build_gene2go(gff_path: Path | None, obo: dict, gaf_path: Path | None = None,
                  transcriptome_fasta: Path | None = None, log=None):
    """Tam annotation birleştirme: GFF otorite → GAF doldurma → propagation.
    Ökaryot: gff_path=None + transcriptome_fasta verilir → sembol köprüsü transkriptomdan,
    gene2go boş başlar, GAF birincil GO kaynağı olur.

    Returns:
        gene2go: dict[locus_tag, set[go_id]]  (propagate edilmiş)
        go_meta: dict[go_id, (namespace, name)]  (obo kanonik, GFF/GAF tamamlar)
        direct: dict[locus_tag, set[go_id]]  (propagation ÖNCESİ, denetim izi için)
        sources: dict[(locus_tag, go_id), "GFF"|"GOA"]
        stats: dict
        gene_symbol: dict[locus_tag, symbol]  (GFF gene= alanı; GFF zaten burada parse
            edildiği için çağıran ikinci kez parse etmesin diye döndürülür)
    """
    gene2go, gff_meta, gene_symbol = parse_annotation_symbols(gff_path, transcriptome_fasta)
    sources: dict[tuple[str, str], str] = {}
    for lt, gos in gene2go.items():
        for g in gos:
            sources[(lt, g)] = "GFF"
    n_gff = len(gene2go)

    n_goa = 0
    if gaf_path is not None:
        additions, gaf_sources = fill_from_gaf(gene2go, gene_symbol, gaf_path, gff_meta)
        for lt, gos in additions.items():
            gene2go.setdefault(lt, set()).update(gos)
        sources.update(gaf_sources)
        n_goa = len(additions)
    elif log is not None:
        log(f"m09: GAF verilmedi — yalnız GFF GO + propagation kullanılıyor")

    direct = {lt: set(gos) for lt, gos in gene2go.items()}
    propagated = propagate(gene2go, obo)

    # go_meta: obo kanonik ad/namespace; obo'da yoksa GFF/GAF'tan gelen.
    go_meta: dict[str, tuple[str, str]] = dict(gff_meta)
    for gid, info in obo.items():
        if info.get("namespace") and info.get("name"):
            go_meta[gid] = (info["namespace"], info["name"])

    stats = {
        "n_annotated": len(propagated),
        "n_gff": n_gff,
        "n_goa": n_goa,
        "n_terms": len({g for gos in propagated.values() for g in gos}),
    }
    return propagated, go_meta, direct, sources, stats, gene_symbol
