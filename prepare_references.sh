#!/usr/bin/env bash
# RNAForge referans-veri hazırlığı (tek seferlik, git-ignore'lu references/ altına).
#
# Kullanım:
#   bash prepare_references.sh --kegg-org eco --string-taxid 511145 \
#        --goa-url https://ftp.ebi.ac.uk/pub/databases/GO/goa/proteomes/18.E_coli_MG1655.goa
#
# Argümanlar (hepsi opsiyonel — verilmeyen blok atlanır, yüksek sesle bildirilir):
#   --kegg-org <kod>       KEGG organizma kodu (eco/hsa/mmu/sty…) → references/kegg/<kod>/
#   --string-taxid <id>    STRING taxid (511145=E.coli K-12) → references/string/<id>/
#   --goa-url <url>        Organizma GO anotasyonu (GAF) kaynağı → references/go/<organizma>.gaf
#   --goa-out <yol>        GAF çıktı yolu (varsayılan references/go/organism.gaf)
#   --skip-obo             go-basic.obo indirmeyi atla (zaten varsa)
#
# NOT (engellenen kaynak): QuickGO indirmesi bazı ağlarda BLOKLU. Bu script GAF için
# doğrudan QuickGO KULLANMAZ; EBI-GOA FTP proteome dosyasını (--goa-url) kullanır — bu
# engel için belgeli fallback yoludur (DURUM 2026-08).
set -euo pipefail
cd "$(dirname "$0")"

KEGG_ORG=""; STRING_TAXID=""; GOA_URL=""; GOA_OUT="references/go/organism.gaf"; SKIP_OBO=0
while [ $# -gt 0 ]; do
  case "$1" in
    --kegg-org) KEGG_ORG="$2"; shift 2;;
    --string-taxid) STRING_TAXID="$2"; shift 2;;
    --goa-url) GOA_URL="$2"; shift 2;;
    --goa-out) GOA_OUT="$2"; shift 2;;
    --skip-obo) SKIP_OBO=1; shift;;
    *) echo "error: bilinmeyen argüman: $1" >&2; exit 2;;
  esac
done

# checksum yaz (indirilen her dosyanın yanına .sha256) — sürümsüz KEGG/STRING
# anlık görüntüleri için tekrarlanabilirlik kaydı.
_stamp() { sha256sum "$1" > "$1.sha256"; echo "  sha256: $(cut -d' ' -f1 "$1.sha256")"; }
_fetch() { echo "  -> $2"; curl -fL --retry 3 -o "$2" "$1"; _stamp "$2"; }

# GO ontolojisi (m09/m12) — organizmadan bağımsız
if [ "$SKIP_OBO" -eq 0 ]; then
  echo "==> GO ontology (go-basic.obo)"
  mkdir -p references/go
  _fetch "http://purl.obolibrary.org/obo/go/go-basic.obo" "references/go/go-basic.obo"
else
  echo "==> go-basic.obo atlandı (--skip-obo)"
fi

# Organizma GO anotasyonu (GAF) — EBI-GOA (QuickGO fallback'i)
if [ -n "$GOA_URL" ]; then
  echo "==> GO annotation (GAF) — EBI-GOA"
  mkdir -p "$(dirname "$GOA_OUT")"
  _fetch "$GOA_URL" "$GOA_OUT"
else
  echo "==> GAF atlandı (--goa-url verilmedi; m09/m12 GAF-doldurma çalışmaz)"
fi

# KEGG (m10) — organizmaya özel REST dosyaları
if [ -n "$KEGG_ORG" ]; then
  echo "==> KEGG ($KEGG_ORG) REST dosyaları"
  d="references/kegg/$KEGG_ORG"; mkdir -p "$d"
  _fetch "https://rest.kegg.jp/link/pathway/$KEGG_ORG" "$d/pathway_links.tsv"
  _fetch "https://rest.kegg.jp/list/pathway/$KEGG_ORG" "$d/pathway_names.tsv"
  _fetch "https://rest.kegg.jp/list/$KEGG_ORG"         "$d/gene_list.tsv"
else
  echo "==> KEGG atlandı (--kegg-org verilmedi; m10 çalışmaz)"
fi

# STRING (m15) — taxon ağı
if [ -n "$STRING_TAXID" ]; then
  echo "==> STRING ($STRING_TAXID) v12.0 ağı"
  d="references/string/$STRING_TAXID"; mkdir -p "$d"
  _fetch "https://stringdb-downloads.org/download/protein.info.v12.0/${STRING_TAXID}.protein.info.v12.0.txt.gz" \
         "$d/protein.info.txt.gz"
  _fetch "https://stringdb-downloads.org/download/protein.links.v12.0/${STRING_TAXID}.protein.links.v12.0.txt.gz" \
         "$d/protein.links.txt.gz"
else
  echo "==> STRING atlandı (--string-taxid verilmedi; m15 çalışmaz)"
fi

echo
echo "Referans hazırlığı bitti. AMR/virulence (m13) abricate'in gömülü CARD/VFDB'sini kullanır (indirme yok)."
