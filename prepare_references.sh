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
# Metatranskriptom kolu (organism_type=metatranscriptome):
#   --kraken2-db-url <url> Kraken2/Bracken DB tarball'ı (.tar.gz) → references/kraken2/<ad>/ (açılır)
#   --kraken2-db-name <ad> Kraken2 DB hedef alt-dizin adı (varsayılan: db)
#   --rrna-db-url <url>    SortMeRNA rRNA referans FASTA'sı → references/rrna/
#   -h, --help             Bu kullanım metnini yazdır ve çık
#
# NOT (engellenen kaynak): QuickGO indirmesi bazı ağlarda BLOKLU. Bu script GAF için
# doğrudan QuickGO KULLANMAZ; EBI-GOA FTP proteome dosyasını (--goa-url) kullanır — bu
# engel için belgeli fallback yoludur (DURUM 2026-08).
set -euo pipefail
cd "$(dirname "$0")"

_usage() { sed -n '2,23p' "$0" | sed 's/^# \{0,1\}//'; }

KEGG_ORG=""; STRING_TAXID=""; GOA_URL=""; GOA_OUT="references/go/organism.gaf"; SKIP_OBO=0
KRAKEN2_DB_URL=""; KRAKEN2_DB_NAME="db"; RRNA_DB_URL=""
while [ $# -gt 0 ]; do
  case "$1" in
    --kegg-org) KEGG_ORG="$2"; shift 2;;
    --string-taxid) STRING_TAXID="$2"; shift 2;;
    --goa-url) GOA_URL="$2"; shift 2;;
    --goa-out) GOA_OUT="$2"; shift 2;;
    --skip-obo) SKIP_OBO=1; shift;;
    --kraken2-db-url) KRAKEN2_DB_URL="$2"; shift 2;;
    --kraken2-db-name) KRAKEN2_DB_NAME="$2"; shift 2;;
    --rrna-db-url) RRNA_DB_URL="$2"; shift 2;;
    -h|--help) _usage; exit 0;;
    *) echo "error: bilinmeyen argüman: $1" >&2; exit 2;;
  esac
done

# checksum yaz (indirilen her dosyanın yanına .sha256) — sürümsüz KEGG/STRING
# anlık görüntüleri için tekrarlanabilirlik kaydı.
_stamp() { sha256sum "$1" > "$1.sha256"; echo "  sha256: $(cut -d' ' -f1 "$1.sha256")"; }
_fetch() { echo "  -> $2"; curl -fL --retry 3 -o "$2" "$1"; _stamp "$2"; }
# tarball indir (kaydını .sha256 ile tut) → hedef dizine aç. Kraken2/Bracken DB'leri
# tek dosya değil tarball geldiği için _fetch yetmez (extract adımı şart).
_fetch_tar() {  # $1=url  $2=hedef_dizin
  local url="$1" dir="$2" tarball
  mkdir -p "$dir"; tarball="$dir/_download.tar.gz"
  echo "  -> $dir (tarball açılıyor)"
  curl -fL --retry 3 -o "$tarball" "$url"; _stamp "$tarball"
  tar -xzf "$tarball" -C "$dir"; rm -f "$tarball"
}

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

# --- Metatranskriptom kolu (organism_type=metatranscriptome) ---------------------

# Kraken2/Bracken DB (m_taxonomy) — tarball indir + aç → config taxonomy.kraken2_db bu dizini gösterir
if [ -n "$KRAKEN2_DB_URL" ]; then
  echo "==> Kraken2/Bracken DB → references/kraken2/$KRAKEN2_DB_NAME"
  _fetch_tar "$KRAKEN2_DB_URL" "references/kraken2/$KRAKEN2_DB_NAME"
else
  echo "==> Kraken2 DB atlandı (--kraken2-db-url verilmedi; m_taxonomy çalışmaz)"
fi

# SortMeRNA rRNA referansı (m_rrna_deplete) — config rrna.db_fasta bu FASTA'yı gösterir
if [ -n "$RRNA_DB_URL" ]; then
  echo "==> rRNA referansı (SortMeRNA) → references/rrna/"
  mkdir -p references/rrna
  _fetch "$RRNA_DB_URL" "references/rrna/$(basename "$RRNA_DB_URL")"
else
  echo "==> rRNA DB atlandı (--rrna-db-url verilmedi; m_rrna_deplete depletion çalışmaz)"
fi

echo
echo "Referans hazırlığı bitti. AMR/virulence (m13) abricate'in gömülü CARD/VFDB'sini kullanır (indirme yok)."
