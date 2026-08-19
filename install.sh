#!/usr/bin/env bash
# RNAForge tek-komut kurulum: tüm conda env'lerini oluşturur + paketi editable kurar.
# Kullanım:  bash install.sh
# Sonra doğrula:  conda run -n rnaforge-core rnaforge doctor
#
# NOT: dorado (ONT ham-sinyal basecalling, m00) conda DIŞINDA, GPU-only bir binary'dir
# ve buraya DAHİL DEĞİLDİR; yalnız FAST5/POD5 girdisi kullanacaksan ayrıca kur ve
# config'te basecall.dorado_bin ile yolunu ver. Referans verileri: prepare_references.sh.
set -euo pipefail

cd "$(dirname "$0")"

if ! command -v conda >/dev/null 2>&1; then
  echo "error: 'conda' PATH'te değil. Miniconda/Miniforge kur ve tekrar dene." >&2
  exit 1
fi

# conda activate'i script içinde kullanabilmek için:
source "$(conda info --base)/etc/profile.d/conda.sh"

echo "==> conda env'leri oluşturuluyor (envs/*.yml)"
for yml in envs/*.yml; do
  name="$(basename "$yml" .yml)"
  if conda env list | awk '{print $1}' | grep -qx "$name"; then
    echo "  [var]    $name — atlanıyor (silmek için: conda env remove -n $name)"
  else
    echo "  [oluştur] $name"
    conda env create -f "$yml"
  fi
done

echo "==> rnaforge paketi editable kuruluyor (rnaforge-core içine)"
conda run -n rnaforge-core pip install -e .

echo
echo "Kurulum bitti. Doğrula:"
echo "  conda run -n rnaforge-core rnaforge doctor"
