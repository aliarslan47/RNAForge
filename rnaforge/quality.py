"""Kalite eşikleri: koda gömülü DEĞİL, profiles/*.yml içinde veri olarak durur.

Yeni bir profil eklemek kod değişikliği değil, dosya eklemektir (spec 2026-07-20).
"""
from __future__ import annotations

import types
from dataclasses import dataclass, field
from pathlib import Path

import yaml

PROFILE_DIR = Path(__file__).parent / "profiles"


class ProfileError(ValueError):
    """Kalite profili yüklenemedi veya geçersiz eşik ezmesi verildi."""


@dataclass(frozen=True)
class Profile:
    # frozen=True yalnizca attribute atamasini engeller; icerideki dict/set
    # yine de mutasyona acik kalirdi. MappingProxyType + frozenset ile
    # gercek degismezlik saglaniyor - bu, is_overridden()'in dogrulugunun
    # (dolayisiyla musteri raporundaki guven kartinin) dayanagi.
    name: str
    permissive: bool
    description: str
    _thresholds: types.MappingProxyType[str, float]
    _overridden: frozenset[str] = field(default_factory=frozenset)

    def threshold(self, gate: str) -> float:
        if gate not in self._thresholds:
            raise ProfileError(
                f"profile {self.name!r} has no threshold for gate {gate!r}; "
                f"known gates: {', '.join(sorted(self._thresholds))}"
            )
        return self._thresholds[gate]

    def is_overridden(self, gate: str) -> bool:
        return gate in self._overridden

    def overrides(self) -> dict[str, float]:
        return {g: self._thresholds[g] for g in sorted(self._overridden)}


def profile_name_for(organism_type: str, read_type: str) -> str:
    """read_type'a göre profil adı. Uzun okuma ONT-permissive profili kullanır
    (`<organism_type>_long`); kısa okuma organizma profilini kullanır. read_type→profil
    eşlemesinin TEK kaynağı (modüller + CLI aynı çözümü kullansın)."""
    return f"{organism_type}_long" if read_type == "long" else organism_type


def load_profile(organism_type: str, overrides: dict | None = None) -> Profile:
    path = PROFILE_DIR / f"{organism_type}.yml"
    if not path.exists():
        available = ", ".join(sorted(p.stem for p in PROFILE_DIR.glob("*.yml")))
        raise ProfileError(
            f"no quality profile for organism_type={organism_type!r} "
            f"(available: {available})"
        )
    raw = yaml.safe_load(path.read_text()) or {}
    # YAML'daki int degerler (ör. read_depth: 1000000) burada float'a
    # sabitleniyor ki threshold()'un -> float imzasi her yolda dogru olsun.
    thresholds = {gate: float(value) for gate, value in (raw.get("thresholds") or {}).items()}

    applied: set[str] = set()
    for gate, value in (overrides or {}).items():
        if gate not in thresholds:
            # Yazım hatasını yutmak, kullanıcıya "eşiği gevşettim" yanılgısı verir.
            raise ProfileError(
                f"quality.{gate}: unknown gate for profile {organism_type!r}; "
                f"known gates: {', '.join(sorted(thresholds))}"
            )
        try:
            thresholds[gate] = float(value)
        except (TypeError, ValueError):
            raise ProfileError(f"quality.{gate}: expected a number, got {value!r}") from None
        applied.add(gate)

    return Profile(
        name=raw.get("name", organism_type),
        permissive=bool(raw.get("permissive", False)),
        description=str(raw.get("description", "")).strip(),
        _thresholds=types.MappingProxyType(thresholds),
        _overridden=frozenset(applied),
    )
