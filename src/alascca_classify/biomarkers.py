"""Informational biomarker detection: BRAF V600E, KRAS, NRAS, MSI status."""

from __future__ import annotations

import re

from .models import InformationalBiomarkers, Variant

_HGVSP_RE = re.compile(r"p\.([A-Z])(\d+)([A-Z*])")


def _is_braf_v600e(variant: Variant) -> bool:
    """Check if a BRAF variant is V600E."""
    if variant.gene.upper() != "BRAF":
        return False
    hgvsp = variant.hgvsp
    if not hgvsp:
        return False
    # Match p.V600E specifically
    m = _HGVSP_RE.search(hgvsp)
    if m and m.group(1) == "V" and m.group(2) == "600" and m.group(3) == "E":
        return True
    return False


def _is_braf_variant(variant: Variant) -> bool:
    """Check if variant is in BRAF."""
    return variant.gene.upper() == "BRAF"


def _is_kras_variant(variant: Variant) -> bool:
    """Check if variant is in KRAS."""
    return variant.gene.upper() == "KRAS"


def _is_nras_variant(variant: Variant) -> bool:
    """Check if variant is in NRAS."""
    return variant.gene.upper() == "NRAS"


def detect_biomarkers(
    variants: list[Variant], msi_status: str = ""
) -> InformationalBiomarkers:
    """Detect informational biomarkers from variant list and metadata."""
    result = InformationalBiomarkers(msi_status=msi_status)

    for v in variants:
        if _is_braf_v600e(v):
            result.braf_v600e = True
        if _is_braf_variant(v) and v.hgvsp:
            result.braf_variants.append(v.hgvsp)

        if _is_kras_variant(v) and v.hgvsp:
            result.kras_mutant = True
            result.kras_variants.append(v.hgvsp)

        if _is_nras_variant(v) and v.hgvsp:
            result.nras_mutant = True
            result.nras_variants.append(v.hgvsp)

    return result
