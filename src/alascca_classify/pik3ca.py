"""PIK3CA variant classification: Group A hotspots and Group B non-hotspot variants."""

from __future__ import annotations

import json
import re
from pathlib import Path

from .models import Variant

_DATA_DIR = Path(__file__).parent / "data"

# Load hotspot data
with open(_DATA_DIR / "pik3ca_hotspots.json") as f:
    _HOTSPOT_DATA = json.load(f)

GROUP_A_POSITIONS: set[int] = set(_HOTSPOT_DATA["group_a_positions"])
POSITION_ANNOTATIONS: dict[str, dict] = _HOTSPOT_DATA["position_annotations"]

# Load consequence types
with open(_DATA_DIR / "consequence_types.json") as f:
    _CONSEQUENCE_DATA = json.load(f)

MAF_TO_VEP: dict[str, str] = _CONSEQUENCE_DATA["maf_to_vep_mapping"]
GROUP_B_CONSEQUENCES: set[str] = set(
    _CONSEQUENCE_DATA["high_impact"] + _CONSEQUENCE_DATA["moderate_impact"]
)

# Regex to extract protein position from HGVSp
# Matches: p.E545K, p.H1047R, p.R130*, p.R130fs, p.E545del, p.K546_E547del, etc.
_HGVSP_POSITION_RE = re.compile(r"p\.([A-Z]?)(\d+)")


def extract_protein_position(hgvsp: str) -> int | None:
    """Extract the first protein position number from an HGVSp string."""
    if not hgvsp:
        return None
    m = _HGVSP_POSITION_RE.search(hgvsp)
    if m:
        return int(m.group(2))
    return None


def map_consequence(maf_classification: str) -> str:
    """Map MAF Variant_Classification to VEP consequence term."""
    return MAF_TO_VEP.get(maf_classification, "")


def is_group_a(variant: Variant) -> bool:
    """Check if a PIK3CA variant is a Group A hotspot mutation."""
    if variant.gene.upper() != "PIK3CA":
        return False

    vep = map_consequence(variant.variant_classification)
    if vep != "missense_variant":
        return False

    pos = extract_protein_position(variant.hgvsp)
    if pos is None:
        return False

    return pos in GROUP_A_POSITIONS


def annotate_group_a(variant: Variant) -> Variant:
    """Annotate a Group A variant with exon and domain information."""
    pos = extract_protein_position(variant.hgvsp)
    if pos is not None:
        pos_str = str(pos)
        if pos_str in POSITION_ANNOTATIONS:
            ann = POSITION_ANNOTATIONS[pos_str]
            variant = variant.model_copy(
                update={"exon": ann["exon"], "domain": ann["domain"]}
            )
    variant = variant.model_copy(update={"vaf": variant.compute_vaf()})
    return variant


def is_group_b_pik3ca(variant: Variant) -> bool:
    """Check if a PIK3CA variant qualifies as Group B (non-hotspot, moderate/high impact)."""
    if variant.gene.upper() != "PIK3CA":
        return False

    # If it's Group A, it's not Group B
    if is_group_a(variant):
        return False

    vep = map_consequence(variant.variant_classification)
    return vep in GROUP_B_CONSEQUENCES


def annotate_variant(variant: Variant) -> Variant:
    """Add VAF computation to a variant."""
    return variant.model_copy(update={"vaf": variant.compute_vaf()})
