"""PIK3R1 variant classification: all moderate/high-impact variants → Group B."""

from __future__ import annotations

from .models import Variant
from .pik3ca import GROUP_B_CONSEQUENCES, map_consequence


def is_group_b_pik3r1(variant: Variant) -> bool:
    """Check if a PIK3R1 variant qualifies as Group B."""
    if variant.gene.upper() != "PIK3R1":
        return False

    vep = map_consequence(variant.variant_classification)
    return vep in GROUP_B_CONSEQUENCES
