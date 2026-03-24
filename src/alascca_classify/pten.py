"""PTEN variant classification: moderate/high-impact variants + deep deletion → Group B."""

from __future__ import annotations

from .models import PTENDeletion, Variant
from .pik3ca import GROUP_B_CONSEQUENCES, map_consequence


def is_group_b_pten(variant: Variant) -> bool:
    """Check if a PTEN variant qualifies as Group B."""
    if variant.gene.upper() != "PTEN":
        return False

    vep = map_consequence(variant.variant_classification)
    return vep in GROUP_B_CONSEQUENCES


def pten_deletion_as_variant(deletion: PTENDeletion) -> Variant:
    """Create a Variant representing PTEN homozygous deletion for reporting."""
    return Variant(
        gene="PTEN",
        hgvsp="homozygous_deletion",
        variant_classification="homozygous_loss",
        chromosome="10",
        domain="deep_deletion",
    )
