"""Main classification orchestrator: combines all modules to produce ALASCCA classification."""

from __future__ import annotations

from . import __version__
from .biomarkers import detect_biomarkers
from .cna import assess_pten_deletion
from .models import ClassificationResult, SampleInput, Variant
from .pik3ca import annotate_group_a, annotate_variant, is_group_a, is_group_b_pik3ca
from .pik3r1 import is_group_b_pik3r1
from .pten import is_group_b_pten, pten_deletion_as_variant


def classify_sample(sample: SampleInput) -> ClassificationResult:
    """Classify a sample's PI3K pathway alterations per ALASCCA criteria.

    Priority: Group A > Group B > none.
    A sample with both Group A and Group B alterations is classified as Group A.
    """
    group_a_variants: list[Variant] = []
    group_b_variants: list[Variant] = []
    clinical_notes: list[str] = []
    flags: list[str] = []

    # Assess each variant
    for variant in sample.variants:
        # PIK3CA Group A
        if is_group_a(variant):
            group_a_variants.append(annotate_group_a(variant))
            continue

        # PIK3CA Group B (non-hotspot)
        if is_group_b_pik3ca(variant):
            group_b_variants.append(annotate_variant(variant))
            continue

        # PIK3R1 Group B
        if is_group_b_pik3r1(variant):
            group_b_variants.append(annotate_variant(variant))
            continue

        # PTEN Group B (somatic variant)
        if is_group_b_pten(variant):
            group_b_variants.append(annotate_variant(variant))
            continue

    # PTEN deep deletion from CNA segments
    pten_deletion = None
    if sample.segments:
        pten_deletion = assess_pten_deletion(
            sample.segments, tumor_purity=sample.metadata.tumor_purity
        )
        if pten_deletion.detected:
            group_b_variants.append(pten_deletion_as_variant(pten_deletion))

    # Determine group (A takes priority)
    has_a = len(group_a_variants) > 0
    has_b = len(group_b_variants) > 0

    if has_a:
        group = "A"
    elif has_b:
        group = "B"
    else:
        group = "none"

    aspirin_eligible = group in ("A", "B")

    # Generate clinical notes
    clinical_notes = _generate_clinical_notes(
        group, group_a_variants, group_b_variants, has_a, has_b, pten_deletion
    )

    # Generate flags
    if has_a and has_b:
        flags.append(
            "Patient has both Group A and Group B alterations; classified as Group A per protocol."
        )

    # Informational biomarkers
    biomarkers = detect_biomarkers(
        sample.variants, msi_status=sample.metadata.msi_status
    )

    return ClassificationResult(
        sample_id=sample.metadata.sample_id,
        alascca_group=group,
        group_a_variants=group_a_variants,
        group_b_variants=group_b_variants,
        has_group_a=has_a,
        has_group_b=has_b,
        aspirin_eligible=aspirin_eligible,
        informational_biomarkers=biomarkers,
        pten_deletion=pten_deletion,
        clinical_notes=clinical_notes,
        flags=flags,
        version=__version__,
    )


def _generate_clinical_notes(
    group: str,
    group_a_variants: list[Variant],
    group_b_variants: list[Variant],
    has_a: bool,
    has_b: bool,
    pten_deletion=None,
) -> list[str]:
    """Generate clinical interpretation notes."""
    notes: list[str] = []

    if group == "A":
        for v in group_a_variants:
            domain_str = f", {v.domain} domain" if v.domain else ""
            exon_str = f"exon {v.exon}" if v.exon else "unknown exon"
            notes.append(
                f"PIK3CA {v.hgvsp} is a Group A hotspot mutation ({exon_str}{domain_str})."
            )
        notes.append(
            "This tumor is eligible for adjuvant aspirin per ALASCCA (HR 0.49 for Group A)."
        )
    elif group == "B":
        for v in group_b_variants:
            if v.variant_classification == "homozygous_loss":
                notes.append("PTEN homozygous deletion detected (Group B).")
            else:
                notes.append(
                    f"{v.gene} {v.hgvsp} ({v.variant_classification}) is a Group B alteration."
                )
        notes.append(
            "This tumor is eligible for adjuvant aspirin per ALASCCA (HR 0.42 for Group B)."
        )
    else:
        notes.append(
            "No PI3K pathway alterations detected."
            " Not eligible for ALASCCA-based aspirin."
        )

    if group in ("A", "B"):
        notes.append(
            "ALASCCA trial: 160 mg aspirin daily for 3 years"
            " post-resection of stage I\u2013III CRC."
        )

    return notes
