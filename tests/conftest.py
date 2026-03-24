"""Shared fixtures for ALASCCA-classify tests."""

from __future__ import annotations

import pytest

from alascca_classify.models import (
    CopyNumberSegment,
    Variant,
)


@pytest.fixture
def pik3ca_e545k() -> Variant:
    return Variant(
        gene="PIK3CA",
        hgvsp="p.E545K",
        variant_classification="Missense_Mutation",
        chromosome="3",
        t_alt_count=50,
        t_ref_count=100,
    )


@pytest.fixture
def pik3ca_h1047r() -> Variant:
    return Variant(
        gene="PIK3CA",
        hgvsp="p.H1047R",
        variant_classification="Missense_Mutation",
        chromosome="3",
    )


@pytest.fixture
def pik3ca_nonhotspot() -> Variant:
    return Variant(
        gene="PIK3CA",
        hgvsp="p.R88Q",
        variant_classification="Missense_Mutation",
        chromosome="3",
    )


@pytest.fixture
def pik3ca_frameshift() -> Variant:
    return Variant(
        gene="PIK3CA",
        hgvsp="p.R108fs",
        variant_classification="Frame_Shift_Del",
        chromosome="3",
    )


@pytest.fixture
def pik3r1_missense() -> Variant:
    return Variant(
        gene="PIK3R1",
        hgvsp="p.R348W",
        variant_classification="Missense_Mutation",
        chromosome="5",
    )


@pytest.fixture
def pik3r1_frameshift() -> Variant:
    return Variant(
        gene="PIK3R1",
        hgvsp="p.N453fs",
        variant_classification="Frame_Shift_Ins",
        chromosome="5",
    )


@pytest.fixture
def pten_nonsense() -> Variant:
    return Variant(
        gene="PTEN",
        hgvsp="p.R130*",
        variant_classification="Nonsense_Mutation",
        chromosome="10",
    )


@pytest.fixture
def pten_frameshift() -> Variant:
    return Variant(
        gene="PTEN",
        hgvsp="p.K267fs",
        variant_classification="Frame_Shift_Del",
        chromosome="10",
    )


@pytest.fixture
def braf_v600e() -> Variant:
    return Variant(
        gene="BRAF",
        hgvsp="p.V600E",
        variant_classification="Missense_Mutation",
        chromosome="7",
    )


@pytest.fixture
def kras_g12d() -> Variant:
    return Variant(
        gene="KRAS",
        hgvsp="p.G12D",
        variant_classification="Missense_Mutation",
        chromosome="12",
    )


@pytest.fixture
def nras_q61k() -> Variant:
    return Variant(
        gene="NRAS",
        hgvsp="p.Q61K",
        variant_classification="Missense_Mutation",
        chromosome="1",
    )


@pytest.fixture
def pten_deletion_segments() -> list[CopyNumberSegment]:
    """Segments showing a clear PTEN deep deletion."""
    return [
        # Upstream flanking (normal)
        CopyNumberSegment(chrom="10", start=86623195, end=89000000, seg_mean=0.0),
        CopyNumberSegment(chrom="10", start=89000000, end=89600000, seg_mean=0.0),
        # PTEN region (deleted)
        CopyNumberSegment(chrom="10", start=89623195, end=89725229, seg_mean=-0.8),
        # Downstream flanking (normal)
        CopyNumberSegment(chrom="10", start=89726000, end=91000000, seg_mean=0.0),
        CopyNumberSegment(chrom="10", start=91000000, end=92725229, seg_mean=0.0),
    ]


@pytest.fixture
def pten_no_deletion_segments() -> list[CopyNumberSegment]:
    """Segments with no PTEN deletion (all near 0)."""
    return [
        CopyNumberSegment(chrom="10", start=86623195, end=89600000, seg_mean=0.0),
        CopyNumberSegment(chrom="10", start=89623195, end=89725229, seg_mean=-0.05),
        CopyNumberSegment(chrom="10", start=89726000, end=92725229, seg_mean=0.0),
    ]


@pytest.fixture
def sample_maf_content() -> str:
    """A simple MAF file content string with mixed variants."""
    return (
        "Hugo_Symbol\tVariant_Classification\tHGVSp_Short\tChromosome\tStart_Position\tt_alt_count\tt_ref_count\n"
        "PIK3CA\tMissense_Mutation\tp.E545K\t3\t178936091\t50\t100\n"
        "KRAS\tMissense_Mutation\tp.G12D\t12\t25398284\t30\t120\n"
        "TP53\tMissense_Mutation\tp.R175H\t17\t7578406\t60\t90\n"
    )
