"""Tests for informational biomarker detection."""

from alascca_classify.biomarkers import detect_biomarkers
from alascca_classify.models import Variant


class TestBiomarkers:
    def test_braf_v600e(self, braf_v600e):
        result = detect_biomarkers([braf_v600e])
        assert result.braf_v600e is True
        assert "p.V600E" in result.braf_variants

    def test_braf_non_v600e(self):
        v = Variant(gene="BRAF", hgvsp="p.K601E", variant_classification="Missense_Mutation")
        result = detect_biomarkers([v])
        assert result.braf_v600e is False
        assert "p.K601E" in result.braf_variants

    def test_kras(self, kras_g12d):
        result = detect_biomarkers([kras_g12d])
        assert result.kras_mutant is True
        assert "p.G12D" in result.kras_variants

    def test_nras(self, nras_q61k):
        result = detect_biomarkers([nras_q61k])
        assert result.nras_mutant is True
        assert "p.Q61K" in result.nras_variants

    def test_msi_status_from_metadata(self):
        result = detect_biomarkers([], msi_status="MSI-H")
        assert result.msi_status == "MSI-H"

    def test_no_biomarkers(self):
        v = Variant(gene="TP53", hgvsp="p.R175H", variant_classification="Missense_Mutation")
        result = detect_biomarkers([v])
        assert result.braf_v600e is False
        assert result.kras_mutant is False
        assert result.nras_mutant is False

    def test_multiple_kras(self):
        v1 = Variant(gene="KRAS", hgvsp="p.G12D", variant_classification="Missense_Mutation")
        v2 = Variant(gene="KRAS", hgvsp="p.G13D", variant_classification="Missense_Mutation")
        result = detect_biomarkers([v1, v2])
        assert result.kras_mutant is True
        assert len(result.kras_variants) == 2
