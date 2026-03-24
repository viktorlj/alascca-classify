"""End-to-end classification tests."""


from alascca_classify.classify import classify_sample
from alascca_classify.models import (
    SampleInput,
    SampleMetadata,
    Variant,
)


class TestClassifySample:
    def test_group_a(self, pik3ca_e545k):
        sample = SampleInput(
            metadata=SampleMetadata(sample_id="TEST_A"),
            variants=[pik3ca_e545k],
        )
        result = classify_sample(sample)
        assert result.alascca_group == "A"
        assert result.has_group_a is True
        assert result.aspirin_eligible is True
        assert len(result.group_a_variants) == 1

    def test_group_b_pten(self, pten_nonsense):
        sample = SampleInput(
            metadata=SampleMetadata(sample_id="TEST_B"),
            variants=[pten_nonsense],
        )
        result = classify_sample(sample)
        assert result.alascca_group == "B"
        assert result.has_group_b is True
        assert result.aspirin_eligible is True

    def test_group_b_pik3r1(self, pik3r1_missense):
        sample = SampleInput(
            metadata=SampleMetadata(sample_id="TEST_B2"),
            variants=[pik3r1_missense],
        )
        result = classify_sample(sample)
        assert result.alascca_group == "B"

    def test_group_b_pik3ca_nonhotspot(self, pik3ca_nonhotspot):
        sample = SampleInput(
            metadata=SampleMetadata(sample_id="TEST_B3"),
            variants=[pik3ca_nonhotspot],
        )
        result = classify_sample(sample)
        assert result.alascca_group == "B"

    def test_priority_a_over_b(self, pik3ca_e545k, pten_nonsense):
        """A + B = A."""
        sample = SampleInput(
            metadata=SampleMetadata(sample_id="TEST_AB"),
            variants=[pik3ca_e545k, pten_nonsense],
        )
        result = classify_sample(sample)
        assert result.alascca_group == "A"
        assert result.has_group_a is True
        assert result.has_group_b is True
        assert len(result.flags) > 0

    def test_no_pi3k_alteration(self, kras_g12d, braf_v600e):
        sample = SampleInput(
            metadata=SampleMetadata(sample_id="TEST_NONE"),
            variants=[kras_g12d, braf_v600e],
        )
        result = classify_sample(sample)
        assert result.alascca_group == "none"
        assert result.aspirin_eligible is False
        # But biomarkers should still be detected
        assert result.informational_biomarkers.kras_mutant is True
        assert result.informational_biomarkers.braf_v600e is True

    def test_empty_variants(self):
        sample = SampleInput(metadata=SampleMetadata(sample_id="TEST_EMPTY"))
        result = classify_sample(sample)
        assert result.alascca_group == "none"
        assert result.aspirin_eligible is False

    def test_pten_deletion_from_segments(self, pten_deletion_segments):
        sample = SampleInput(
            metadata=SampleMetadata(sample_id="TEST_DEL", tumor_purity=0.50),
            segments=pten_deletion_segments,
        )
        result = classify_sample(sample)
        assert result.alascca_group == "B"
        assert result.pten_deletion is not None
        assert result.pten_deletion.detected is True

    def test_all_hotspot_positions(self):
        """Verify each of the 7 hotspot positions classifies as Group A."""
        hotspots = [
            ("p.E542K", 542),
            ("p.E545K", 545),
            ("p.Q546R", 546),
            ("p.H1021R", 1021),
            ("p.N1043K", 1043),
            ("p.H1044R", 1044),
            ("p.H1047R", 1047),
        ]
        for hgvsp, pos in hotspots:
            v = Variant(gene="PIK3CA", hgvsp=hgvsp, variant_classification="Missense_Mutation")
            sample = SampleInput(variants=[v])
            result = classify_sample(sample)
            assert result.alascca_group == "A", f"Failed for {hgvsp}"

    def test_clinical_notes_group_a(self, pik3ca_e545k):
        sample = SampleInput(variants=[pik3ca_e545k])
        result = classify_sample(sample)
        assert any("Group A" in n for n in result.clinical_notes)
        assert any("aspirin" in n.lower() for n in result.clinical_notes)

    def test_clinical_notes_no_alteration(self):
        sample = SampleInput()
        result = classify_sample(sample)
        assert any("No PI3K" in n for n in result.clinical_notes)

    def test_biomarkers_included(self, pik3ca_e545k, kras_g12d):
        sample = SampleInput(
            metadata=SampleMetadata(msi_status="MSS"),
            variants=[pik3ca_e545k, kras_g12d],
        )
        result = classify_sample(sample)
        assert result.informational_biomarkers.kras_mutant is True
        assert result.informational_biomarkers.msi_status == "MSS"
