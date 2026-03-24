"""Tests for PTEN variant classification."""

from alascca_classify.models import PTENDeletion, Variant
from alascca_classify.pten import is_group_b_pten, pten_deletion_as_variant


class TestGroupBPTEN:
    def test_nonsense(self, pten_nonsense):
        assert is_group_b_pten(pten_nonsense) is True

    def test_frameshift(self, pten_frameshift):
        assert is_group_b_pten(pten_frameshift) is True

    def test_missense(self):
        v = Variant(gene="PTEN", hgvsp="p.C124S", variant_classification="Missense_Mutation")
        assert is_group_b_pten(v) is True

    def test_wrong_gene(self):
        v = Variant(gene="TP53", hgvsp="p.R130*", variant_classification="Nonsense_Mutation")
        assert is_group_b_pten(v) is False

    def test_silent_not_group_b(self):
        v = Variant(gene="PTEN", hgvsp="p.=", variant_classification="Silent")
        assert is_group_b_pten(v) is False


class TestPTENDeletionVariant:
    def test_creates_variant(self):
        deletion = PTENDeletion(detected=True)
        v = pten_deletion_as_variant(deletion)
        assert v.gene == "PTEN"
        assert v.hgvsp == "homozygous_deletion"
        assert v.variant_classification == "homozygous_loss"
