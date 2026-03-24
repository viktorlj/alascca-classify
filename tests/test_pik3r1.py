"""Tests for PIK3R1 Group B classification."""

from alascca_classify.models import Variant
from alascca_classify.pik3r1 import is_group_b_pik3r1


class TestGroupBPIK3R1:
    def test_missense(self, pik3r1_missense):
        assert is_group_b_pik3r1(pik3r1_missense) is True

    def test_frameshift(self, pik3r1_frameshift):
        assert is_group_b_pik3r1(pik3r1_frameshift) is True

    def test_nonsense(self):
        v = Variant(gene="PIK3R1", hgvsp="p.R386*", variant_classification="Nonsense_Mutation")
        assert is_group_b_pik3r1(v) is True

    def test_inframe_del(self):
        v = Variant(gene="PIK3R1", hgvsp="p.K567del", variant_classification="In_Frame_Del")
        assert is_group_b_pik3r1(v) is True

    def test_wrong_gene(self):
        v = Variant(gene="PIK3CA", hgvsp="p.R348W", variant_classification="Missense_Mutation")
        assert is_group_b_pik3r1(v) is False

    def test_silent_not_group_b(self):
        v = Variant(gene="PIK3R1", hgvsp="p.=", variant_classification="Silent")
        assert is_group_b_pik3r1(v) is False

    def test_start_lost(self):
        v = Variant(gene="PIK3R1", hgvsp="p.M1?", variant_classification="Translation_Start_Site")
        assert is_group_b_pik3r1(v) is True
