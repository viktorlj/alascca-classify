"""Tests for PIK3CA Group A hotspot and Group B non-hotspot classification."""

import pytest

from alascca_classify.models import Variant
from alascca_classify.pik3ca import (
    annotate_group_a,
    extract_protein_position,
    is_group_a,
    is_group_b_pik3ca,
    map_consequence,
)


class TestExtractProteinPosition:
    def test_standard_missense(self):
        assert extract_protein_position("p.E545K") == 545

    def test_four_digit_position(self):
        assert extract_protein_position("p.H1047R") == 1047

    def test_nonsense(self):
        assert extract_protein_position("p.R130*") == 130

    def test_frameshift(self):
        assert extract_protein_position("p.R108fs") == 108

    def test_empty(self):
        assert extract_protein_position("") is None

    def test_no_match(self):
        assert extract_protein_position("c.1234A>G") is None

    def test_deletion(self):
        assert extract_protein_position("p.K546del") == 546


class TestMapConsequence:
    def test_missense(self):
        assert map_consequence("Missense_Mutation") == "missense_variant"

    def test_frameshift_del(self):
        assert map_consequence("Frame_Shift_Del") == "frameshift_variant"

    def test_frameshift_ins(self):
        assert map_consequence("Frame_Shift_Ins") == "frameshift_variant"

    def test_nonsense(self):
        assert map_consequence("Nonsense_Mutation") == "stop_gained"

    def test_splice(self):
        assert map_consequence("Splice_Site") == "splice_site"

    def test_in_frame_del(self):
        assert map_consequence("In_Frame_Del") == "inframe_deletion"

    def test_in_frame_ins(self):
        assert map_consequence("In_Frame_Ins") == "inframe_insertion"

    def test_start_lost(self):
        assert map_consequence("Translation_Start_Site") == "start_lost"

    def test_unknown(self):
        assert map_consequence("Silent") == ""


class TestGroupA:
    """Test all 7 Group A hotspot positions."""

    @pytest.mark.parametrize(
        "hgvsp,position",
        [
            ("p.E542K", 542),
            ("p.E545K", 545),
            ("p.Q546K", 546),
            ("p.H1021R", 1021),
            ("p.N1043K", 1043),
            ("p.H1044R", 1044),
            ("p.H1047R", 1047),
        ],
    )
    def test_all_hotspot_positions(self, hgvsp, position):
        v = Variant(gene="PIK3CA", hgvsp=hgvsp, variant_classification="Missense_Mutation")
        assert is_group_a(v) is True

    def test_hotspot_alternate_change(self):
        """Different amino acid change at same position should still be Group A."""
        v = Variant(gene="PIK3CA", hgvsp="p.H1047L", variant_classification="Missense_Mutation")
        assert is_group_a(v) is True

    def test_non_hotspot_position(self):
        v = Variant(gene="PIK3CA", hgvsp="p.R88Q", variant_classification="Missense_Mutation")
        assert is_group_a(v) is False

    def test_wrong_gene(self):
        v = Variant(gene="PIK3R1", hgvsp="p.E545K", variant_classification="Missense_Mutation")
        assert is_group_a(v) is False

    def test_frameshift_at_hotspot_is_not_group_a(self):
        """Frameshift at a hotspot position should not be Group A (only missense qualifies)."""
        v = Variant(gene="PIK3CA", hgvsp="p.E545fs", variant_classification="Frame_Shift_Del")
        assert is_group_a(v) is False

    def test_case_insensitive_gene(self):
        v = Variant(gene="pik3ca", hgvsp="p.E545K", variant_classification="Missense_Mutation")
        assert is_group_a(v) is True


class TestGroupBPIK3CA:
    def test_nonhotspot_missense(self, pik3ca_nonhotspot):
        assert is_group_b_pik3ca(pik3ca_nonhotspot) is True

    def test_frameshift(self, pik3ca_frameshift):
        assert is_group_b_pik3ca(pik3ca_frameshift) is True

    def test_hotspot_is_not_group_b(self, pik3ca_e545k):
        assert is_group_b_pik3ca(pik3ca_e545k) is False

    def test_nonsense(self):
        v = Variant(gene="PIK3CA", hgvsp="p.Q546*", variant_classification="Nonsense_Mutation")
        assert is_group_b_pik3ca(v) is True

    def test_inframe_del(self):
        v = Variant(gene="PIK3CA", hgvsp="p.E110del", variant_classification="In_Frame_Del")
        assert is_group_b_pik3ca(v) is True

    def test_inframe_ins(self):
        v = Variant(gene="PIK3CA", hgvsp="p.D350_N351insD", variant_classification="In_Frame_Ins")
        assert is_group_b_pik3ca(v) is True

    def test_splice_site(self):
        v = Variant(gene="PIK3CA", hgvsp="", variant_classification="Splice_Site")
        # splice_site maps to "splice_site" which is NOT in GROUP_B_CONSEQUENCES
        # (only splice_acceptor_variant and splice_donor_variant are)
        assert is_group_b_pik3ca(v) is False

    def test_start_lost(self):
        v = Variant(gene="PIK3CA", hgvsp="p.M1?", variant_classification="Translation_Start_Site")
        assert is_group_b_pik3ca(v) is True


class TestAnnotateGroupA:
    def test_annotate_e545k(self, pik3ca_e545k):
        annotated = annotate_group_a(pik3ca_e545k)
        assert annotated.exon == 9
        assert annotated.domain == "helical"

    def test_annotate_h1047r(self, pik3ca_h1047r):
        annotated = annotate_group_a(pik3ca_h1047r)
        assert annotated.exon == 20
        assert annotated.domain == "kinase"

    def test_annotate_computes_vaf(self, pik3ca_e545k):
        annotated = annotate_group_a(pik3ca_e545k)
        assert annotated.vaf == pytest.approx(50 / 150, rel=1e-3)
