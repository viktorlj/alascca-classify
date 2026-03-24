"""Tests for MAF and SEG file parsing."""


from alascca_classify.io import parse_maf, parse_seg


class TestParseMAF:
    def test_basic_maf(self, sample_maf_content):
        variants = parse_maf(sample_maf_content)
        assert len(variants) == 3
        assert variants[0].gene == "PIK3CA"
        assert variants[0].hgvsp == "p.E545K"
        assert variants[0].variant_classification == "Missense_Mutation"
        assert variants[0].t_alt_count == 50
        assert variants[0].t_ref_count == 100

    def test_comment_lines_skipped(self):
        content = (
            "#version 2.4\n"
            "#Some comment\n"
            "Hugo_Symbol\tVariant_Classification\tHGVSp_Short\n"
            "PIK3CA\tMissense_Mutation\tp.H1047R\n"
        )
        variants = parse_maf(content)
        assert len(variants) == 1
        assert variants[0].hgvsp == "p.H1047R"

    def test_hgvsp_column_fallback(self):
        """Use HGVSp if HGVSp_Short is missing."""
        content = (
            "Hugo_Symbol\tVariant_Classification\tHGVSp\n"
            "PTEN\tNonsense_Mutation\tp.R130*\n"
        )
        variants = parse_maf(content)
        assert len(variants) == 1
        assert variants[0].hgvsp == "p.R130*"

    def test_missing_optional_columns(self):
        content = (
            "Hugo_Symbol\tVariant_Classification\tHGVSp_Short\n"
            "PIK3R1\tMissense_Mutation\tp.R348W\n"
        )
        variants = parse_maf(content)
        assert len(variants) == 1
        assert variants[0].chromosome == ""
        assert variants[0].t_alt_count is None

    def test_empty_maf(self):
        content = "Hugo_Symbol\tVariant_Classification\tHGVSp_Short\n"
        variants = parse_maf(content)
        assert len(variants) == 0


class TestParseSEG:
    def test_basic_seg(self):
        content = (
            "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n"
            "SAMPLE_001\t10\t89623195\t89725229\t50\t-0.8\n"
            "SAMPLE_001\t10\t86623195\t89600000\t200\t0.0\n"
        )
        segments = parse_seg(content)
        assert len(segments) == 2
        assert segments[0].chrom == "10"
        assert segments[0].seg_mean == -0.8
        assert segments[0].sample_id == "SAMPLE_001"

    def test_chr_prefix(self):
        content = (
            "ID\tchrom\tloc.start\tloc.end\tnum.mark\tseg.mean\n"
            "S1\tchr10\t89623195\t89725229\t50\t-0.5\n"
        )
        segments = parse_seg(content)
        assert segments[0].chrom == "10"

    def test_alternative_column_names(self):
        content = (
            "Sample\tChromosome\tStart\tEnd\tnum_mark\tsegment_mean\n"
            "S1\t10\t100000\t200000\t30\t0.1\n"
        )
        segments = parse_seg(content)
        assert len(segments) == 1
        assert segments[0].chrom == "10"
