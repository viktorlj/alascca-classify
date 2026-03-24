"""Tests for PTEN deep deletion detection from SEG files."""


from alascca_classify.cna import assess_pten_deletion
from alascca_classify.models import CopyNumberSegment


class TestPTENDeletion:
    def test_clear_deletion(self, pten_deletion_segments):
        result = assess_pten_deletion(pten_deletion_segments)
        assert result.detected is True
        assert result.passes_threshold is True
        assert result.seg_mean_pten is not None
        assert result.seg_mean_pten < -0.2

    def test_no_deletion(self, pten_no_deletion_segments):
        result = assess_pten_deletion(pten_no_deletion_segments)
        assert result.detected is False

    def test_empty_segments(self):
        result = assess_pten_deletion([])
        assert result.detected is False
        assert "No segments provided" in result.note

    def test_low_tumor_purity(self, pten_deletion_segments):
        result = assess_pten_deletion(pten_deletion_segments, tumor_purity=0.20)
        assert result.detected is False
        assert "purity" in result.note.lower()

    def test_adequate_tumor_purity(self, pten_deletion_segments):
        result = assess_pten_deletion(pten_deletion_segments, tumor_purity=0.45)
        assert result.detected is True

    def test_no_chr10_segments(self):
        segments = [
            CopyNumberSegment(chrom="1", start=1000, end=2000, seg_mean=-0.5),
        ]
        result = assess_pten_deletion(segments)
        assert result.detected is False
        assert "chromosome 10" in result.note.lower()

    def test_chr_prefix_handled(self):
        """Segments with 'chr10' should work the same as '10'."""
        segments = [
            CopyNumberSegment(chrom="chr10", start=86623195, end=89600000, seg_mean=0.0),
            CopyNumberSegment(chrom="chr10", start=89623195, end=89725229, seg_mean=-0.8),
            CopyNumberSegment(chrom="chr10", start=89726000, end=92725229, seg_mean=0.0),
        ]
        result = assess_pten_deletion(segments)
        assert result.detected is True

    def test_borderline_deletion(self):
        """Log-ratio diff exactly at threshold."""
        segments = [
            CopyNumberSegment(chrom="10", start=86623195, end=89600000, seg_mean=0.0),
            CopyNumberSegment(chrom="10", start=89623195, end=89725229, seg_mean=-0.2),
            CopyNumberSegment(chrom="10", start=89726000, end=92725229, seg_mean=0.0),
        ]
        result = assess_pten_deletion(segments)
        # -0.2 - 0.0 = -0.2, which equals threshold
        assert result.detected is True
