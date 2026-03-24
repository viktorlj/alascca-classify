"""Copy number analysis: PTEN deep deletion detection from SEG files."""

from __future__ import annotations

import json
from pathlib import Path

from .models import CopyNumberSegment, PTENDeletion

_DATA_DIR = Path(__file__).parent / "data"

with open(_DATA_DIR / "pten_coordinates.json") as f:
    _PTEN_COORDS = json.load(f)

PTEN_CHROM = _PTEN_COORDS["chromosome"]
PTEN_GENE_START = _PTEN_COORDS["gene_start"]
PTEN_GENE_END = _PTEN_COORDS["gene_end"]
PTEN_EXONS = _PTEN_COORDS["exons"]
FLANKING_UPSTREAM_START = _PTEN_COORDS["flanking_upstream_3mb"]
FLANKING_DOWNSTREAM_END = _PTEN_COORDS["flanking_downstream_3mb"]
DELETION_THRESHOLD = _PTEN_COORDS["deletion_log_ratio_threshold"]


def _normalize_chrom(chrom: str) -> str:
    """Normalize chromosome string to bare number/letter."""
    return chrom.replace("chr", "").strip()


def _segments_on_chrom(segments: list[CopyNumberSegment], chrom: str) -> list[CopyNumberSegment]:
    """Filter segments for a specific chromosome."""
    target = _normalize_chrom(chrom)
    return [s for s in segments if _normalize_chrom(s.chrom) == target]


def _overlaps_any_exon(seg: CopyNumberSegment) -> bool:
    """Check if a segment overlaps at least one PTEN exon."""
    for exon in PTEN_EXONS:
        if seg.start <= exon["end"] and seg.end >= exon["start"]:
            return True
    return False


def _weighted_mean_seg(segments: list[CopyNumberSegment], start: int, end: int) -> float | None:
    """Compute weighted mean seg_mean for segments overlapping a region.

    Weights by the length of overlap within the region.
    """
    total_weight = 0.0
    weighted_sum = 0.0

    for seg in segments:
        overlap_start = max(seg.start, start)
        overlap_end = min(seg.end, end)
        if overlap_start < overlap_end:
            weight = overlap_end - overlap_start
            weighted_sum += seg.seg_mean * weight
            total_weight += weight

    if total_weight == 0:
        return None

    return weighted_sum / total_weight


def assess_pten_deletion(
    segments: list[CopyNumberSegment],
    tumor_purity: float | None = None,
) -> PTENDeletion:
    """Assess PTEN deep deletion from SEG segments.

    Per ALASCCA supplementary methods:
    - Tumor cell fraction must be >= 0.30
    - A deletion is called if the PTEN region seg_mean is at least 0.2 (absolute log-ratio)
      lower than BOTH the upstream 3Mb and downstream 3Mb flanking regions.
    """
    result = PTENDeletion()

    if not segments:
        result.note = "No segments provided"
        return result

    # Check tumor purity if provided
    if tumor_purity is not None and tumor_purity < 0.30:
        result.note = f"Tumor purity {tumor_purity:.2f} below threshold (0.30)"
        return result

    chr10_segments = _segments_on_chrom(segments, PTEN_CHROM)
    if not chr10_segments:
        result.note = "No segments on chromosome 10"
        return result

    # Find segments overlapping PTEN exons
    pten_segments = [s for s in chr10_segments if _overlaps_any_exon(s)]
    if not pten_segments:
        result.note = "No segments overlapping PTEN exons"
        return result

    # Compute mean log-ratio for PTEN, upstream flanking, and downstream flanking
    seg_mean_pten = _weighted_mean_seg(chr10_segments, PTEN_GENE_START, PTEN_GENE_END)
    seg_mean_upstream = _weighted_mean_seg(
        chr10_segments, FLANKING_UPSTREAM_START, PTEN_GENE_START
    )
    seg_mean_downstream = _weighted_mean_seg(
        chr10_segments, PTEN_GENE_END, FLANKING_DOWNSTREAM_END
    )

    if seg_mean_pten is None:
        result.note = "Could not compute PTEN region seg_mean"
        return result

    result.seg_mean_pten = round(seg_mean_pten, 4)

    # Check against upstream flanking
    diff_upstream = None
    if seg_mean_upstream is not None:
        diff_upstream = seg_mean_pten - seg_mean_upstream
        result.seg_mean_upstream = round(seg_mean_upstream, 4)
        result.log_ratio_diff_upstream = round(diff_upstream, 4)

    # Check against downstream flanking
    diff_downstream = None
    if seg_mean_downstream is not None:
        diff_downstream = seg_mean_pten - seg_mean_downstream
        result.seg_mean_downstream = round(seg_mean_downstream, 4)
        result.log_ratio_diff_downstream = round(diff_downstream, 4)

    # Both flanking regions must be available and both differences must exceed threshold
    if diff_upstream is not None and diff_downstream is not None:
        passes = (diff_upstream <= DELETION_THRESHOLD) and (
            diff_downstream <= DELETION_THRESHOLD
        )
        result.passes_threshold = passes
        result.detected = passes
        if passes:
            result.note = (
                f"PTEN deep deletion detected: seg_mean PTEN={seg_mean_pten:.3f}, "
                f"upstream={seg_mean_upstream:.3f} (diff={diff_upstream:.3f}), "
                f"downstream={seg_mean_downstream:.3f} (diff={diff_downstream:.3f})"
            )
        else:
            result.note = (
                f"PTEN deletion not detected: log-ratio differences "
                f"(upstream={diff_upstream:.3f}, downstream={diff_downstream:.3f}) "
                f"do not both exceed threshold ({DELETION_THRESHOLD})"
            )
    else:
        missing = []
        if diff_upstream is None:
            missing.append("upstream")
        if diff_downstream is None:
            missing.append("downstream")
        result.note = f"Missing flanking region data: {', '.join(missing)}"

    return result
