"""MAF and SEG file parsing with stdlib csv."""

from __future__ import annotations

import csv
from io import StringIO
from pathlib import Path

from .models import CopyNumberSegment, Variant

# Column name normalization for MAF files
_MAF_COLUMN_MAP = {
    "hugo_symbol": "Hugo_Symbol",
    "gene": "Hugo_Symbol",
    "variant_classification": "Variant_Classification",
    "variant_type": "Variant_Type",
    "hgvsp_short": "HGVSp_Short",
    "hgvsp": "HGVSp",
    "chromosome": "Chromosome",
    "chr": "Chromosome",
    "chrom": "Chromosome",
    "start_position": "Start_Position",
    "start_pos": "Start_Position",
    "reference_allele": "Reference_Allele",
    "ref": "Reference_Allele",
    "tumor_seq_allele2": "Tumor_Seq_Allele2",
    "alt": "Tumor_Seq_Allele2",
    "t_alt_count": "t_alt_count",
    "t_ref_count": "t_ref_count",
}

# SEG column name normalization
_SEG_COLUMN_MAP = {
    "id": "ID",
    "sample": "ID",
    "sample_id": "ID",
    "chrom": "chrom",
    "chromosome": "chrom",
    "chr": "chrom",
    "loc.start": "loc.start",
    "start": "loc.start",
    "loc.end": "loc.end",
    "end": "loc.end",
    "num.mark": "num.mark",
    "num_mark": "num.mark",
    "num.probes": "num.mark",
    "seg.mean": "seg.mean",
    "seg_mean": "seg.mean",
    "segment_mean": "seg.mean",
}


def _normalize_columns(
    fieldnames: list[str], column_map: dict[str, str]
) -> dict[str, str]:
    """Build a rename mapping from raw fieldnames using column_map (case-insensitive)."""
    rename = {}
    used_targets: set[str] = set()
    for col in fieldnames:
        lower = col.lower().replace(" ", "_")
        if lower in column_map:
            target = column_map[lower]
            if target not in used_targets:
                rename[col] = target
                used_targets.add(target)
        elif col in column_map.values():
            rename[col] = col  # already normalized
            used_targets.add(col)
    return rename


def _read_tsv(source: str | Path) -> list[dict[str, str]]:
    """Read a TSV from a file path or string content, returning list of row dicts."""
    is_file = isinstance(source, Path) or (
        isinstance(source, str) and "\t" not in source and "\n" not in source
    )
    if is_file:
        text = Path(source).read_text()
    else:
        text = source

    # Strip comment lines (MAF headers)
    lines = [line for line in text.splitlines() if not line.startswith("#")]
    reader = csv.DictReader(StringIO("\n".join(lines)), delimiter="\t")
    return list(reader)


def _remap_row(row: dict[str, str], rename: dict[str, str]) -> dict[str, str]:
    """Apply column rename mapping to a single row dict."""
    return {rename.get(k, k): v for k, v in row.items()}


def parse_maf(source: str | Path) -> list[Variant]:
    """Parse a MAF file (path or string content) into a list of Variants."""
    rows = _read_tsv(source)
    if not rows:
        return []

    rename = _normalize_columns(list(rows[0].keys()), _MAF_COLUMN_MAP)

    # Determine protein change column after normalization
    sample_renamed = {rename.get(k, k) for k in rows[0].keys()}
    hgvsp_col = "HGVSp_Short" if "HGVSp_Short" in sample_renamed else (
        "HGVSp" if "HGVSp" in sample_renamed else None
    )

    variants: list[Variant] = []
    for raw_row in rows:
        row = _remap_row(raw_row, rename)

        gene = (row.get("Hugo_Symbol") or "").strip()
        if not gene:
            continue

        hgvsp = ""
        if hgvsp_col:
            raw = row.get(hgvsp_col) or ""
            hgvsp = raw.strip()

        var_class = (row.get("Variant_Classification") or "").strip()
        chrom = (row.get("Chromosome") or "").strip()

        start_pos = row.get("Start_Position")
        if start_pos is not None:
            try:
                start_pos = int(start_pos)
            except (ValueError, TypeError):
                start_pos = None

        ref_allele = (row.get("Reference_Allele") or "").strip()
        tumor_allele = (row.get("Tumor_Seq_Allele2") or "").strip()

        t_alt = row.get("t_alt_count")
        t_ref = row.get("t_ref_count")
        try:
            t_alt = int(t_alt) if t_alt else None
        except (ValueError, TypeError):
            t_alt = None
        try:
            t_ref = int(t_ref) if t_ref else None
        except (ValueError, TypeError):
            t_ref = None

        variant = Variant(
            gene=gene,
            hgvsp=hgvsp,
            variant_classification=var_class,
            chromosome=chrom,
            start_position=start_pos,
            reference_allele=ref_allele,
            tumor_allele=tumor_allele,
            t_alt_count=t_alt,
            t_ref_count=t_ref,
        )
        variants.append(variant)

    return variants


def parse_seg(source: str | Path) -> list[CopyNumberSegment]:
    """Parse a SEG file (path or string content) into a list of CopyNumberSegments."""
    rows = _read_tsv(source)
    if not rows:
        return []

    rename = _normalize_columns(list(rows[0].keys()), _SEG_COLUMN_MAP)

    segments: list[CopyNumberSegment] = []
    for raw_row in rows:
        row = _remap_row(raw_row, rename)

        chrom_raw = row.get("chrom") or ""
        chrom = chrom_raw.replace("chr", "").strip()

        try:
            start = int(row.get("loc.start") or 0)
            end = int(row.get("loc.end") or 0)
        except (ValueError, TypeError):
            continue

        try:
            seg_mean = float(row.get("seg.mean") or 0.0)
        except (ValueError, TypeError):
            seg_mean = 0.0

        try:
            num_marks = int(row.get("num.mark") or 0)
        except (ValueError, TypeError):
            num_marks = 0

        sample_id = (row.get("ID") or "").strip()

        segments.append(
            CopyNumberSegment(
                sample_id=sample_id,
                chrom=chrom,
                start=start,
                end=end,
                num_marks=num_marks,
                seg_mean=seg_mean,
            )
        )

    return segments
