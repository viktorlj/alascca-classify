"""MAF and SEG file parsing with polars."""

from __future__ import annotations

from io import StringIO
from pathlib import Path

import polars as pl

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


def _normalize_columns(df: pl.DataFrame, column_map: dict[str, str]) -> pl.DataFrame:
    """Normalize column names using a mapping (case-insensitive)."""
    rename_map = {}
    for col in df.columns:
        lower = col.lower().replace(" ", "_")
        if lower in column_map:
            target = column_map[lower]
            if target not in rename_map.values():
                rename_map[col] = target
        elif col in column_map.values():
            pass  # Already normalized
    if rename_map:
        df = df.rename(rename_map)
    return df


def _skip_comment_lines(text: str) -> str:
    """Remove lines starting with # (MAF comment headers)."""
    lines = text.split("\n")
    filtered = [line for line in lines if not line.startswith("#")]
    return "\n".join(filtered)


def _read_tsv(source: str | Path) -> pl.DataFrame:
    """Read a TSV from a file path or string content."""
    is_file = isinstance(source, Path) or (
        isinstance(source, str) and "\t" not in source and "\n" not in source
    )
    if is_file:
        path = Path(source)
        text = path.read_text()
        text = _skip_comment_lines(text)
        return pl.read_csv(StringIO(text), separator="\t", infer_schema_length=10000)
    else:
        text = _skip_comment_lines(source)
        return pl.read_csv(StringIO(text), separator="\t", infer_schema_length=10000)


def parse_maf(source: str | Path) -> list[Variant]:
    """Parse a MAF file (path or string content) into a list of Variants."""
    df = _read_tsv(source)
    df = _normalize_columns(df, _MAF_COLUMN_MAP)

    # Determine protein change column
    hgvsp_col = None
    if "HGVSp_Short" in df.columns:
        hgvsp_col = "HGVSp_Short"
    elif "HGVSp" in df.columns:
        hgvsp_col = "HGVSp"

    variants: list[Variant] = []
    for row in df.iter_rows(named=True):
        gene = str(row.get("Hugo_Symbol", "")).strip()
        if not gene:
            continue

        hgvsp = ""
        if hgvsp_col:
            raw = row.get(hgvsp_col, "")
            hgvsp = str(raw).strip() if raw is not None else ""

        var_class = str(row.get("Variant_Classification", "")).strip()
        chrom = str(row.get("Chromosome", "")).strip()

        start_pos = row.get("Start_Position")
        if start_pos is not None:
            try:
                start_pos = int(start_pos)
            except (ValueError, TypeError):
                start_pos = None

        ref_allele = str(row.get("Reference_Allele", "")).strip()
        tumor_allele = str(row.get("Tumor_Seq_Allele2", "")).strip()

        t_alt = row.get("t_alt_count")
        t_ref = row.get("t_ref_count")
        try:
            t_alt = int(t_alt) if t_alt is not None else None
        except (ValueError, TypeError):
            t_alt = None
        try:
            t_ref = int(t_ref) if t_ref is not None else None
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
    df = _read_tsv(source)
    df = _normalize_columns(df, _SEG_COLUMN_MAP)

    segments: list[CopyNumberSegment] = []
    for row in df.iter_rows(named=True):
        chrom_raw = row.get("chrom", "")
        chrom = str(chrom_raw).replace("chr", "").strip()

        try:
            start = int(row.get("loc.start", 0))
            end = int(row.get("loc.end", 0))
        except (ValueError, TypeError):
            continue

        try:
            seg_mean = float(row.get("seg.mean", 0.0))
        except (ValueError, TypeError):
            seg_mean = 0.0

        try:
            num_marks = int(row.get("num.mark", 0))
        except (ValueError, TypeError):
            num_marks = 0

        sample_id = str(row.get("ID", "")).strip()

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
