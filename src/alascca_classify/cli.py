"""CLI entrypoints using typer."""

from __future__ import annotations

import json
from pathlib import Path

import typer
from rich.console import Console
from rich.table import Table

from . import __version__
from .classify import classify_sample
from .io import parse_maf, parse_seg
from .models import SampleInput, SampleMetadata
from .pik3ca import (
    GROUP_A_POSITIONS,
    POSITION_ANNOTATIONS,
    extract_protein_position,
    map_consequence,
)
from .report import result_to_human, result_to_json, result_to_pdf, results_to_tsv

app = typer.Typer(
    name="alascca-classify",
    help="Classify PI3K pathway alterations in colorectal cancer per the ALASCCA trial.",
)
console = Console()


@app.command()
def classify(
    maf_file: Path = typer.Argument(..., help="Path to MAF file"),
    seg_file: Path | None = typer.Option(None, "--seg", help="Path to SEG file for PTEN deletion"),
    metadata_file: Path | None = typer.Option(
        None, "--metadata", help="Path to JSON metadata file"
    ),
    sample_id: str = typer.Option("", "--sample-id", help="Sample identifier"),
    output: Path | None = typer.Option(None, "--output", "-o", help="Output file path"),
    human: bool = typer.Option(False, "--human", help="Human-readable output"),
    pdf: bool = typer.Option(False, "--pdf", help="Generate PDF report"),
) -> None:
    """Classify a single sample from a MAF file."""
    if not maf_file.exists():
        console.print(f"[red]Error:[/red] MAF file not found: {maf_file}")
        raise typer.Exit(1)

    variants = parse_maf(maf_file)

    # Parse SEG if provided
    segments = []
    if seg_file:
        if not seg_file.exists():
            console.print(f"[red]Error:[/red] SEG file not found: {seg_file}")
            raise typer.Exit(1)
        segments = parse_seg(seg_file)

    # Load metadata
    metadata = SampleMetadata()
    if metadata_file and metadata_file.exists():
        with open(metadata_file) as f:
            meta_data = json.load(f)
        metadata = SampleMetadata(**meta_data)
    if sample_id:
        metadata.sample_id = sample_id
    elif metadata.sample_id == "UNKNOWN":
        metadata.sample_id = maf_file.stem

    sample = SampleInput(metadata=metadata, variants=variants, segments=segments)
    result = classify_sample(sample)

    if pdf:
        pdf_bytes = result_to_pdf(result)
        out_path = output or Path(f"{result.sample_id}_alascca.pdf")
        out_path.write_bytes(pdf_bytes)
        console.print(f"PDF report written to {out_path}")
    elif human:
        text = result_to_human(result)
        if output:
            output.write_text(text)
        else:
            console.print(text)
    else:
        text = result_to_json(result)
        if output:
            output.write_text(text)
        else:
            console.print(text)


@app.command("classify-batch")
def classify_batch(
    input_dir: Path = typer.Option(..., "--input-dir", help="Directory containing MAF files"),
    output: Path = typer.Option("results.tsv", "--output", "-o", help="Output TSV file"),
    metadata: Path | None = typer.Option(
        None, "--metadata", help="Metadata TSV with sample_id column"
    ),
    seg_dir: Path | None = typer.Option(None, "--seg-dir", help="Directory containing SEG files"),
) -> None:
    """Classify a batch of samples from a directory of MAF files."""
    if not input_dir.is_dir():
        console.print(f"[red]Error:[/red] Input directory not found: {input_dir}")
        raise typer.Exit(1)

    maf_files = sorted(input_dir.glob("*.maf"))
    if not maf_files:
        console.print(f"[yellow]Warning:[/yellow] No MAF files found in {input_dir}")
        raise typer.Exit(0)

    # Load metadata if provided
    meta_map: dict[str, dict] = {}
    if metadata and metadata.exists():
        import polars as pl

        meta_df = pl.read_csv(metadata, separator="\t")
        for row in meta_df.iter_rows(named=True):
            sid = str(row.get("sample_id", ""))
            if sid:
                meta_map[sid] = row

    results = []
    for maf_path in maf_files:
        sample_id = maf_path.stem
        variants = parse_maf(maf_path)

        segments = []
        if seg_dir:
            seg_path = seg_dir / f"{sample_id}.seg"
            if seg_path.exists():
                segments = parse_seg(seg_path)

        meta = SampleMetadata(sample_id=sample_id)
        if sample_id in meta_map:
            meta = SampleMetadata(**{k: v for k, v in meta_map[sample_id].items() if v is not None})

        sample = SampleInput(metadata=meta, variants=variants, segments=segments)
        result = classify_sample(sample)
        results.append(result)

        group_color = {"A": "green", "B": "blue", "none": "dim"}.get(result.alascca_group, "white")
        console.print(
            f"  {sample_id}: [{group_color}]Group {result.alascca_group}[/{group_color}]"
        )

    tsv = results_to_tsv(results)
    output.write_text(tsv)
    console.print(f"\n[green]Results written to {output}[/green] ({len(results)} samples)")


@app.command("check-variant")
def check_variant(
    gene: str = typer.Argument(..., help="Gene symbol (PIK3CA, PIK3R1, PTEN)"),
    hgvsp: str = typer.Argument(..., help="Protein change (e.g., p.E545K)"),
    variant_classification: str = typer.Option(
        "Missense_Mutation", "--type", help="MAF Variant_Classification"
    ),
) -> None:
    """Check the ALASCCA classification of a specific variant."""
    from .models import Variant
    from .pik3ca import is_group_a, is_group_b_pik3ca
    from .pik3r1 import is_group_b_pik3r1
    from .pten import is_group_b_pten

    variant = Variant(
        gene=gene, hgvsp=hgvsp, variant_classification=variant_classification
    )

    vep = map_consequence(variant_classification)
    pos = extract_protein_position(hgvsp)

    table = Table(title=f"ALASCCA Classification: {gene} {hgvsp}")
    table.add_column("Property", style="bold")
    table.add_column("Value")

    table.add_row("Gene", gene.upper())
    table.add_row("Protein change", hgvsp)
    table.add_row("MAF classification", variant_classification)
    table.add_row("VEP consequence", vep or "(unmapped)")
    table.add_row("Protein position", str(pos) if pos else "(not parsed)")

    if gene.upper() == "PIK3CA" and pos:
        is_hotspot = pos in GROUP_A_POSITIONS
        table.add_row("PIK3CA hotspot", "[green]Yes[/green]" if is_hotspot else "No")
        if is_hotspot and str(pos) in POSITION_ANNOTATIONS:
            ann = POSITION_ANNOTATIONS[str(pos)]
            table.add_row("Exon", str(ann["exon"]))
            table.add_row("Domain", ann["domain"])

    group_a = is_group_a(variant)
    group_b = (
        is_group_b_pik3ca(variant)
        or is_group_b_pik3r1(variant)
        or is_group_b_pten(variant)
    )

    if group_a:
        table.add_row("ALASCCA Group", "[green bold]A[/green bold]")
    elif group_b:
        table.add_row("ALASCCA Group", "[blue bold]B[/blue bold]")
    else:
        table.add_row("ALASCCA Group", "[dim]none[/dim]")

    eligible_str = (
        "[green]Yes[/green]" if (group_a or group_b) else "[dim]No[/dim]"
    )
    table.add_row("Aspirin eligible", eligible_str)

    console.print(table)


@app.command()
def serve(
    port: int = typer.Option(8000, "--port", "-p", help="Port to serve on"),
    host: str = typer.Option("127.0.0.1", "--host", help="Host to bind to"),
) -> None:
    """Start the web server."""
    import uvicorn

    from .web.app import app as web_app

    console.print(f"Starting ALASCCA-classify web server on http://{host}:{port}")
    uvicorn.run(web_app, host=host, port=port)


@app.command()
def version() -> None:
    """Show version."""
    console.print(f"alascca-classify v{__version__}")


if __name__ == "__main__":
    app()
