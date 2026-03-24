"""Output formatting: JSON, TSV, human-readable text, and PDF reports."""

from __future__ import annotations

import json
from io import BytesIO

from .models import ClassificationResult


def result_to_dict(result: ClassificationResult) -> dict:
    """Convert ClassificationResult to a JSON-serializable dict."""
    return result.model_dump(mode="json")


def result_to_json(result: ClassificationResult, indent: int = 2) -> str:
    """Render classification result as formatted JSON."""
    return json.dumps(result_to_dict(result), indent=indent)


def result_to_human(result: ClassificationResult) -> str:
    """Render classification result as human-readable text."""
    lines: list[str] = []
    lines.append("=" * 60)
    lines.append("ALASCCA PI3K Pathway Classification Report")
    lines.append("=" * 60)
    lines.append(f"Sample ID:     {result.sample_id}")
    lines.append(f"ALASCCA Group: {_format_group(result.alascca_group)}")
    lines.append(f"Aspirin Eligible: {'Yes' if result.aspirin_eligible else 'No'}")
    lines.append("")

    if result.group_a_variants:
        lines.append("--- Group A Variants (PIK3CA hotspots) ---")
        for v in result.group_a_variants:
            vaf_str = f"  VAF={v.vaf:.2f}" if v.vaf is not None else ""
            lines.append(f"  {v.gene} {v.hgvsp} ({v.variant_classification}){vaf_str}")
        lines.append("")

    if result.group_b_variants:
        lines.append("--- Group B Variants ---")
        for v in result.group_b_variants:
            vaf_str = f"  VAF={v.vaf:.2f}" if v.vaf is not None else ""
            lines.append(f"  {v.gene} {v.hgvsp} ({v.variant_classification}){vaf_str}")
        lines.append("")

    if result.pten_deletion and result.pten_deletion.detected:
        lines.append("--- PTEN Deep Deletion ---")
        lines.append(f"  {result.pten_deletion.note}")
        lines.append("")

    bio = result.informational_biomarkers
    lines.append("--- Informational Biomarkers ---")
    lines.append(f"  MSI status: {bio.msi_status or 'not provided'}")
    lines.append(f"  BRAF V600E: {'Yes' if bio.braf_v600e else 'No'}")
    kras_str = ", ".join(bio.kras_variants) if bio.kras_variants else "None"
    lines.append(f"  KRAS:       {'Mutant' if bio.kras_mutant else 'Wild-type'} ({kras_str})")
    nras_str = ", ".join(bio.nras_variants) if bio.nras_variants else "None"
    lines.append(f"  NRAS:       {'Mutant' if bio.nras_mutant else 'Wild-type'} ({nras_str})")
    lines.append("")

    lines.append("--- Clinical Notes ---")
    for note in result.clinical_notes:
        lines.append(f"  - {note}")
    lines.append("")

    if result.flags:
        lines.append("--- Flags ---")
        for flag in result.flags:
            lines.append(f"  ! {flag}")
        lines.append("")

    lines.append(f"Version: {result.version}")
    lines.append("=" * 60)
    return "\n".join(lines)


def _format_group(group: str) -> str:
    if group == "A":
        return "Group A (PIK3CA hotspot)"
    elif group == "B":
        return "Group B (other PI3K pathway alteration)"
    return "No PI3K pathway alteration"


def results_to_tsv(results: list[ClassificationResult]) -> str:
    """Render batch results as TSV."""
    header = [
        "sample_id",
        "alascca_group",
        "aspirin_eligible",
        "group_a_variants",
        "group_b_variants",
        "braf_v600e",
        "kras_mutant",
        "kras_variants",
        "nras_mutant",
        "msi_status",
    ]
    lines = ["\t".join(header)]
    for r in results:
        a_vars = "; ".join(f"{v.gene} {v.hgvsp}" for v in r.group_a_variants)
        b_vars = "; ".join(f"{v.gene} {v.hgvsp}" for v in r.group_b_variants)
        bio = r.informational_biomarkers
        row = [
            r.sample_id,
            r.alascca_group,
            str(r.aspirin_eligible),
            a_vars,
            b_vars,
            str(bio.braf_v600e),
            str(bio.kras_mutant),
            "; ".join(bio.kras_variants),
            str(bio.nras_mutant),
            bio.msi_status,
        ]
        lines.append("\t".join(row))
    return "\n".join(lines)


def _sanitize_text(text: str) -> str:
    """Replace problematic Unicode characters for PDF rendering."""
    replacements = {
        "\u2013": "-",  # en-dash
        "\u2014": "--",  # em-dash
        "\u2018": "'",
        "\u2019": "'",
        "\u201c": '"',
        "\u201d": '"',
        "\u2265": ">=",
        "\u2264": "<=",
    }
    for orig, repl in replacements.items():
        text = text.replace(orig, repl)
    return text


GROUP_COLORS = {
    "A": (0, 128, 0),      # green
    "B": (0, 100, 180),    # blue
    "none": (128, 128, 128),  # gray
}


def result_to_pdf(result: ClassificationResult) -> bytes:
    """Generate a PDF report and return as bytes."""
    from fpdf import FPDF

    pdf = FPDF()
    pdf.set_auto_page_break(auto=True, margin=20)
    pdf.add_page()

    # Title
    pdf.set_font("Helvetica", "B", 16)
    pdf.cell(0, 10, "ALASCCA PI3K Pathway Classification Report", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(4)

    # Sample info
    pdf.set_font("Helvetica", "", 11)
    pdf.cell(0, 7, f"Sample ID: {result.sample_id}", new_x="LMARGIN", new_y="NEXT")

    # Group badge
    color = GROUP_COLORS.get(result.alascca_group, (128, 128, 128))
    pdf.set_font("Helvetica", "B", 14)
    pdf.set_text_color(*color)
    group_text = f"ALASCCA Group: {_format_group(result.alascca_group)}"
    pdf.cell(0, 10, group_text, new_x="LMARGIN", new_y="NEXT")
    pdf.set_text_color(0, 0, 0)
    pdf.ln(2)

    eligibility = "ELIGIBLE" if result.aspirin_eligible else "NOT ELIGIBLE"
    pdf.set_font("Helvetica", "B", 12)
    pdf.cell(0, 8, f"Aspirin Eligibility: {eligibility}", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(4)

    # Variants
    pdf.set_font("Helvetica", "B", 12)
    if result.group_a_variants:
        pdf.cell(0, 8, "Group A Variants (PIK3CA hotspots)", new_x="LMARGIN", new_y="NEXT")
        pdf.set_font("Helvetica", "", 10)
        for v in result.group_a_variants:
            vaf_str = f"  VAF={v.vaf:.2f}" if v.vaf is not None else ""
            line = f"  {v.gene} {v.hgvsp} ({v.variant_classification}){vaf_str}"
            pdf.cell(0, 6, _sanitize_text(line), new_x="LMARGIN", new_y="NEXT")
        pdf.ln(3)

    pdf.set_font("Helvetica", "B", 12)
    if result.group_b_variants:
        pdf.cell(0, 8, "Group B Variants", new_x="LMARGIN", new_y="NEXT")
        pdf.set_font("Helvetica", "", 10)
        for v in result.group_b_variants:
            vaf_str = f"  VAF={v.vaf:.2f}" if v.vaf is not None else ""
            line = f"  {v.gene} {v.hgvsp} ({v.variant_classification}){vaf_str}"
            pdf.cell(0, 6, _sanitize_text(line), new_x="LMARGIN", new_y="NEXT")
        pdf.ln(3)

    # Biomarkers
    pdf.set_font("Helvetica", "B", 12)
    pdf.cell(0, 8, "Informational Biomarkers", new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("Helvetica", "", 10)
    bio = result.informational_biomarkers
    msi_text = f"  MSI status: {bio.msi_status or 'not provided'}"
    pdf.cell(0, 6, msi_text, new_x="LMARGIN", new_y="NEXT")
    braf_text = f"  BRAF V600E: {'Yes' if bio.braf_v600e else 'No'}"
    pdf.cell(0, 6, braf_text, new_x="LMARGIN", new_y="NEXT")
    kras_vars = ", ".join(bio.kras_variants) or "None"
    kras_status = "Mutant" if bio.kras_mutant else "Wild-type"
    pdf.cell(0, 6, f"  KRAS: {kras_status} ({kras_vars})", new_x="LMARGIN", new_y="NEXT")
    nras_vars = ", ".join(bio.nras_variants) or "None"
    nras_status = "Mutant" if bio.nras_mutant else "Wild-type"
    pdf.cell(0, 6, f"  NRAS: {nras_status} ({nras_vars})", new_x="LMARGIN", new_y="NEXT")
    pdf.ln(3)

    # Clinical notes
    pdf.set_font("Helvetica", "B", 12)
    pdf.cell(0, 8, "Clinical Notes", new_x="LMARGIN", new_y="NEXT")
    pdf.set_font("Helvetica", "", 10)
    for note in result.clinical_notes:
        pdf.multi_cell(0, 6, _sanitize_text(f"  - {note}"))
    pdf.ln(3)

    # Flags
    if result.flags:
        pdf.set_font("Helvetica", "B", 12)
        pdf.cell(0, 8, "Flags", new_x="LMARGIN", new_y="NEXT")
        pdf.set_font("Helvetica", "", 10)
        for flag in result.flags:
            pdf.multi_cell(0, 6, _sanitize_text(f"  ! {flag}"))
        pdf.ln(3)

    # Disclaimer
    pdf.ln(6)
    pdf.set_font("Helvetica", "I", 8)
    pdf.multi_cell(0, 4, _sanitize_text(
        "Disclaimer: This report is for research and clinical decision support only. "
        "Final treatment decisions must be made by the treating physician. "
        "Based on Martling et al., NEJM 2025;393:1051-64. "
        f"alascca-classify v{result.version}"
    ))

    buf = BytesIO()
    pdf.output(buf)
    return buf.getvalue()
