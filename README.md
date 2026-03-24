# alascca-classify

Classification of PI3K pathway alterations in colorectal cancer to determine eligibility for adjuvant aspirin therapy per the **ALASCCA trial** (Martling et al., NEJM 2025).

## What it does

Ingests MAF files (and optionally SEG files for copy number) from targeted panel sequencing of colorectal cancer and classifies somatic alterations in **PIK3CA**, **PIK3R1**, and **PTEN** into Group A and Group B categories.

| Group | Definition | Aspirin HR |
|-------|-----------|-----------|
| **A** | PIK3CA exon 9/20 hotspot missense (7 positions) | 0.49 |
| **B** | Other moderate/high-impact variants in PIK3CA, PIK3R1, PTEN; PTEN deep deletion | 0.42 |
| none | No qualifying PI3K pathway alterations | — |

Priority: A + B → A. Approximately 37% of CRC patients have PI3K pathway alterations.

## Installation

```bash
# Clone and install with uv
git clone https://github.com/viktorlj/alascca-classify.git
cd alascca-classify
uv venv && source .venv/bin/activate
uv sync
```

## Usage

### Web interface

```bash
alascca-classify serve --port 8000
# Open http://localhost:8000
```

### CLI — Single sample

```bash
# JSON output
alascca-classify classify sample.maf

# With SEG file for PTEN deletion
alascca-classify classify sample.maf --seg sample.seg

# Human-readable output
alascca-classify classify sample.maf --human

# PDF report
alascca-classify classify sample.maf --pdf
```

### CLI — Batch

```bash
alascca-classify classify-batch --input-dir ./samples/ --output results.tsv
```

### CLI — Check a variant

```bash
alascca-classify check-variant PIK3CA "p.E545K"
alascca-classify check-variant PTEN "p.R130*" --type Nonsense_Mutation
```

## Input format

**MAF file** (tab-separated, minimum columns):

| Column | Required |
|--------|----------|
| Hugo_Symbol | Yes |
| Variant_Classification | Yes |
| HGVSp_Short | Yes |

**SEG file** (tab-separated, for PTEN deletion):

| Column | Description |
|--------|-------------|
| ID | Sample identifier |
| chrom | Chromosome |
| loc.start | Segment start |
| loc.end | Segment end |
| num.mark | Number of probes |
| seg.mean | Log2 ratio |

## Group A — PIK3CA hotspot positions

| Position | Exon | Domain | Common changes |
|----------|------|--------|----------------|
| E542 | 9 | Helical | E542K, E542V |
| E545 | 9 | Helical | E545K, E545G, E545A |
| Q546 | 9 | Helical | Q546K, Q546R |
| H1021 | 20 | Kinase | H1021R |
| N1043 | 20 | Kinase | N1043K |
| H1044 | 20 | Kinase | H1044R, H1044L |
| H1047 | 20 | Kinase | H1047R, H1047L, H1047Y |

## Demo samples

| File | Classification |
|------|---------------|
| `demo/group_a_pik3ca_e545k.maf` | Group A — PIK3CA E545K |
| `demo/group_b_pten_frameshift.maf` | Group B — PTEN frameshift |
| `demo/group_b_pik3r1_missense.maf` | Group B — PIK3R1 missense |
| `demo/group_ab_combined.maf` | A + B → Group A |
| `demo/no_pi3k_alteration.maf` | No PI3K alteration |
| `demo/pten_deletion.maf` + `.seg` | Group B — PTEN deep deletion |

## Project structure

```
src/alascca_classify/
├── models.py        # Pydantic data models
├── io.py            # MAF/SEG parsing (polars)
├── classify.py      # Main orchestrator
├── pik3ca.py        # PIK3CA Group A/B
├── pik3r1.py        # PIK3R1 Group B
├── pten.py          # PTEN Group B
├── cna.py           # PTEN deep deletion from SEG
├── biomarkers.py    # BRAF/KRAS/NRAS/MSI
├── report.py        # JSON/TSV/text/PDF output
├── cli.py           # Typer CLI
├── data/            # Classification criteria JSON
└── web/             # FastAPI web interface
    ├── app.py
    ├── templates/
    └── static/
```

## Caveats

1. **Research and clinical decision support only.** Final treatment decisions must be made by the treating physician.
2. PTEN deep deletion requires tumor purity ≥ 30%.
3. All missense variants are included per ALASCCA protocol (VEP moderate impact threshold), without pathogenicity assessment.
4. Germline filtering is assumed — input should be somatic variants only.
5. This does not replace the full ALASCCA genomic pipeline.

## References

- Martling A, et al. Low-Dose Aspirin for PI3K-Altered Localized Colorectal Cancer. *N Engl J Med* 2025;393:1051-64.
- Liao X, et al. Aspirin use, tumor PIK3CA mutation, and colorectal-cancer survival. *N Engl J Med* 2012;367:1596-606.
- Yang J, et al. Aspirin prevents metastasis by limiting platelet TXA2 suppression of T cell immunity. *Nature* 2025;640:1052-61.

## License

MIT
