# CLAUDE.md — alascca-classify

## Project overview

A Python CLI and web app for **classifying PI3K pathway alterations in colorectal cancer** to determine eligibility for adjuvant aspirin therapy per the **ALASCCA trial** (Martling et al., NEJM 2025; DOI: 10.1056/NEJMoa2504650).

The tool ingests MAF files (and optionally SEG files for copy number) from targeted panel sequencing of colorectal cancer, classifies somatic alterations in PIK3CA, PIK3R1, and PTEN into Group A and Group B categories, and produces a structured clinical report indicating whether the patient's tumor harbors alterations that predict benefit from adjuvant low-dose aspirin (160 mg daily for 3 years).

## Clinical background

### The ALASCCA trial

ALASCCA was a double-blind, randomized, placebo-controlled trial of adjuvant aspirin (160 mg/day × 3 years) in resected stage I–III colorectal cancer with somatic PI3K pathway alterations. Key findings:

- **Group A patients** (PIK3CA hotspot mutations in exon 9/20): 3-year recurrence 7.7% aspirin vs 14.1% placebo (HR 0.49, p=0.04)
- **Group B patients** (other PI3K pathway alterations): 3-year recurrence 7.7% aspirin vs 16.8% placebo (HR 0.42)
- Combined: HR 0.45 (95% CI 0.28–0.74)
- ~37% of CRC patients have PI3K pathway alterations (17.3% Group A, 19.7% Group B)

### Mutation classification (from ALASCCA Supplementary Table S3)

**Group A — PIK3CA exon 9/20 hotspot mutations:**

Missense variants at the following protein positions in PIK3CA (ENSG00000121879, transcript ENST00000263967):
- Exon 9 (helical domain): **E542**, **E545**, **Q546**
- Exon 20 (kinase domain): **H1021**, **H1043**, **H1044**, **H1047**

These are the canonical PIK3CA gain-of-function hotspots. Any missense at these 7 positions → Group A.

**Group B — Other moderate/high-impact PI3K pathway variants:**

Three genes, four alteration types:

| Gene | Transcript | Consequence types |
|------|-----------|-------------------|
| PIK3CA | ENST00000263967 | missense_variant (non-hotspot), frameshift_variant, inframe_insertion, inframe_deletion, start_lost, stop_gained, splice_acceptor_variant, splice_donor_variant |
| PIK3R1 | ENST00000521381 | missense_variant, frameshift_variant, inframe_insertion, inframe_deletion, start_lost, stop_gained, splice_acceptor_variant, splice_donor_variant |
| PTEN | ENST00000371953 | missense_variant, frameshift_variant, inframe_insertion, inframe_deletion, start_lost, stop_gained, splice_acceptor_variant, splice_donor_variant, **homozygous_loss** (deep deletion) |

**Priority rule:** If a patient has BOTH Group A and Group B alterations, they are assigned to **Group A**.

### PTEN deep deletion detection (requires SEG file)

From the ALASCCA supplementary methods:
- Requires estimated tumor cell fraction ≥ 0.30
- A deletion is called if a CNV segment overlapping PTEN (including at least one exon) shows at least 15% lower DNA abundance (absolute log-ratio difference of 0.2) compared to flanking regions both 3 Mb upstream and downstream of PTEN
- This is the only CNV-based alteration in the classification

### Additional biomarkers (informational, not part of classification)

The ALASCCA genomic report also included:
- **MSI status** (MSI-H, MSS, not determined) — informational for clinical context
- **BRAF V600E** mutation status — informational
- **KRAS** mutation status — informational (note: ~50% of randomized patients were KRAS-mutant)
- **NRAS** mutation status — informational

These should be detected and reported but do NOT affect the Group A/B classification.

### Clinical context for the report

Key numbers-needed-to-treat from the trial (combined cohort):
- pTNM stage III rectal cancer: NNT = 6
- pTNM stage III colon cancer: NNT = 9
- pTNM stage II rectal cancer: NNT = 12
- pTNM stage I rectal cancer: NNT = 22
- pTNM stage II colon cancer: NNT = 42

The report should note that benefit was observed across stages but absolute benefit is greatest in higher-stage disease.

## Technical conventions

### Language & tooling

- **Python only** (≥3.11)
- **uv** for package management (NOT pip, NOT conda)
- **Marimo** notebooks exclusively (NEVER Jupyter). Follow Marimo conventions: functional reactive cells, `mo.md()` for markdown, no global mutable state.
- **polars** for dataframes (NOT pandas)
- **pydantic** for data models and validation
- **typer** for CLI
- **FastAPI** for web backend
- **pytest** for testing
- **ruff** for linting/formatting

### Project structure

```
alascca-classify/
├── CLAUDE.md                    # this file
├── pyproject.toml               # uv project, all deps here
├── Dockerfile                   # for Cloud Run deployment
├── .dockerignore
├── .gitignore
├── LICENSE                      # MIT
├── README.md
├── logbook.md                   # experiment log (append-only)
├── findings.md                  # key findings and decisions
├── src/
│   └── alascca_classify/
│       ├── __init__.py
│       ├── models.py            # pydantic: SampleInput, ClassificationResult, Variant, etc.
│       ├── classify.py          # main orchestrator: Group A/B classification
│       ├── pik3ca.py            # PIK3CA hotspot lookup + non-hotspot variant assessment
│       ├── pik3r1.py            # PIK3R1 variant assessment
│       ├── pten.py              # PTEN variant + deep deletion assessment
│       ├── biomarkers.py        # MSI, BRAF, KRAS, NRAS detection (informational)
│       ├── cna.py               # SEG file parsing + PTEN deep deletion logic
│       ├── io.py                # MAF parsing (file and in-memory), SEG parsing
│       ├── report.py            # JSON/TSV/human-readable output + PDF generation
│       ├── cli.py               # typer CLI entrypoints
│       ├── data/
│       │   ├── pik3ca_hotspots.json     # Group A positions with exon/domain annotation
│       │   ├── consequence_types.json   # VEP consequence terms for Group B
│       │   ├── gene_transcripts.json    # canonical transcript IDs per gene
│       │   └── pten_coordinates.json    # PTEN exon coordinates + flanking region definitions
│       └── web/
│           ├── app.py           # FastAPI application
│           ├── templates/       # Jinja2 templates
│           │   ├── base.html
│           │   ├── classify.html
│           │   ├── result.html
│           │   └── methods.html
│           └── static/
│               └── style.css
├── demo/
│   ├── README.md
│   ├── group_a_pik3ca_e545k.maf        # clear Group A case
│   ├── group_b_pten_frameshift.maf     # Group B: PTEN truncating
│   ├── group_b_pik3r1_missense.maf     # Group B: PIK3R1
│   ├── group_ab_combined.maf           # both Group A + B → classified as A
│   ├── no_pi3k_alteration.maf          # wild-type for PI3K pathway
│   ├── pten_deletion.seg               # SEG file with PTEN deep deletion
│   └── metadata.tsv                    # sample metadata for batch mode
├── scripts/
│   ├── validate_genie.py        # validate against AACR GENIE data
│   └── generate_demo_data.py    # create synthetic demo MAF/SEG files
├── notebooks/
│   ├── 01_explore_classification.py    # marimo: interactive walkthrough
│   └── 02_validate_cohort.py           # marimo: validate on real data
└── tests/
    ├── conftest.py              # shared fixtures
    ├── test_classify.py         # end-to-end classification
    ├── test_pik3ca.py           # PIK3CA hotspot + non-hotspot
    ├── test_pik3r1.py           # PIK3R1 assessment
    ├── test_pten.py             # PTEN variant + deletion
    ├── test_cna.py              # SEG file + PTEN deletion logic
    ├── test_biomarkers.py       # informational biomarker detection
    ├── test_io.py               # MAF/SEG parsing
    └── data/                    # test fixture files
```

### CLI interface

```bash
# Classify a single sample
alascca-classify classify sample.maf --output result.json

# With SEG file for PTEN deletion
alascca-classify classify sample.maf --seg sample.seg --output result.json

# Human-readable output
alascca-classify classify sample.maf --human

# Classify a batch
alascca-classify classify-batch --input-dir ./samples/ --output results.tsv

# Check a specific variant
alascca-classify check-variant PIK3CA "p.E545K"

# Start web server
alascca-classify serve --port 8000
```

### Input format

**MAF file** (minimum required columns):

| Column | Required | Description |
|--------|----------|-------------|
| Hugo_Symbol | Yes | Gene symbol (PIK3CA, PIK3R1, PTEN, etc.) |
| Variant_Classification | Yes | VEP consequence (Missense_Mutation, Frame_Shift_Del, etc.) |
| HGVSp_Short | Yes | Protein change (e.g., p.E545K) |
| Chromosome | No | For positional validation |
| Start_Position | No | Genomic coordinate |
| Reference_Allele | No | Reference allele |
| Tumor_Seq_Allele2 | No | Alternate allele |
| t_alt_count | No | Tumor alt read count (for VAF) |
| t_ref_count | No | Tumor ref read count (for VAF) |

**SEG file** (for PTEN deletion, tab-delimited):

| Column | Description |
|--------|-------------|
| ID / Sample | Sample identifier |
| chrom | Chromosome |
| loc.start | Segment start |
| loc.end | Segment end |
| num.mark | Number of probes |
| seg.mean | Log2 ratio |

**Sample metadata** (JSON sidecar, optional):

```json
{
  "sample_id": "SAMPLE_001",
  "tumor_purity": 0.45,
  "ptnm_stage": "III",
  "tumor_location": "colon_right",
  "msi_status": "MSS"
}
```
### Example data
The AACR Genie V18 dataset can be found here: /Users/viklj600/Bioinformatik/Data/Genie_V18

For some panels there are both mutation data in MAF format and CNA data in .seg-format that can be used to test out the algorithm.

### Output format

JSON per sample:

```json
{
  "sample_id": "SAMPLE_001",
  "alascca_group": "A",
  "group_a_variants": [
    {
      "gene": "PIK3CA",
      "hgvsp": "p.E545K",
      "exon": 9,
      "domain": "helical",
      "variant_classification": "Missense_Mutation",
      "vaf": 0.35
    }
  ],
  "group_b_variants": [],
  "has_group_a": true,
  "has_group_b": false,
  "aspirin_eligible": true,
  "informational_biomarkers": {
    "msi_status": "MSS",
    "braf_v600e": false,
    "kras_mutant": true,
    "kras_variants": ["p.G12D"],
    "nras_mutant": false
  },
  "pten_deletion": null,
  "clinical_notes": [
    "PIK3CA E545K is a Group A hotspot mutation (exon 9, helical domain).",
    "This tumor is eligible for adjuvant aspirin per ALASCCA (HR 0.49 for recurrence).",
    "ALASCCA trial: 160 mg aspirin daily for 3 years post-resection."
  ],
  "version": "0.1.0"
}
```

## Web app design

### Architecture

- FastAPI backend with Jinja2 templates (same pattern as ec-molsubtype)
- Single-page classify form with file upload (MAF required, SEG optional)
- Optional metadata fields (sample ID, stage, location, MSI status)
- Result page with classification, variant table, and clinical interpretation
- Methods page documenting the ALASCCA trial and classification algorithm
- Demo samples available on the landing page
- PDF report export

### Cloud Run deployment

- Dockerfile: python:3.12-slim base, uv install, gunicorn + uvicorn workers
- PORT from environment variable (Cloud Run convention)
- Stateless — no persistent storage needed
- Health check endpoint at /health

### Pages

1. **Classify** (landing page): MAF upload, optional SEG upload, optional metadata fields, demo sample buttons
2. **Result**: classification badge (Group A / Group B / No PI3K alteration), variant table, informational biomarkers, clinical notes, PDF download
3. **Methods**: ALASCCA trial summary, classification algorithm description, Group A/B definitions, PTEN deletion criteria, limitations and caveats

## Development workflow

1. `models.py` — pydantic models for all data structures
2. `io.py` — MAF and SEG file parsing with polars
3. `pik3ca.py` — Group A hotspot detection + Group B non-hotspot variants
4. `pik3r1.py` — Group B PIK3R1 variant assessment
5. `pten.py` — Group B PTEN variant assessment
6. `cna.py` — SEG file parsing + PTEN deep deletion logic
7. `biomarkers.py` — MSI/BRAF/KRAS/NRAS detection
8. `classify.py` — orchestrator combining all modules
9. `report.py` — output formatting
10. `cli.py` — typer CLI
11. `web/` — FastAPI web interface
12. `demo/` — synthetic demo samples
13. Tests throughout (>90% coverage on classification logic)
14. Dockerfile + Cloud Run configuration

## Variant classification mapping

MAF `Variant_Classification` values → ALASCCA consequence types:

| MAF Variant_Classification | ALASCCA consequence |
|---------------------------|-------------------|
| Missense_Mutation | missense_variant |
| Frame_Shift_Del | frameshift_variant |
| Frame_Shift_Ins | frameshift_variant |
| In_Frame_Del | inframe_deletion |
| In_Frame_Ins | inframe_insertion |
| Nonsense_Mutation | stop_gained |
| Splice_Site | splice_acceptor_variant or splice_donor_variant |
| Translation_Start_Site | start_lost |

## Important caveats to document in the tool

1. **Research and clinical decision support only.** Final treatment decisions must be made by the treating physician with full clinical context.
2. **PTEN deep deletion detection from SEG files** requires adequate tumor purity (≥30%) and uses a simplified algorithm compared to the ALASCCA trial's full CNVkit pipeline.
3. **The tool classifies variants but does not assess pathogenicity.** All missense variants in PIK3R1 and PTEN are included per the ALASCCA protocol, which used VEP "moderate impact" as the threshold.
4. **Germline filtering is assumed.** The input MAF should contain somatic variants only (tumor-normal paired calling). The ALASCCA trial used paired tumor-blood sequencing.
5. **This tool does not replace the full ALASCCA genomic analysis pipeline** which included custom library prep, targeted sequencing at ~30M read pairs, and manual QC review of all results.

## Key references

- Martling A, et al. Low-Dose Aspirin for PI3K-Altered Localized Colorectal Cancer. N Engl J Med 2025;393:1051-64.
- ALASCCA Supplementary Appendix (Table S3: somatic alteration definitions)
- Cathomas G. PIK3CA in colorectal cancer. Front Oncol 2014;4:35.
- Hall DCN, Benndorf RA. Aspirin sensitivity of PIK3CA-mutated CRC. Cell Mol Life Sci 2022;79:393.
- Liao X, et al. Aspirin use, tumor PIK3CA mutation, and colorectal-cancer survival. N Engl J Med 2012;367:1596-606.
- Yang J, et al. Aspirin prevents metastasis by limiting platelet TXA2 suppression of T cell immunity. Nature 2025;640:1052-61.
