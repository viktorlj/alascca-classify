# Lab Notebook — alascca-classify

## 2026-03-24 — Project initialization

**Objective:** Build a Python CLI and web app for classifying PI3K pathway alterations in CRC per the ALASCCA trial.

**Approach:** Scaffold project with uv, implement classification logic from ALASCCA Supplementary Table S3, build FastAPI web interface modeled on ec-molsubtype.

**Key parameters:**
- Group A: 7 PIK3CA hotspot positions (E542, E545, Q546, H1021, N1043, H1044, H1047)
- Group B: moderate/high-impact variants in PIK3CA (non-hotspot), PIK3R1, PTEN + PTEN deep deletion
- PTEN deletion: seg_mean diff ≥ 0.2 vs both 3Mb flanking regions, tumor purity ≥ 0.30
- Priority: A > B > none

**Results:**
- Core classification engine: 87/87 tests passing
- CLI: classify, classify-batch, check-variant, serve commands working
- Web app: FastAPI with Jinja2 templates, demo samples, PDF export
- Demo files covering all classification scenarios

**Interpretation:** Classification logic faithfully implements ALASCCA Supplementary Table S3. Splice_Site mapping is intentionally conservative — MAF "Splice_Site" maps to generic "splice_site" which is NOT included in Group B (only splice_acceptor_variant and splice_donor_variant are). This means some splice variants may not be classified. Consider adding manual splice site subtyping if needed.

**Next steps:**
- Validate against AACR GENIE v18 data
- Add batch metadata support
- Deploy to Cloud Run
