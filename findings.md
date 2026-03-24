# Findings — alascca-classify

## Classification engine validated against ALASCCA criteria

**Date:** 2026-03-24
**Confidence:** High
**Evidence:** 87/87 unit tests pass covering all Group A hotspot positions, Group B consequence types for all three genes, priority rule, PTEN deep deletion, and informational biomarkers.

**Open questions:**
- How does splice site classification compare to ALASCCA trial in practice? Our mapping is conservative (generic Splice_Site from MAF does not auto-classify).
- What is the concordance with the original ALASCCA genomic pipeline on real data?

---

## Emerging Hypotheses

### H1: Splice site variants may be under-classified
**Falsifiable test:** Run against GENIE v18 CRC cohort, compare splice site variant counts in PI3K genes to expected rates from literature.
**Priority:** Medium

### H2: PTEN deep deletion threshold may need panel-specific calibration
**Falsifiable test:** Compare PTEN deletion calls across multiple panel types in GENIE (e.g., MSK-IMPACT vs DFCI-OncoPanel) to check for systematic differences.
**Priority:** Medium
