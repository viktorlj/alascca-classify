"""Pydantic data models for ALASCCA classification."""

from __future__ import annotations

from typing import Literal

from pydantic import BaseModel, Field


class Variant(BaseModel):
    """A somatic variant from a MAF file."""

    gene: str
    hgvsp: str = ""
    variant_classification: str = ""
    chromosome: str = ""
    start_position: int | None = None
    reference_allele: str = ""
    tumor_allele: str = ""
    t_alt_count: int | None = None
    t_ref_count: int | None = None
    vaf: float | None = None
    exon: int | None = None
    domain: str = ""

    def compute_vaf(self) -> float | None:
        if self.vaf is not None:
            return self.vaf
        if self.t_alt_count is not None and self.t_ref_count is not None:
            total = self.t_alt_count + self.t_ref_count
            if total > 0:
                return self.t_alt_count / total
        return None


class CopyNumberSegment(BaseModel):
    """A segment from a SEG file."""

    sample_id: str = ""
    chrom: str
    start: int
    end: int
    num_marks: int = 0
    seg_mean: float = 0.0


class PTENDeletion(BaseModel):
    """Result of PTEN deep deletion assessment."""

    detected: bool = False
    seg_mean_pten: float | None = None
    seg_mean_upstream: float | None = None
    seg_mean_downstream: float | None = None
    log_ratio_diff_upstream: float | None = None
    log_ratio_diff_downstream: float | None = None
    passes_threshold: bool = False
    note: str = ""


class SampleMetadata(BaseModel):
    """Optional metadata for a sample."""

    sample_id: str = "UNKNOWN"
    tumor_purity: float | None = None
    ptnm_stage: str = ""
    tumor_location: str = ""
    msi_status: str = ""


class InformationalBiomarkers(BaseModel):
    """Informational biomarkers (not part of Group A/B classification)."""

    msi_status: str = ""
    braf_v600e: bool = False
    braf_variants: list[str] = Field(default_factory=list)
    kras_mutant: bool = False
    kras_variants: list[str] = Field(default_factory=list)
    nras_mutant: bool = False
    nras_variants: list[str] = Field(default_factory=list)


class ClassificationResult(BaseModel):
    """Final ALASCCA classification result for a sample."""

    sample_id: str = "UNKNOWN"
    alascca_group: Literal["A", "B", "none"] = "none"
    group_a_variants: list[Variant] = Field(default_factory=list)
    group_b_variants: list[Variant] = Field(default_factory=list)
    has_group_a: bool = False
    has_group_b: bool = False
    aspirin_eligible: bool = False
    informational_biomarkers: InformationalBiomarkers = Field(
        default_factory=InformationalBiomarkers
    )
    pten_deletion: PTENDeletion | None = None
    clinical_notes: list[str] = Field(default_factory=list)
    flags: list[str] = Field(default_factory=list)
    version: str = "0.1.0"


class SampleInput(BaseModel):
    """Input bundle for classification: variants + optional metadata and CNA segments."""

    metadata: SampleMetadata = Field(default_factory=SampleMetadata)
    variants: list[Variant] = Field(default_factory=list)
    segments: list[CopyNumberSegment] = Field(default_factory=list)
