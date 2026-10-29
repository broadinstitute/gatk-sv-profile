"""Run provenance and the label-independent summary.json written by analyze."""

from __future__ import annotations

import hashlib
import json
import logging
import re
from pathlib import Path
from typing import Dict, List, Mapping, Optional

import pandas as pd

from .config import AnalysisConfig, AnalysisMode

logger = logging.getLogger(__name__)

SUMMARY_SCHEMA_VERSION = 1
PROVENANCE_FILE_NAME = "provenance.json"
SUMMARY_FILE_NAME = "summary.json"

_GATK_VERSION_PATTERN = re.compile(r"\(GATK\)\s+v?(\S+)")


def file_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for chunk in iter(lambda: handle.read(1 << 20), b""):
            digest.update(chunk)
    return digest.hexdigest()


def parse_gatk_version(raw_output: Optional[str]) -> Optional[str]:
    """Extract the version from `gatk --version` output, falling back to the trimmed text."""
    if not raw_output:
        return None
    match = _GATK_VERSION_PATTERN.search(raw_output)
    if match:
        return match.group(1)
    text = raw_output.strip()
    return text or None


def clustering_config_entry(clustering_config: Optional[Path]) -> Optional[Dict[str, str]]:
    """Describe the SVConcordance clustering config, or None when GATK's built-in defaults are used."""
    if clustering_config is None:
        return None
    return {"path": str(clustering_config.resolve()), "sha256": file_sha256(clustering_config)}


def write_preprocess_provenance(
    preprocess_dir: Path,
    *,
    config: AnalysisConfig,
    gatk_version: Optional[str],
    annotated_a: Path,
    annotated_b: Optional[Path],
) -> Path:
    """Record how preprocess produced the annotated VCFs so analyze can report it."""
    paired = config.mode == AnalysisMode.PAIRED
    payload = {
        "schema_version": SUMMARY_SCHEMA_VERSION,
        "svconcordance_run": paired,
        "clustering_config": clustering_config_entry(config.clustering_config) if paired else None,
        "gatk_version": gatk_version,
        "annotated_a": annotated_a.name,
        "annotated_b": annotated_b.name if annotated_b is not None else None,
    }
    path = preprocess_dir / PROVENANCE_FILE_NAME
    path.write_text(json.dumps(payload, indent=2) + "\n")
    return path


def load_preprocess_provenance(vcf_a_path: Optional[Path], vcf_b_path: Optional[Path] = None) -> Optional[dict]:
    """Return preprocess provenance stored next to the annotated VCF(s), if it describes these files."""
    if vcf_a_path is None:
        return None
    path = vcf_a_path.parent / PROVENANCE_FILE_NAME
    if not path.exists():
        return None
    try:
        payload = json.loads(path.read_text())
    except (OSError, ValueError):
        logger.warning("Ignoring unreadable provenance file %s", path)
        return None
    if payload.get("annotated_a") != vcf_a_path.name:
        return None
    if vcf_b_path is not None and (payload.get("annotated_b") != vcf_b_path.name or vcf_b_path.parent != vcf_a_path.parent):
        return None
    return payload


def build_labels_block(config: AnalysisConfig, provenance: Optional[Mapping[str, object]]) -> dict:
    """Describe the run itself: labels, matching config, GATK version and bucketing mode."""
    paired = config.mode == AnalysisMode.PAIRED
    return {
        "schema_version": SUMMARY_SCHEMA_VERSION,
        "mode": "paired" if paired else "single",
        "label_a": config.vcf_a_label,
        "label_b": config.vcf_b_label if paired else None,
        # False when analyze could not find preprocess provenance for these VCFs, in which case
        # clustering_config and gatk_version are unknown rather than "GATK defaults".
        "clustering_config_recorded": provenance is not None,
        "clustering_config": provenance.get("clustering_config") if provenance is not None else None,
        "gatk_version": provenance.get("gatk_version") if provenance is not None else None,
        "algorithm_dedup": bool(config.dedup_algorithms),
        "pass_only": bool(config.pass_only),
        "context_overlap": float(config.context_overlap),
    }


def _table_records(frame: pd.DataFrame) -> List[dict]:
    if frame.empty:
        return []
    return json.loads(frame.to_json(orient="records"))


def write_summary_json(
    config: AnalysisConfig,
    summary_tables: Mapping[str, pd.DataFrame],
    modules_ran: List[str],
    modules_skipped: List[str],
) -> Path:
    """Write summary.json; table columns use fixed _a/_b suffixes, never run labels."""
    provenance = load_preprocess_provenance(config.vcf_a_path, config.vcf_b_path if config.mode == AnalysisMode.PAIRED else None)
    payload = {
        "labels": build_labels_block(config, provenance),
        "modules_ran": list(modules_ran),
        "modules_skipped": list(modules_skipped),
        "tables": {key: _table_records(frame) for key, frame in sorted(summary_tables.items())},
    }
    config.output_dir.mkdir(parents=True, exist_ok=True)
    path = config.output_dir / SUMMARY_FILE_NAME
    path.write_text(json.dumps(payload, indent=2) + "\n")
    return path
