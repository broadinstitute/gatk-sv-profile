from __future__ import annotations

import json
from pathlib import Path

from gatk_sv_profile.config import AnalysisConfig, AnalysisMode
from gatk_sv_profile.provenance import (
    build_labels_block,
    load_preprocess_provenance,
    parse_gatk_version,
    write_preprocess_provenance,
    write_summary_json,
)


def test_parse_gatk_version_extracts_version_line() -> None:
    raw = "Using GATK jar /gatk.jar\nThe Genome Analysis Toolkit (GATK) v4.6.1.0\nHTSJDK Version: 4.1.3\n"
    assert parse_gatk_version(raw) == "4.6.1.0"
    assert parse_gatk_version("4.5.0.0") == "4.5.0.0"
    assert parse_gatk_version("") is None


def test_provenance_round_trip_records_clustering_config(tmp_path) -> None:
    clustering_config = tmp_path / "clustering.yaml"
    clustering_config.write_text("x: 1\n")
    config = AnalysisConfig(mode=AnalysisMode.PAIRED, clustering_config=clustering_config)
    write_preprocess_provenance(
        tmp_path,
        config=config,
        gatk_version="4.6.1.0",
        annotated_a=tmp_path / "annotated_a.vcf.gz",
        annotated_b=tmp_path / "annotated_b.vcf.gz",
    )

    provenance = load_preprocess_provenance(tmp_path / "annotated_a.vcf.gz", tmp_path / "annotated_b.vcf.gz")

    assert provenance is not None
    assert provenance["clustering_config"]["path"] == str(clustering_config.resolve())
    assert len(provenance["clustering_config"]["sha256"]) == 64
    assert provenance["gatk_version"] == "4.6.1.0"


def test_provenance_ignored_for_other_vcfs(tmp_path) -> None:
    config = AnalysisConfig(mode=AnalysisMode.PAIRED)
    write_preprocess_provenance(
        tmp_path,
        config=config,
        gatk_version=None,
        annotated_a=tmp_path / "annotated_a.vcf.gz",
        annotated_b=tmp_path / "annotated_b.vcf.gz",
    )
    assert load_preprocess_provenance(tmp_path / "other.vcf.gz", tmp_path / "annotated_b.vcf.gz") is None
    assert load_preprocess_provenance(tmp_path / "annotated_a.vcf.gz", tmp_path / "other.vcf.gz") is None


def test_labels_block_distinguishes_defaults_from_unknown() -> None:
    config = AnalysisConfig(vcf_a_label="branch", vcf_b_label="baseline", dedup_algorithms=True)

    unknown = build_labels_block(config, None)
    defaults = build_labels_block(config, {"clustering_config": None, "gatk_version": "4.6.1.0"})

    assert unknown["clustering_config_recorded"] is False
    assert defaults["clustering_config_recorded"] is True
    assert defaults["clustering_config"] is None
    assert defaults["algorithm_dedup"] is True
    assert defaults["label_a"] == "branch" and defaults["label_b"] == "baseline"
    assert defaults["schema_version"] == 1


def test_write_summary_json_uses_fixed_keys(tmp_path) -> None:
    import pandas as pd

    config = AnalysisConfig(output_dir=tmp_path, vcf_a_path=tmp_path / "a.vcf.gz", vcf_b_path=tmp_path / "b.vcf.gz")
    tables = {"site_overlap/overlap_metrics": pd.DataFrame({"svtype": ["DEL"], "n_total_a": [3], "n_total_b": [None]})}

    path = write_summary_json(config, tables, ["site_overlap"], [])

    payload = json.loads(Path(path).read_text())
    assert payload["labels"]["clustering_config_recorded"] is False
    assert payload["tables"]["site_overlap/overlap_metrics"] == [{"svtype": "DEL", "n_total_a": 3.0, "n_total_b": None}]


def test_provenance_that_is_not_an_object_is_ignored(tmp_path) -> None:
    for content in ("null", "[]", "\"text\""):
        (tmp_path / "provenance.json").write_text(content)
        assert load_preprocess_provenance(tmp_path / "annotated_a.vcf.gz") is None
