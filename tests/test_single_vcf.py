"""Tests for single-VCF mode (Steps 1-4, 10-12 of the refactor plan)."""

from __future__ import annotations

import pytest

from gatk_sv_profile.aggregate import aggregate
from gatk_sv_profile.cli import build_parser, main
from gatk_sv_profile.config import AnalysisConfig, AnalysisMode
from gatk_sv_profile.modules.allele_freq import AlleleFreqModule
from gatk_sv_profile.modules.genotype_concordance import GenotypeConcordanceModule
from gatk_sv_profile.modules.genotype_exact_match import GenotypeExactMatchModule
from gatk_sv_profile.modules.overall_counts import OverallCountsModule
from gatk_sv_profile.modules.site_overlap import SiteOverlapModule


# ---------------------------------------------------------------------------
# Step 1 — AnalysisMode enum
# ---------------------------------------------------------------------------

class TestAnalysisMode:
    def test_default_is_paired(self):
        config = AnalysisConfig()
        assert config.mode is AnalysisMode.PAIRED

    def test_single_mode_roundtrip(self):
        config = AnalysisConfig(mode=AnalysisMode.SINGLE)
        assert config.mode is AnalysisMode.SINGLE


# ---------------------------------------------------------------------------
# Step 2 — AggregatedData.is_paired
# ---------------------------------------------------------------------------

class TestAggregatedDataIsPaired:
    def test_is_paired_true_when_sites_b_present(self, tmp_path, make_vcf):
        extra_headers = [
            "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Calling algorithms\">",
            "##INFO=<ID=EVIDENCE,Number=.,Type=String,Description=\"Evidence types\">",
            "##INFO=<ID=STATUS,Number=1,Type=String,Description=\"Match status\">",
            "##INFO=<ID=TRUTH_VID,Number=1,Type=String,Description=\"Truth variant id\">",
        ]
        vcf_a = make_vcf(file_name="a.vcf", extra_header_lines=extra_headers,
                         records=["chr1\t100\ta1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=100;ALGORITHMS=manta;EVIDENCE=RD\tGT:GQ:ECN\t0/1:30:2\t0/0:20:2"])
        vcf_b = make_vcf(file_name="b.vcf", extra_header_lines=extra_headers,
                         records=["chr1\t110\tb1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=110;ALGORITHMS=manta;EVIDENCE=RD\tGT:GQ:ECN\t0/0:30:2\t0/1:20:2"])
        config = AnalysisConfig(
            vcf_a_path=vcf_a, vcf_b_path=vcf_b, contigs=["chr1"],
            output_dir=tmp_path / "out", n_workers=1,
        )
        data = aggregate(config)
        assert data.is_paired is True
        assert data.sites_b is not None
        assert data.matched_pairs is not None
        assert data.sample_names_b is not None

    def test_is_paired_false_in_single_mode(self, tmp_path, make_vcf):
        extra_headers = [
            "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Calling algorithms\">",
            "##INFO=<ID=EVIDENCE,Number=.,Type=String,Description=\"Evidence types\">",
        ]
        vcf = make_vcf(file_name="a.vcf", extra_header_lines=extra_headers,
                       records=["chr1\t100\ta1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=100;ALGORITHMS=manta;EVIDENCE=RD\tGT:GQ:ECN\t0/1:30:2\t0/0:20:2"])
        config = AnalysisConfig(
            mode=AnalysisMode.SINGLE,
            vcf_a_path=vcf, contigs=["chr1"],
            output_dir=tmp_path / "out", n_workers=1,
        )
        data = aggregate(config)
        assert data.is_paired is False
        assert data.sites_b is None
        assert data.matched_pairs is None
        assert data.sample_names_b is None
        assert data.label_b is None
        assert data.shared_samples == []


# ---------------------------------------------------------------------------
# Step 3 — aggregate() single-VCF path produces correct shape
# ---------------------------------------------------------------------------

class TestAggregateSingle:
    @pytest.fixture
    def single_data(self, tmp_path, make_vcf):
        extra_headers = [
            "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Calling algorithms\">",
            "##INFO=<ID=EVIDENCE,Number=.,Type=String,Description=\"Evidence types\">",
            "##INFO=<ID=STATUS,Number=1,Type=String,Description=\"Match status\">",
            "##INFO=<ID=TRUTH_VID,Number=1,Type=String,Description=\"Truth variant id\">",
        ]
        vcf = make_vcf(
            file_name="single.vcf",
            sample_names=["S1", "S2", "S3"],
            extra_header_lines=extra_headers,
            records=[
                "chr1\t100\tv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=200;ALGORITHMS=manta;EVIDENCE=RD,PE\tGT:GQ:ECN\t0/1:50:2\t0/0:60:2\t1/1:70:2",
                "chr1\t400\tv2\tN\t<INS>\t.\tPASS\tSVTYPE=INS;SVLEN=80;ALGORITHMS=wham;EVIDENCE=SR\tGT:GQ:ECN\t0/0:30:2\t0/1:40:2\t0/0:50:2",
            ],
        )
        config = AnalysisConfig(
            mode=AnalysisMode.SINGLE,
            vcf_a_path=vcf, vcf_a_label="MyVCF",
            contigs=["chr1"], output_dir=tmp_path / "out", n_workers=1,
        )
        return aggregate(config)

    def test_sites_a_shape(self, single_data):
        assert single_data.sites_a.shape[0] == 2

    def test_label_a(self, single_data):
        assert single_data.label_a == "MyVCF"

    def test_sample_names_a(self, single_data):
        assert single_data.sample_names_a == ["S1", "S2", "S3"]

    def test_parquet_written(self, single_data):
        assert (single_data.site_table_dir / "sites_a.all.parquet").exists()

    def test_no_b_side_files(self, single_data):
        b_files = list(single_data.site_table_dir.glob("sites_b*"))
        assert b_files == []

    def test_raises_without_vcf_a(self, tmp_path):
        config = AnalysisConfig(mode=AnalysisMode.SINGLE, contigs=["chr1"], output_dir=tmp_path)
        from gatk_sv_profile.aggregate import aggregate as agg
        with pytest.raises(ValueError, match="vcf_a_path"):
            agg(config)

    def test_paired_mode_raises_without_vcf_b(self, tmp_path, make_vcf):
        vcf = make_vcf(file_name="a.vcf", records=["chr1\t100\tv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=50\tGT:GQ:ECN\t0/1:30:2\t0/0:20:2"])
        config = AnalysisConfig(
            mode=AnalysisMode.PAIRED,
            vcf_a_path=vcf, contigs=["chr1"], output_dir=tmp_path,
        )
        from gatk_sv_profile.aggregate import aggregate as agg
        with pytest.raises(ValueError, match="vcf_b_path"):
            agg(config)


# ---------------------------------------------------------------------------
# Step 4 — requires_paired_input gating
# ---------------------------------------------------------------------------

class TestRequiresPairedInput:
    def test_paired_only_modules_declare_flag(self):
        assert SiteOverlapModule().requires_paired_input is True
        assert AlleleFreqModule().requires_paired_input is True
        assert GenotypeConcordanceModule().requires_paired_input is True
        assert GenotypeExactMatchModule().requires_paired_input is True

    def test_per_vcf_module_does_not_require_paired(self):
        assert OverallCountsModule().requires_paired_input is False

    def test_paired_only_modules_skipped_in_single_mode(self, tmp_path, make_vcf, capsys, monkeypatch):
        extra_headers = [
            "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Calling algorithms\">",
            "##INFO=<ID=EVIDENCE,Number=.,Type=String,Description=\"Evidence types\">",
        ]
        vcf = make_vcf(
            file_name="single_skip.vcf",
            extra_header_lines=extra_headers,
            records=["chr1\t100\tv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=100;ALGORITHMS=manta;EVIDENCE=RD\tGT:GQ:ECN\t0/1:30:2\t0/0:20:2"],
        )
        monkeypatch.setattr(
            "gatk_sv_profile.cli.ALL_MODULES",
            [SiteOverlapModule, AlleleFreqModule, GenotypeConcordanceModule, GenotypeExactMatchModule],
        )
        exit_code = main([
            "analyze", "--vcf", str(vcf),
            "--output-dir", str(tmp_path / "out"),
            "--modules", "site_overlap,allele_freq,genotype_concordance,genotype_exact_match",
        ])
        out = capsys.readouterr().out
        assert exit_code == 0
        assert "modules_ran=" in out
        assert "modules_skipped=site_overlap,allele_freq,genotype_concordance,genotype_exact_match" in out


# ---------------------------------------------------------------------------
# Step 10 — CLI argument parsing for single-VCF mode
# ---------------------------------------------------------------------------

class TestCliSingleVcfArgParsing:
    def test_analyze_accepts_single_vcf_flag(self, tmp_path):
        parser = build_parser()
        args = parser.parse_args([
            "analyze", "--vcf", "/some/file.vcf.gz",
            "--output-dir", str(tmp_path),
        ])
        assert args.vcf is not None
        assert args.vcf_a is None

    def test_analyze_accepts_paired_flags(self, tmp_path):
        parser = build_parser()
        args = parser.parse_args([
            "analyze",
            "--vcf-a", "/some/a.vcf.gz",
            "--vcf-b", "/some/b.vcf.gz",
            "--output-dir", str(tmp_path),
        ])
        assert args.vcf is None
        assert args.vcf_a is not None

    def test_analyze_vcf_and_vcf_a_are_mutually_exclusive(self, tmp_path):
        parser = build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args([
                "analyze",
                "--vcf", "/some/file.vcf.gz",
                "--vcf-a", "/some/a.vcf.gz",
                "--output-dir", str(tmp_path),
            ])

    def test_analyze_requires_at_least_one_vcf_flag(self, tmp_path):
        parser = build_parser()
        with pytest.raises(SystemExit):
            parser.parse_args(["analyze", "--output-dir", str(tmp_path)])

    def test_analyze_single_mode_sets_label_default(self, tmp_path):
        parser = build_parser()
        args = parser.parse_args([
            "analyze", "--vcf", "/some/file.vcf.gz",
            "--output-dir", str(tmp_path),
        ])
        assert args.label == "VCF"

    def test_analyze_single_mode_custom_label(self, tmp_path):
        parser = build_parser()
        args = parser.parse_args([
            "analyze", "--vcf", "/some/file.vcf.gz",
            "--label", "MyCallset",
            "--output-dir", str(tmp_path),
        ])
        assert args.label == "MyCallset"

    def test_preprocess_single_vcf_flag(self, tmp_path):
        parser = build_parser()
        args = parser.parse_args([
            "preprocess", "--vcf", "/some/file.vcf.gz",
            "--reference-dict", str(tmp_path / "ref.dict"),
            "--contig-list", str(tmp_path / "contigs.list"),
            "--output-dir", str(tmp_path),
        ])
        assert args.vcf is not None
        assert args.vcf_a is None

    def test_run_single_vcf_flag(self, tmp_path):
        parser = build_parser()
        args = parser.parse_args([
            "run", "--vcf", "/some/file.vcf.gz",
            "--reference-dict", str(tmp_path / "ref.dict"),
            "--contig-list", str(tmp_path / "contigs.list"),
            "--output-dir", str(tmp_path),
        ])
        assert args.vcf is not None


# ---------------------------------------------------------------------------
# Step 10 — CLI handler builds SINGLE-mode AnalysisConfig
# ---------------------------------------------------------------------------

class TestCliAnalyzeSingleMode:
    def test_analyze_single_builds_correct_config(self, make_vcf, tmp_path, monkeypatch, capsys):
        extra_headers = [
            "##INFO=<ID=ALGORITHMS,Number=.,Type=String,Description=\"Calling algorithms\">",
            "##INFO=<ID=EVIDENCE,Number=.,Type=String,Description=\"Evidence types\">",
        ]
        vcf = make_vcf(
            file_name="s.vcf",
            extra_header_lines=extra_headers,
            records=["chr1\t100\tv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=100;ALGORITHMS=manta;EVIDENCE=RD\tGT:GQ:ECN\t0/1:30:2\t0/0:20:2"],
        )
        captured = {}

        def fake_aggregate(config):
            captured["config"] = config
            import pandas as pd
            from types import SimpleNamespace
            return SimpleNamespace(
                sites_a=pd.DataFrame([{"variant_id": "v1"}]),
                sites_b=None,
                matched_pairs=None,
                shared_samples=[],
                label_a=config.vcf_a_label,
                label_b=None,
            )

        from gatk_sv_profile.modules.overall_counts import OverallCountsModule as OCM
        monkeypatch.setattr("gatk_sv_profile.cli.ALL_MODULES", [])
        monkeypatch.setattr("gatk_sv_profile.cli.aggregate", fake_aggregate)

        exit_code = main([
            "analyze", "--vcf", str(vcf),
            "--label", "Pilot",
            "--output-dir", str(tmp_path / "out"),
        ])
        assert exit_code == 0
        config = captured["config"]
        assert config.mode is AnalysisMode.SINGLE
        assert config.vcf_a_path == vcf
        assert config.vcf_b_path is None
        assert config.vcf_a_label == "Pilot"
        assert "chr1" in config.contigs

    def test_analyze_paired_mode_missing_vcf_b_returns_error(self, make_vcf, tmp_path, capsys):
        vcf = make_vcf(
            file_name="a.vcf",
            records=["chr1\t100\tv1\tN\t<DEL>\t.\tPASS\tSVTYPE=DEL;SVLEN=50\tGT:GQ:ECN\t0/1:30:2\t0/0:20:2"],
        )
        exit_code = main([
            "analyze", "--vcf-a", str(vcf),
            "--output-dir", str(tmp_path / "out"),
        ])
        assert exit_code == 2
