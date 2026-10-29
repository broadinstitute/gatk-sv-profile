from __future__ import annotations

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

from gatk_sv_profile.modules.site_overlap import SiteOverlapModule, _plot_heatmap, build_overlap_metrics


def test_build_overlap_metrics_counts_matches_for_both_callsets(module_test_context) -> None:
    metrics = build_overlap_metrics(module_test_context.data.sites_a, module_test_context.data.sites_b, pass_only=True)

    del_metrics = metrics.loc[(metrics["svtype"] == "DEL") & (metrics["algorithm"] == "manta")]
    assert int(del_metrics["n_total_a"].sum()) == 1
    assert int(del_metrics["n_matched_a"].sum()) == 1
    assert int(del_metrics["n_total_b"].sum()) == 1
    assert int(del_metrics["n_matched_b"].sum()) == 1


def test_site_overlap_module_writes_tables_and_plots(module_test_context) -> None:
    module = SiteOverlapModule()

    module.run(module_test_context.data, module_test_context.config)

    output_dir = module_test_context.config.output_dir / "site_overlap"
    table = pd.read_csv(output_dir / "tables" / "overlap_metrics.tsv.gz", sep="\t")
    assert not table.empty
    assert (output_dir / "tables" / "overlap_metrics.parquet").exists()
    assert (output_dir / "overlap.by_class.CallsetA.png").exists()
    assert (output_dir / "overlap.by_algorithm.CallsetA.png").exists()
    assert (output_dir / "overlap.by_evidence.CallsetA.png").exists()
    assert (output_dir / "overlap.by_context.CallsetA.png").exists()
    assert (output_dir / "overlap.by_size.CallsetB.png").exists()
    assert (output_dir / "heatmap.size_x_class.CallsetA.png").exists()
    assert (output_dir / "heatmap.freq_x_class.CallsetB.png").exists()
    assert (output_dir / "heatmap.size_x_algorithm.CallsetA.png").exists()
    assert (output_dir / "heatmap.freq_x_algorithm.CallsetB.png").exists()
    assert (output_dir / "heatmap.size_x_evidence.CallsetA.png").exists()
    assert (output_dir / "heatmap.freq_x_evidence.CallsetB.png").exists()
    assert (output_dir / "heatmap.size_x_context.CallsetA.png").exists()
    assert (output_dir / "heatmap.freq_x_context.CallsetA.png").exists()


def test_plot_heatmap_includes_bnd_and_excludes_ctx_columns(tmp_path, monkeypatch) -> None:
    sites = pd.DataFrame(
        [
            {"variant_id": "a", "size_bucket": "<100bp", "af_bucket": "AC=1", "svtype": "DEL", "truth_vid": "b", "status": "MATCHED"},
            {"variant_id": "b", "size_bucket": "<100bp", "af_bucket": "AC=1", "svtype": "BND", "truth_vid": "c", "status": "MATCHED"},
            {"variant_id": "c", "size_bucket": "<100bp", "af_bucket": "AC=1", "svtype": "CTX", "truth_vid": "d", "status": "MATCHED"},
        ]
    )
    captured = {}
    original_savefig = Figure.savefig

    def spy_savefig(self, *args, **kwargs):
        captured["xlabels"] = [tick.get_text() for tick in self.axes[0].get_xticklabels()]
        return original_savefig(self, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", spy_savefig)

    _plot_heatmap(sites, "size_bucket", "svtype", tmp_path / "heatmap.png", "CallsetA")

    assert (tmp_path / "heatmap.png").exists()
    assert "BND" in captured["xlabels"]
    assert "CTX" not in captured["xlabels"]
    plt.close("all")


def test_plot_heatmap_labels_colorbar(tmp_path, monkeypatch) -> None:
    sites = pd.DataFrame(
        [
            {"variant_id": "a", "size_bucket": "<100bp", "af_bucket": "AC=1", "svtype": "DEL", "truth_vid": "b", "status": "MATCHED"},
            {"variant_id": "b", "size_bucket": "<100bp", "af_bucket": "AC=1", "svtype": "DEL", "truth_vid": None, "status": "UNMATCHED"},
        ]
    )
    captured = {}
    original_colorbar = Figure.colorbar

    def spy_colorbar(self, *args, **kwargs):
        colorbar = original_colorbar(self, *args, **kwargs)
        captured["colorbar"] = colorbar
        return colorbar

    monkeypatch.setattr(Figure, "colorbar", spy_colorbar)

    _plot_heatmap(sites, "size_bucket", "svtype", tmp_path / "heatmap_labeled.png", "CallsetA")

    assert captured["colorbar"].ax.get_ylabel() == "Matched fraction"
    plt.close("all")


def test_plot_heatmap_includes_count_labels(tmp_path, monkeypatch) -> None:
    sites = pd.DataFrame(
        [
            {"variant_id": "a", "size_bucket": "<100bp", "af_bucket": "AC=1", "svtype": "DEL", "truth_vid": "b", "status": "MATCHED"},
            {"variant_id": "b", "size_bucket": "<100bp", "af_bucket": "AC=1", "svtype": "DEL", "truth_vid": None, "status": "UNMATCHED"},
        ]
    )
    captured = {}
    original_savefig = Figure.savefig

    def spy_savefig(self, *args, **kwargs):
        captured["texts"] = [text.get_text() for axis in self.axes for text in axis.texts]
        return original_savefig(self, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", spy_savefig)

    _plot_heatmap(sites, "size_bucket", "svtype", tmp_path / "heatmap_counts.png", "CallsetA")

    assert "n=2" in captured["texts"]
    plt.close("all")


def test_build_site_match_table_both_directions(module_test_context) -> None:
    from gatk_sv_profile.modules.site_overlap import SITE_MATCH_COLUMNS, build_site_match_table

    data = module_test_context.data
    a_to_b = build_site_match_table(data.sites_a, data.sites_b, "a")
    b_to_a = build_site_match_table(data.sites_b, data.sites_a, "b")

    assert list(a_to_b.columns) == SITE_MATCH_COLUMNS
    assert list(b_to_a.columns) == SITE_MATCH_COLUMNS
    assert len(a_to_b) == len(data.sites_a)
    assert len(b_to_a) == len(data.sites_b)

    matched = a_to_b.set_index("vid_a").loc["a_del"]
    assert matched["vid_b"] == "b_del"
    assert matched["filter_a"] == "PASS"
    assert matched["filter_b"] == "PASS"
    assert matched["algorithms_a"] == "manta,wham"
    assert matched["algorithms_b"] == "manta"
    assert pd.isna(a_to_b.set_index("vid_a").loc["a_dup", "vid_b"])
    assert a_to_b.set_index("vid_a").loc["a_dup", "filter_a"] == "UNRESOLVED"
    assert pd.isna(b_to_a.set_index("vid_b").loc["b_inv", "vid_a"])


def test_site_overlap_module_writes_site_match_tables_only_when_enabled(module_test_context) -> None:
    tables_dir = module_test_context.config.output_dir / "site_overlap" / "tables"
    SiteOverlapModule().run(module_test_context.data, module_test_context.config)
    assert not (tables_dir / "site_match.a_to_b.tsv.gz").exists()

    module_test_context.config.enable_site_match_table = True
    SiteOverlapModule().run(module_test_context.data, module_test_context.config)
    a_to_b = pd.read_csv(tables_dir / "site_match.a_to_b.tsv.gz", sep="\t")
    b_to_a = pd.read_csv(tables_dir / "site_match.b_to_a.tsv.gz", sep="\t")
    assert "vid_a" in a_to_b.columns and "filter_b" in b_to_a.columns
    assert not any("CallsetA" in column for column in list(a_to_b.columns) + list(b_to_a.columns))


def test_site_overlap_records_label_free_summary_table(module_test_context) -> None:
    SiteOverlapModule().run(module_test_context.data, module_test_context.config)

    summary = module_test_context.data.summary_tables["site_overlap/overlap_metrics"]
    assert {"n_total_a", "n_matched_b"}.issubset(summary.columns)
    labeled = pd.read_csv(module_test_context.config.output_dir / "site_overlap" / "tables" / "overlap_metrics.tsv.gz", sep="\t")
    assert "n_total_CallsetA" in labeled.columns
