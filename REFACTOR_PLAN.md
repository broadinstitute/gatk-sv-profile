# `gatk-sv-profile` — One-VCF / Two-VCF Refactor Plan

## Goal

Generalize the package so the same CLI commands (`validate`, `preprocess`, `analyze`, `run`)
can operate on either a single VCF (profiling mode) or a pair of VCFs (comparison mode,
current behavior). The mode is selected purely by which input flags the user provides.

## Design Principles

1. **Single CLI surface.** No new subcommands. `--vcf-a/--vcf-b` (pair) and `--vcf` (single)
   are mutually exclusive; the chosen pair determines the mode.
2. **Mode is a first-class enum**, threaded through `AnalysisConfig`. Modules read it; they
   do not infer mode from `data.sites_b is None`.
3. **Panel reuse over module forking.** Per-callset visual content (one figure for one VCF)
   is implemented once as a `Panel` class. Two reports — `SingleCallsetReport` and
   `PairedCallsetReport` — drive the same panels with different inputs.
4. **Comparison modules are gated, not duplicated.** Modules that fundamentally require two
   callsets (overlap, AF correlation, genotype concordance, genotype exact match) declare
   `requires_paired_input = True` and are skipped in single-VCF mode with a logged reason.
5. **No behavior change for paired runs.** All existing outputs, file names, and table
   schemas must remain byte-identical for the two-VCF path. Add tests asserting this.
6. **Small, mergeable steps.** Each numbered step below should land as one PR/commit with
   tests green. Steps are ordered; later steps depend on earlier ones.

---

## Module Classification

This drives gating and report wiring. Lock this table before coding.

| Module | Single-VCF | Paired | Notes |
|---|---|---|---|
| `overall_counts` | yes | yes (runs per VCF) | Already iterates per VCF |
| `genotype_dist` (HWE, carrier freq) | yes | yes (runs per VCF) | Per-VCF only |
| `genotype_quality` | yes | yes (runs per VCF) | Per-VCF only |
| `counts_per_genome` | yes | yes (runs per VCF) | Per-VCF only |
| `upset` | yes | yes (runs per VCF) | Per-VCF only |
| `size_signatures` | yes | yes (runs per VCF) | Per-VCF only |
| `family_analysis` | yes | yes (runs per VCF) | Per-VCF; pedigree intersected with each VCF independently |
| `binned_counts` | **single-VCF variant**: emit `counts.tsv` for one VCF, no `_a/_b` combo | yes | Refactor combiner to wrap a per-VCF summarizer |
| `site_overlap` | **skipped** | yes | Requires `truth_vid`/`STATUS` |
| `allele_freq` | **skipped** | yes | Requires matched pairs |
| `genotype_concordance` | **skipped** | yes | Requires SVConcordance across two callsets |
| `genotype_exact_match` | **skipped** | yes | Requires shared samples across two callsets |

---

## Step 1 — Introduce `AnalysisMode` and thread it through `AnalysisConfig`

- Add `AnalysisMode` enum (`SINGLE`, `PAIRED`) in [src/gatk_sv_profile/config.py](src/gatk_sv_profile/config.py).
- Add `mode: AnalysisMode = AnalysisMode.PAIRED` to `AnalysisConfig` (default keeps current behavior).
- Do not touch CLI or modules yet. Just expose the field.
- Unit test: dataclass defaults round-trip.

## Step 2 — Make VCF-B optional on the data model

- In `AggregatedData` ([src/gatk_sv_profile/aggregate.py](src/gatk_sv_profile/aggregate.py)),
  make every B-side field `Optional` (`sites_b`, `sample_names_b`, `sample_indices_b`,
  `label_b`) and `matched_pairs` optional. Keep paired callers passing real values.
- Add helper `AggregatedData.is_paired -> bool` returning `sites_b is not None`.
- Update type hints; no functional change yet. All existing tests must still pass.

## Step 3 — Split `aggregate()` into per-VCF + pairing stages

- Extract a private `_aggregate_single(vcf_path, label, config) -> SingleAggregated` returning
  `(sites_df, sample_names, site_table_dir)`.
- Rewrite `aggregate(config)` to call `_aggregate_single` once or twice based on
  `config.mode` (or based on whether `vcf_b_path` is set, but prefer the explicit mode).
- In `SINGLE` mode return `AggregatedData` with `sites_b=None`, `matched_pairs=None`,
  `sample_names_b=None`, `shared_samples=[]`, etc.
- Add `_aggregate_single` unit test against a small fixture VCF.

## Step 4 — Add `requires_paired_input` to `AnalysisModule`

- In [src/gatk_sv_profile/modules/base.py](src/gatk_sv_profile/modules/base.py) add property
  `requires_paired_input: bool = False`.
- Set `True` on: `SiteOverlapModule`, `AlleleFreqModule`, `GenotypeConcordanceModule`,
  `GenotypeExactMatchModule`.
- Extend `_skip_reason` in [src/gatk_sv_profile/cli.py](src/gatk_sv_profile/cli.py) to skip
  with reason `"requires two VCFs"` when `config.mode == SINGLE` and module requires paired.
- Existing paired-mode behavior unchanged. Tests: assert all four modules report
  `requires_paired_input is True`; assert skip reason wiring in single mode.

---

## Step 5 — Plotting class hierarchy: introduce `Panel` and `Report`

Goal: a single per-VCF panel implementation runs once or twice without if/else inside
modules.

- New file `src/gatk_sv_profile/plotting/__init__.py` exposing `Panel`, `SingleCallsetReport`,
  `PairedCallsetReport`.
- New file `src/gatk_sv_profile/plotting/panels.py` with:
  ```python
  class Panel(ABC):
      """One reusable figure unit operating on a single callset."""
      name: str  # filename stem, e.g. "size_distribution"
      @abstractmethod
      def render(self, sites: pd.DataFrame, *, label: str, output_path: Path, config: AnalysisConfig) -> None: ...
  ```
- New file `src/gatk_sv_profile/plotting/reports.py` with:
  ```python
  class CallsetReport(ABC):
      panels: list[Panel]
      comparison_panels: list[ComparisonPanel]  # only used by paired
      @abstractmethod
      def run(self, data: AggregatedData, config: AnalysisConfig, output_dir: Path) -> None: ...

  class SingleCallsetReport(CallsetReport):
      """Calls each Panel exactly once with sites_a/label_a."""

  class PairedCallsetReport(CallsetReport):
      """Calls each Panel twice (A, B) and each ComparisonPanel once."""
  ```
- `ComparisonPanel` is a separate ABC whose `render` takes both callset frames
  (e.g. matched-pair scatter, AF correlation).
- Add no panels yet. Just the scaffold + interface tests.

## Step 6 — Port `overall_counts` panels (pilot module)

This is the template all subsequent module ports follow. Do it carefully.

- Convert each private `_plot_*` helper in
  [src/gatk_sv_profile/modules/overall_counts.py](src/gatk_sv_profile/modules/overall_counts.py)
  into a `Panel` subclass living in `src/gatk_sv_profile/plotting/panels/overall_counts.py`:
  `TypeCountsPanel`, `SizeDistributionPanel`, `AfDistributionPanel`, `ContextByTypePanel`.
- Each `Panel.__init__` takes its configuration (e.g. `field`, `group_field`, filename
  template, title template).
- Rewrite `OverallCountsModule.run` to instantiate a `SingleCallsetReport` or
  `PairedCallsetReport` (selected from `config.mode`) wired with the panel list and call
  `report.run(data, config, output_dir)`.
- The `PairedCallsetReport` must produce the same file names as today
  (`*.{label}.png` per VCF). Snapshot existing output filenames before and after with a
  test that runs the module on a tiny fixture and asserts the set of emitted paths.

## Step 7 — Port remaining per-VCF modules to `Panel`s

Repeat the Step 6 procedure for each. Land one PR per module:

1. `genotype_dist` → panels under `plotting/panels/genotype_dist.py`.
2. `genotype_quality` → panels under `plotting/panels/genotype_quality.py`.
3. `counts_per_genome` → panels under `plotting/panels/counts_per_genome.py`.
4. `upset` → panels under `plotting/panels/upset.py`.
5. `size_signatures` → panels (or table-only "Report") under `plotting/panels/size_signatures.py`.
6. `family_analysis` → panels under `plotting/panels/family_analysis.py`.

Acceptance per port: existing module tests in `tests/test_modules/` pass unchanged, and the
set of emitted file paths and table schemas in paired mode is unchanged.

## Step 8 — Port `binned_counts` to support both modes

- Keep `summarize_binned_counts(sites, pass_only)` as the per-VCF primitive (already exists).
- In `BinnedCountsModule.run`:
  - SINGLE: write `counts.tsv(.gz)` and `counts.parquet` from `summarize_binned_counts(data.sites_a, ...)`. Columns: bare `n_variants`, `n_matched`, `n_unmatched` (no `_a/_b` suffix).
  - PAIRED: existing `build_combined_binned_counts` path. No change.
- Add tests for both branches.

## Step 9 — Port paired-only modules as `ComparisonPanel`s

These still execute, but only in paired mode. Convert their plotting helpers into
`ComparisonPanel`s so they live in the same plotting taxonomy.

1. `site_overlap`: per-VCF bar/heatmap helpers become regular `Panel`s (they consume
   `sites_a` and `sites_b` independently, but rely on `truth_vid`/STATUS — they remain
   paired-only). Keep them grouped under `PairedCallsetReport` and require
   `requires_paired_input = True` on the module.
2. `allele_freq`: scatter plots over matched pairs → `ComparisonPanel`.
3. `genotype_concordance`: per-VCF and joint panels.
4. `genotype_exact_match`: per-pair panels.

No CLI behavior change yet.

---

## Step 10 — CLI: add `--vcf` single-VCF flags and validation

Edit [src/gatk_sv_profile/cli.py](src/gatk_sv_profile/cli.py).

- For each subcommand that currently takes `--vcf-a/--vcf-b`:
  - Add `--vcf` (single input). Make `--vcf-a` and `--vcf-b` no longer `required=True`.
  - Add `--label` for the single-VCF label (default `VCF`).
  - Add a post-parse validator that enforces exactly one of:
    - `--vcf` (SINGLE mode), OR
    - both `--vcf-a` and `--vcf-b` (PAIRED mode).
    Reject any other combination with a clear error.
- Affected subcommands: `preprocess`, `analyze`, `run`. `validate` already takes a single
  `--vcf` and is unaffected.
- For `preprocess` in SINGLE mode: run SVRegionOverlap only (skip SVConcordance, since it
  needs a truth callset). Document this clearly in `--help`.
- For `analyze` in SINGLE mode: set `config.mode = SINGLE`, populate `vcf_a_path` from
  `--vcf`, leave `vcf_b_path = None`.
- Update `_infer_common_contigs` callers: in SINGLE mode call a new
  `_read_single_vcf_contigs(vcf_path)` helper instead.
- Update `_build_analysis_config` signature to accept either single or paired inputs and
  set `mode` accordingly. Prefer two separate builder functions over a polymorphic one.

## Step 11 — `preprocess` single-VCF path

Edit [src/gatk_sv_profile/preprocess.py](src/gatk_sv_profile/preprocess.py).

- Refactor `run_preprocess(config)` so the SVConcordance step is conditional on
  `config.mode == PAIRED` (or equivalently on `vcf_b_path is not None`).
- SINGLE-mode `run_preprocess` runs only SVRegionOverlap on `vcf_a_path` and returns
  `(annotated_a, None)`.
- Update CLI handler `_handle_preprocess` / `_handle_run` to print only `annotated_a=...`
  in SINGLE mode.

## Step 12 — `_run_analysis` single-mode wiring

- Update `_run_analysis` logging to handle a `None` `vcf_b_label`: log "%s across %s
  contig(s)..." instead of "%s vs %s...".
- The skip logic added in Step 4 already prevents paired-only modules from running.
- `aggregate(config)` (updated in Step 3) returns a single-VCF `AggregatedData`.
- Ensure the post-run print statements (`analysis_output_dir=`, `modules_ran=`) are mode-agnostic.

---

## Step 13 — Tests

- New `tests/test_cli_single_vcf.py`: argument-parsing validation (mutex of flag groups,
  error messages, label defaults).
- New `tests/test_aggregate_single.py`: single-VCF aggregation returns the expected
  `AggregatedData` shape.
- New `tests/test_modules/test_<module>_single.py` for every per-VCF module that runs in
  single mode: assert it emits the same per-VCF outputs as the A side of a paired run on
  the same fixture.
- New `tests/test_modules/test_paired_only_skipped.py`: in SINGLE mode, the four
  paired-only modules are skipped with the correct reason.
- New `tests/test_plotting/test_reports.py`: `SingleCallsetReport` runs each panel once;
  `PairedCallsetReport` runs each panel twice and each comparison panel once. Use a fake
  `Panel` that records its calls.
- Regression: add a paired-mode "golden file list" test that runs the full analyze
  pipeline on the existing fixture and asserts the set of emitted output paths matches a
  checked-in manifest. Run this before and after each Step 6/7/8/9 port.

## Step 14 — Documentation

- Update [README.md](README.md):
  - Replace the "two VCFs" framing with "one or two VCFs"; show both invocation styles.
  - Document the SINGLE-mode skip list (paired-only modules).
  - Document `preprocess` SINGLE-mode behavior (SVRegionOverlap only).
- Add a short note at the top of [PLAN.md](PLAN.md) pointing to this refactor and marking
  any superseded sections.
- Add an `Architecture` section in README describing the `Panel` / `Report` split, with a
  one-paragraph example.

## Step 15 — Cleanup

- Remove now-dead per-module plotting helpers that have been fully replaced by `Panel`s.
- Delete any compatibility shims introduced in Step 6 once all modules are ported.
- Run `pytest --cov=gatk_sv_profile --cov-report=term-missing` and confirm coverage has
  not regressed.

---

## Out of Scope

- Changing output file names or table schemas for paired mode.
- New analysis content (new metrics, new plots). The single-VCF mode only exposes the
  per-VCF panels we already render in paired mode.
- Multi-VCF (>2) support. Design leaves room for it (the `Report` interface generalizes),
  but no implementation here.

## Risk Notes for the Implementing Agent

- The `Panel` refactor (Steps 5–9) is the largest source of risk. Land a snapshot test of
  emitted file paths and parquet/tsv schemas **before** Step 6 so every later port is
  validated against it.
- Several modules touch the source VCF directly (`genotype_concordance`,
  `genotype_exact_match`, `counts_per_genome`, `genotype_quality`). When porting, be
  careful that VCF I/O happens once per callset, not once per panel.
- `family_analysis` is per-VCF but consumes a single shared PED file. In single-VCF mode,
  it should still work and intersect the PED with the one callset's samples.
- The plotting style (`plot_utils.py`) and dimension helpers (`dimensions.py`) are pure
  utilities and need no changes.
