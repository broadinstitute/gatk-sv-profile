# gatk-sv-profile

`gatk-sv-profile` is a Python command-line package for comparing two GATK-SV VCF callsets.
It provides a single workflow for validating inputs, optionally preprocessing the two VCFs
with GATK structural-variant concordance annotations, and generating comparison tables and
plots across site-level and genotype-level quality-control dimensions.

The package is intended for side-by-side evaluation of two callsets from the same cohort or
from closely related cohorts where overlap, allele frequency, genotype behavior, and family
inheritance patterns are useful comparison signals.

## Workflow

The CLI exposes four subcommands:

```text
validate -> preprocess -> analyze
        \-> run
```

- `validate` checks whether a VCF looks like a compatible GATK-SV-style input.
- `preprocess` runs GATK `SVConcordance`-based annotation on both callsets.
- `analyze` reads two annotated VCFs and writes module-specific plots and tables.
- `run` executes `preprocess` and `analyze` end to end.

## Installation

Python 3.9 or newer is required.

Install the package from the repository root:

```bash
python -m pip install .
```

For development:

```bash
python -m pip install -e '.[dev]'
```

### Primary Dependencies

- `pysam` for VCF I/O and indexing
- `pandas` and `numpy` for data aggregation
- `matplotlib` for plots
- `pyarrow` for parquet outputs
- `scipy` for statistical summaries
- `upsetplot` for set-overlap visualizations

### External Requirement for Preprocessing

The `preprocess` and `run` subcommands shell out to GATK. A compatible `gatk` executable
must be available on `PATH`, or passed explicitly with `--gatk-path`.

## Input Expectations

The tool expects two structural-variant VCFs with GATK-SV-style symbolic alleles and
per-sample genotype fields. In practice:

- `SVTYPE` should be present.
- `GT` should be present.
- `GQ` and `ECN` are strongly preferred.
- Final annotated GATK-SV VCFs are the intended input shape.

The `validate` command reports structural issues before you run the rest of the workflow.

## Usage

### 1. Validate a VCF

```bash
gatk-sv-profile validate --vcf callset_a.vcf.gz
```

If the file only has fixable issues, you can request a rewritten output VCF:

```bash
gatk-sv-profile validate \
  --vcf callset_a.vcf.gz \
  --fix \
  --out callset_a.fixed.vcf.gz
```

Optional repair helpers are available for some cases:

- `--ploidy-table` for repairing missing `ECN`
- `--drop-bnd` to discard all `BND` records during fix mode
- `--drop-ctx` to discard all `CTX` records during fix mode

### 2. Preprocess Two Callsets

`preprocess` runs concordance-style annotation across the requested contigs and produces two
annotated VCFs for downstream comparison.

```bash
gatk-sv-profile preprocess \
  --vcf-a callset_a.vcf.gz \
  --vcf-b callset_b.vcf.gz \
  --reference-dict reference.dict \
  --contig-list contigs.list \
  --output-dir compare_out
```

Useful options:

- `--contig chr1` to preprocess a single contig
- `--num-workers 4` to control contig-level parallelism
- `--seg-dup-track`, `--simple-repeat-track`, `--repeatmasker-track` for extra interval annotations
- `--gatk-path` and `--java-options` to control the GATK invocation

This command prints the resolved output paths as:

```text
annotated_a=...
annotated_b=...
```

### 3. Analyze Annotated VCFs

```bash
gatk-sv-profile analyze \
  --vcf-a compare_out/preprocess/concordance_a.vcf.gz \
  --vcf-b compare_out/preprocess/concordance_b.vcf.gz \
  --label-a CallsetA \
  --label-b CallsetB \
  --output-dir compare_out
```

Common options:

- `--modules site_overlap,overall_counts,allele_freq`
- `--pass-only`
- `--context-overlap 0.5`
- `--per-chrom`
- `--enable-site-match-table`
- `--per-sample-counts-table`
- `--ped cohort.ped`
- `--num-workers 4`

The command reports the output root and executed module names:

```text
analysis_output_dir=...
modules_ran=...
modules_skipped=...
```

### 4. Run the End-to-End Workflow

```bash
gatk-sv-profile run \
  --vcf-a callset_a.vcf.gz \
  --vcf-b callset_b.vcf.gz \
  --label-a CallsetA \
  --label-b CallsetB \
  --reference-dict reference.dict \
  --contig-list contigs.list \
  --output-dir compare_out \
  --ped cohort.ped
```

This is the simplest entry point when you want concordance annotation and downstream
comparison in one command.

## Analysis Modules

By default, `analyze` runs all available modules. You can restrict execution with `--modules`.

Site-level modules:

- `binned_counts`
- `overall_counts`
- `site_overlap`
- `allele_freq`
- `genotype_dist`
- `genotype_quality`
- `counts_per_genome`
- `size_signatures`
- `upset`

Genotype-level modules:

- `genotype_exact_match`
- `genotype_concordance`
- `family_analysis`

`family_analysis` requires a PED file via `--ped`. Some genotype-oriented modules may be
skipped automatically when the inputs do not share samples.

## Output Layout

Outputs are organized into subdirectories under `--output-dir`.

- `preprocess/` contains intermediate and annotated concordance VCFs.
- Each analysis module writes into its own directory, such as `site_overlap/` or
  `genotype_quality/`.
- Many modules write both `tables/` outputs and PNG figures.
- Some modules also write parquet tables for downstream programmatic reuse.

Typical examples include:

- `site_overlap/tables/overlap_metrics.tsv.gz`
- `overall_counts/sv_count_by_type.<label>.png`
- `allele_freq/tables/af_correlation_stats.tsv.gz`
- `genotype_quality/tables/gq_summary.<label>.tsv.gz`
- `family_analysis/tables/inheritance_stats.trios.tsv.gz`

## Development

Run the test suite from the repository root:

```bash
pytest
```

The project includes tests for CLI behavior, preprocessing orchestration, validation logic,
and each individual analysis module.

## Repository Layout

```text
src/gatk_sv_profile/
tests/
pyproject.toml
PLAN.md
```

`PLAN.md` contains the implementation plan and design notes. `README.md` is the user-facing
entry point for installation and usage.