# gatk-sv-profile

`gatk-sv-profile` is a Python command-line package for profiling one GATK-SV VCF callset or
comparing two. It validates inputs, optionally preprocesses the VCF(s) with GATK structural-variant
annotations, and generates quality-control tables and plots across site-level and genotype-level
dimensions.

- **Single-callset mode** (`--vcf`): profile one VCF — counts, size signatures, genotype distributions,
  family inheritance, and more.
- **Paired-callset mode** (`--vcf-a` / `--vcf-b`): everything above for each callset, plus
  site-overlap, allele-frequency correlation, and genotype-concordance comparisons between the two.

## Workflow

The CLI exposes four subcommands:

```text
validate -> preprocess -> analyze
        \-> run
```

- `validate` checks whether a VCF looks like a compatible GATK-SV-style input.
- `preprocess` runs GATK `SVRegionOverlap` (and `SVConcordance` in paired mode) to annotate the
  VCF(s) before analysis.
- `analyze` reads annotated VCF(s) and writes module-specific plots and tables.
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

The tool expects one or two structural-variant VCFs with GATK-SV-style symbolic alleles and
per-sample genotype fields. In practice:

- `SVTYPE` should be present.
- `GT` should be present.
- `GQ` and `ECN` are strongly preferred.
- Final annotated GATK-SV VCFs are the intended input shape.

The `validate` command reports structural issues before you run the rest of the workflow.

## Usage

### 1. Validate a VCF

```bash
gatk-sv-profile validate --vcf callset.vcf.gz
```

If the file only has fixable issues, you can request a rewritten output VCF:

```bash
gatk-sv-profile validate \
  --vcf callset.vcf.gz \
  --fix \
  --out callset.fixed.vcf.gz
```

Optional repair helpers:

- `--ploidy-table` for repairing missing `ECN`
- `--drop-bnd` to discard all `BND` records during fix mode
- `--drop-ctx` to discard all `CTX` records during fix mode

### 2. Preprocess

`preprocess` annotates VCF(s) with genomic-context and (in paired mode) concordance information.

**Single callset** — runs `SVRegionOverlap` only:

```bash
gatk-sv-profile preprocess \
  --vcf callset.vcf.gz \
  --reference-dict reference.dict \
  --contig-list contigs.list \
  --output-dir profile_out
```

**Two callsets** — runs `SVConcordance` + `SVRegionOverlap`:

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

This command prints the resolved annotated VCF path(s):

```text
annotated_a=...
annotated_b=...   # paired mode only
```

### 3. Analyze

**Single callset:**

```bash
gatk-sv-profile analyze \
  --vcf profile_out/preprocess/annotated_a.vcf.gz \
  --label MyCallset \
  --output-dir profile_out
```

**Two callsets:**

```bash
gatk-sv-profile analyze \
  --vcf-a compare_out/preprocess/annotated_a.vcf.gz \
  --vcf-b compare_out/preprocess/annotated_b.vcf.gz \
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

### 4. Run End to End

**Single callset:**

```bash
gatk-sv-profile run \
  --vcf callset.vcf.gz \
  --label MyCallset \
  --reference-dict reference.dict \
  --contig-list contigs.list \
  --output-dir profile_out
```

**Two callsets:**

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

`run` is the simplest entry point when you want preprocessing and analysis in one command.

## Analysis Modules

By default, `analyze` runs all available modules. Restrict with `--modules`.

Modules that run in both single and paired mode:

| Module | Description |
|---|---|
| `overall_counts` | SV counts by type, size, evidence, and genomic context |
| `genotype_dist` | Genotype and Hardy–Weinberg distributions |
| `genotype_quality` | GQ and per-variant quality summaries |
| `counts_per_genome` | Per-sample SV burden |
| `upset` | Algorithm and evidence-type combination membership |
| `size_signatures` | Size distributions and MEI subtype summaries |
| `binned_counts` | Cross-tabulated counts by size, AF, type, and context |
| `family_analysis` | Transmission and de novo rates (requires `--ped`) |

Modules that require two callsets (automatically skipped in single-callset mode):

| Module | Description |
|---|---|
| `site_overlap` | Per-type site-level recall and precision |
| `allele_freq` | Allele-frequency correlation across matched pairs |
| `genotype_concordance` | Genotype-level concordance across shared samples |
| `genotype_exact_match` | Exact genotype-match rates across shared samples |

`family_analysis` requires a PED file via `--ped`. Genotype-oriented modules may be skipped
automatically when the inputs do not share samples.

## Output Layout

Outputs are organized into subdirectories under `--output-dir`.

- `preprocess/` contains intermediate and annotated VCFs.
- Each analysis module writes into its own directory, e.g. `site_overlap/` or `genotype_quality/`.
- Most modules write both `tables/` outputs (`.tsv.gz`) and PNG figures.
- Some modules also write parquet tables for programmatic reuse.

Examples:

- `overall_counts/sv_count_by_type.<label>.png`
- `site_overlap/tables/overlap_metrics.tsv.gz`
- `allele_freq/tables/af_correlation_stats.tsv.gz`
- `genotype_quality/tables/gq_summary.<label>.tsv.gz`
- `family_analysis/tables/inheritance_stats.trios.tsv.gz`

## Development

Run the test suite from the repository root:

```bash
pytest
```

The project includes tests for CLI behaviour, preprocessing orchestration, validation logic,
and each individual analysis module.

## Repository Layout

```text
src/gatk_sv_profile/
tests/
pyproject.toml
```
