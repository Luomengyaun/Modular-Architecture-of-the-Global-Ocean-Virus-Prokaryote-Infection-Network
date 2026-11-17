# Modular Architecture of the Global Ocean Virus–Prokaryote Infection Network

Analysis code for OVPIN (Ocean Virus–Prokaryote Infection Network), providing tools for recursive BRIM-based multiscale modularity analysis, nestedness assessment, network perturbation, and circular permutation tests for functional enrichment.

The repository is organized into three components, each with a dedicated README describing its inputs, parameters, and outputs.

## Repository structure

```text
.
├── README.md
├── multilayer_modularity_core/
│   ├── README.md
│   ├── analysis_config.json
│   ├── recursive_bimat_worker.m
│   ├── root_null_worker.m
│   ├── core_rules.py
│   ├── root_null_exports.py
│   ├── NullAdaptiveBrim.m
│   ├── bimat_null_modularity.m
│   ├── bimat_nestedness_fast.m
│   ├── bimat_validate_fast_null.m
│   ├── CODE_INTEGRITY.json
│   └── examples/
│       ├── network.tsv
│       └── prepare_jobs.py
├── perturbation_OVPIN/
│   ├── README.md
│   ├── expanding_host_ranges.py
│   ├── replacing_erroneous_links.py
│   └── inferring_hosts_for_unrecalled_viruses.py
└── circular_permutation_functional_enrichment/
    ├── README.md
    ├── requirements.txt
    └── circular_permutation_functional_enrichment.py
```

## Analysis components

### 1. Multiscale modularity and nestedness

[`multilayer_modularity_core`](multilayer_modularity_core/README.md) combines MATLAB workers using BiMat with Python utilities for validating results and selecting child modules for recursive analysis.

- BRIM-based bipartite modularity optimization and hierarchical module selection.
- Null-model evaluation using EQUIPROBABLE, AVERAGE, AVERAGE_COLS, AVERAGE_ROWS, and FIXED models.
- Observed and null-network NODF calculations, with optional child-network NTC.
- Export of node and edge module assignments, statistics, and diagnostic records.

The default configuration uses 20 optimization trials and 100 replicates per null model. Module significance requires upper-tail Qb P values strictly below 0.05 under at least two null models. Root null testing during the observed calculation is disabled by default and can be configured or run separately.

For recursive analysis, prepare and run each selected child network, validate its results, and select eligible modules for the next level. The caller manages this iteration using the supplied workers and Python utilities.

| File | Purpose |
| --- | --- |
| `recursive_bimat_worker.m` | Analyze one network or module and export its partition and statistics. |
| `root_null_worker.m` | Test selected null models against a saved root partition, with checkpoints for completed replicates. |
| `core_rules.py` | Read networks, validate results, assess significance, and select child modules. |
| `root_null_exports.py` | Validate saved root labels and partitions. |
| `NullAdaptiveBrim.m`, `bimat_null_modularity.m` | Optimize null networks without final display sorting. |
| `bimat_nestedness_fast.m` | Calculate null-network NODF in blocks. |
| `bimat_validate_fast_null.m` | Check optimization and nestedness helpers against BiMat at runtime. |
| `examples/prepare_jobs.py` | Prepare JSON jobs for the MATLAB workers. |

Settings are stored in `analysis_config.json`. The `examples/network.tsv` file provides a small artificial network for trying the workflow.

### 2. Perturbation networks for robustness analysis

[`perturbation_OVPIN`](perturbation_OVPIN/README.md) contains three standalone generators:

| Script | Purpose |
| --- | --- |
| `expanding_host_ranges.py` | Add links between existing viruses and hosts under three taxonomic-distance scenarios. |
| `replacing_erroneous_links.py` | Replace host assignments for selected viruses under four taxonomic-distance scenarios. |
| `inferring_hosts_for_unrecalled_viruses.py` | Add viruses absent from the baseline with inferred host links under three cumulative host-pool scenarios. |

The scripts generate networks and associated metrics/audit outputs for downstream analysis. They do not run host prediction tools or automatically execute downstream modularity analysis.

Each generator accepts `--repeats` and `--seed-base`, with defaults of 100 and 1000, respectively. The number of generated networks depends on the number of scenarios and the requested repetitions. Perturbation magnitudes and additional inputs are specified through each script's command-line options.

### 3. Circular permutation tests for functional enrichment

[`circular_permutation_functional_enrichment`](circular_permutation_functional_enrichment/README.md) evaluates functional-gene capture by virus–host matching regions against 1,000 circular permutations by default.

Matching intervals are shifted together along each circularized contig while gene coordinates remain fixed. One relative shift is shared across contigs belonging to the same vOTU in each replicate. Analyses include direct CDS overlap and an extended definition with up to three flanking genes on each side, evaluated at contig and vOTU aggregation levels. Empirical tests use plus-one correction, with Benjamini–Hochberg adjustment across the shared family of testable category, region-definition, and aggregation-level combinations.

## Requirements

| Component | Requirements |
| --- | --- |
| Modularity core | Python 3.8+, MATLAB with JSON support, and an external BiMat installation |
| Perturbation generators | Python 3.8+; standard library only |
| Functional enrichment | Python 3.10+ and NumPy |

Use Python 3.10 or newer when using all three components in one environment. Install enrichment dependencies from the repository root:

```sh
python -m pip install -r circular_permutation_functional_enrichment/requirements.txt
```

BiMat is not bundled. Keep its version and settings consistent across comparisons; reference dependency hashes are recorded in `multilayer_modularity_core/CODE_INTEGRITY.json`.

## Getting started

Each component has its own README with input schemas, commands, parameter definitions, and output descriptions. Run component examples from the corresponding subdirectory.

For the modularity example:

```sh
cd multilayer_modularity_core
python examples/prepare_jobs.py module --bimat-path /path/to/BiMat
```

Then, with MATLAB's current directory set to `multilayer_modularity_core`, run:

```matlab
addpath(pwd);
recursive_bimat_worker(fullfile(pwd, 'example_output', 'observed', 'job.json'));
```

Replace `/path/to/BiMat` with your installation path. See the component README for separate root null tests and subsequent recursion levels.

For perturbation options:

```sh
cd perturbation_OVPIN
python expanding_host_ranges.py --help
python replacing_erroneous_links.py --help
python inferring_hosts_for_unrecalled_viruses.py --help
```

For functional enrichment, place the required study inputs in the component directory and supply your protein FASTA directory:

```sh
cd circular_permutation_functional_enrichment
python circular_permutation_functional_enrichment.py --faa-dirs prodigal_faa --workers 8
```

These commands assume a fresh terminal at the repository root for each component. Supply the required input files as described in the component READMEs.

## Inputs and reproducibility

Network analyses use binary host–virus edge lists with `pOTU` and `vOTU` identifiers. Perturbation analyses additionally require host taxonomy and, for newly inferred viruses, the supporting catalogues and prediction tables described in their README. Enrichment analyses require matching intervals, contig-to-vOTU assignments, selected high-quality contigs, functional annotations, and protein FASTA files with gene coordinates.

The repository contains analysis code and illustrative examples. Supply analysis datasets separately. For reproducibility, retain input versions, the code revision, software versions, parameter values, random seeds, and realized perturbation metrics.

## Outputs

| Component | Main outputs |
| --- | --- |
| Modularity core | Observed statistics in `result.json`; edge and node module assignments; ordered node labels; per-model statistics and null values; checkpoints and diagnostic records where applicable. |
| Perturbation generators | Perturbed edge lists, realized perturbation metrics, and method-specific audit files. The new-virus generator also exports added edges, degree summaries, and run metadata. |
| Functional enrichment | `D_global_random_tests.tsv`, `D_group_null_distributions.npz`, `D_function_boxplot.tsv`, `D_function_significance.tsv`, and `D_random_summary.tsv`. |

See the component READMEs for output schemas and interpretation.

## Citation

When using this code, cite the archived software release and its version-specific DOI once available. For analyses using BiMat, also cite:

Flores CO, Poisot T, Valverde S, et al. (2016). BiMat: a MATLAB package to facilitate the analysis of bipartite networks. *Methods in Ecology and Evolution*, 7(1), 127–132. https://doi.org/10.1111/2041-210X.12458

