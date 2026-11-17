# Multilayer bipartite modularity

Tools for detecting and testing hierarchical modules in bipartite networks using
[BiMat](https://bimat.github.io/). MATLAB functions calculate modularity,
nestedness and null-model statistics. Python functions validate results and
select modules for analysis at subsequent levels.

## Requirements

- Python 3.8 or later. Only the Python standard library is required.
- MATLAB with JSON support and BiMat installed.

BiMat is an external dependency. Its source files and global options affect
network optimization and null-model generation. Use a consistent BiMat version
when comparing analyses. Reference dependency hashes are recorded in
`CODE_INTEGRITY.json`.

## Main files

| File | Purpose |
| --- | --- |
| `recursive_bimat_worker.m` | Calculate the partition and statistics of one network or module. |
| `root_null_worker.m` | Test selected null models against a saved root partition. |
| `core_rules.py` | Read networks, validate outputs, assess significance and select child modules. |
| `root_null_exports.py` | Validate saved root labels and partition. |
| `NullAdaptiveBrim.m`, `bimat_null_modularity.m` | Optimize null networks without final display sorting. |
| `bimat_nestedness_fast.m` | Calculate NODF for null networks in blocks. |
| `bimat_validate_fast_null.m` | Check optimization and NODF helpers against BiMat at runtime. |
| `analysis_config.json` | Analysis settings. |
| `examples/prepare_jobs.py` | Prepare JSON inputs for the MATLAB functions. |
| `examples/network.tsv` | Artificial example network. |

## Network input

Provide a tab-separated edge list with the header `pOTU` and `vOTU`. The two
columns identify nodes in the two parts of the network:

```text
pOTU	vOTU
p1	v1
p1	v2
p2	v2
```

Node identifiers are strings without whitespace. Edges are binary; repeated
edges are removed, and normalized edges are sorted before calculation.

## Analysis settings

Edit `analysis_config.json` before preparing an analysis.

| Setting | Default | Meaning |
| --- | --- | --- |
| `optimization_trials` | `20` | Optimization trials per observed or null network. |
| `null_replicates` | `100` | Random networks per null model. |
| `alpha` | `0.05` | Significance threshold; the comparison is strict. |
| `min_significant_models` | `2` | Significant null models required for a modularity decision. |
| `run_root_null` | `false` | Include root null tests in the observed-network calculation and recursion decision. |
| `run_nestedness` | `true` | Calculate NODF. |
| `run_ntc` | `false` | Calculate optional child-network NTC when nestedness is enabled. |
| `max_depth` | `null` | Maximum recursion depth; `null` imposes no depth limit. |
| `seed` | `2026` | Base seed used to derive network-specific seeds. |

Keep all five supported models in `null_models`: `EQUIPROBABLE`, `AVERAGE`,
`AVERAGE_COLS`, `AVERAGE_ROWS` and `FIXED`. Selected root models can be specified
independently with the input-preparation script's `--models` option.

## Usage

Run the examples from the repository directory. Replace `../BiMat` with the
location of your BiMat installation. Relative paths supplied to the preparation
script are resolved against the repository directory.

### Calculate a network

Prepare an input using the supplied example network:

```sh
python examples/prepare_jobs.py module --bimat-path ../BiMat
```

Run the calculation in MATLAB:

```matlab
addpath(pwd);
recursive_bimat_worker(fullfile(pwd, 'example_output', 'observed', 'job.json'));
```

Use `--network` to supply another edge list and `--output` to choose a new output
directory. Preparation requires a new or empty directory. `--cpus` sets the
MATLAB thread limit.

### Test a saved root partition

After the root calculation finishes, prepare the required null models:

```sh
python examples/prepare_jobs.py root --bimat-path ../BiMat --models EQUIPROBABLE AVERAGE AVERAGE_COLS AVERAGE_ROWS FIXED
```

```matlab
addpath(pwd);
root_null_worker(fullfile(pwd, 'example_output', 'root_null', 'job.json'));
```

List one or more model names after `--models`. The script reads the saved root
labels, partition and observed statistics. Use `--observed` to select another
saved root directory, `--replicates` to override the number of random networks,
and `--output` to choose the destination.

To resume a root null calculation with identical inputs and settings, rerun the
MATLAB function with the same `job.json`. Successful replicates are recovered
from checkpoints. Use a separate output directory when changing the input or
analysis settings.

### Continue through subsequent levels

Validate a completed calculation and select eligible children in Python:

```python
from core_rules import read_json, validate_node_result, plan_children, write_network

config = read_json("analysis_config.json")
result, edges, groups = validate_node_result("example_output/observed", config)
node = {"node_id": "root", "level": 0}
records, children, stop_reason = plan_children(node, result, edges, groups, config)

for node_id, child_edges in children.items():
    write_network(f"example_output/child_{node_id}.tsv", child_edges)
```

Prepare each selected child with the `module` command, supplying its `--network`,
`--node-id`, `--level` and a distinct `--output`. Also supply `--root-network`
with the root calculation's normalized `network.tsv`: the exact file bytes are
used to derive child seeds. At deeper levels, node IDs describe the module path,
such as `1_2` at level 2.

Run the MATLAB worker for each prepared child, then validate and select its
children using that child's node ID and level. Callers manage this iteration;
`plan_children` returns selections without executing further calculations.

## Module selection and significance

A tested network is significant when at least `min_significant_models` null
models have Qb upper-tail p values strictly below `alpha`. Under the defaults,
this requires at least two of the five models to have `p_upper < 0.05`.

Child modules are selected when:

1. The parent satisfies the significance rule. With `run_root_null: false`, the
   initial root can be split without a root significance decision.
2. The partition contains at least two modules with internal edges.
3. The child has fewer edges than its parent and at least
   `sqrt(2 * parent_edge_count)` internal edges.
4. The child is within `max_depth`, when a depth limit is set.

There is no minimum row-count or column-count filter. Cross-module edges are
recorded with module ID `-1` and are excluded from child networks.

Failed null tests remain unknown votes. The significance decision is `NA` when
the unresolved votes could change the outcome. For single-row or single-column
networks, the FIXED generator returns the unchanged matrix because the degree
sequence determines all edges. A separately calculated subset of root models
provides model-specific statistics; it does not replace the five-model recursion
decision.

The Qb p value uses an upper-tail Monte Carlo test with a plus-one correction:

```text
tolerance = 8 * eps(max(1, abs(observed_Qb)))
p_upper = (1 + count(null_Qb >= observed_Qb - tolerance)) / (B + 1)
```

Here, `B` is the number of null replicates. Observed and null networks use the
same optimization-trial count. NODF is tested separately when enabled. Qr is
reported as a descriptive statistic without a Qr p value. Optional nestedness
failures are recorded as diagnostics and do not invalidate valid Qb results.

## Outputs

| Output | Contents |
| --- | --- |
| `result.json` | Observed statistics, null-model summaries and statuses. |
| `Module.txt` | Edge assignments, including cross-module edges. |
| `node_modules.tsv` | Node-to-module assignments. |
| `label_row.txt`, `label_col.txt` | Ordered labels for the two node sets. |
| `result_null_model.tsv` | Per-model statistics; `p_upper` tests Qb and `NODF_p_upper` tests NODF. |
| `null_values.tsv` | One row per completed random network, with seed, Qb, Qr and NODF. |
| `bp_for_null.mat` | BiMat object saved by the network worker. |
| `checkpoints/` | Completed replicates saved by the root null worker. |
| `null_model_failures.tsv`, `nestedness_errors.tsv` | Diagnostic records. |

Output availability depends on the calculation type and enabled settings.

## Reference

This project uses BiMat for bipartite network analysis. Please cite the BiMat
publication when using these analyses:

Flores C O, Poisot T, Valverde S, et al. BiMat: a MATLAB package to facilitate
the analysis of bipartite networks[J]. *Methods in Ecology and Evolution*,
2016, **7**(1): 127–132. [doi:10.1111/2041-210X.12458](https://doi.org/10.1111/2041-210X.12458).

BiMat documentation and source: [https://bimat.github.io/](https://bimat.github.io/).
