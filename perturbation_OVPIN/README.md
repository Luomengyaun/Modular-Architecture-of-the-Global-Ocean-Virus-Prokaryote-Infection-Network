# OVPIN perturbation networks

Three standalone Python scripts generate perturbed host–virus networks for sensitivity analyses. Networks are represented as unique `pOTU`–`vOTU` edges, where `pOTU` identifies a host and `vOTU` identifies a virus. The baseline network is the reference for defining retained, added, and replaced links.

| Script | Perturbation | Main control | Scenarios per replicate |
|---|---|---|---|
| `expanding_host_ranges.py` | Add links between existing hosts and viruses | Added edges as a percentage of baseline edges | 3 taxonomic distances |
| `replacing_erroneous_links.py` | Replace original hosts of selected viruses | Percentage of viruses retaining original hosts | 4 taxonomic distances |
| `inferring_hosts_for_unrecalled_viruses.py` | Add viruses absent from the baseline, with inferred host links | Added viruses as a percentage of baseline viruses, or an exact count | 3 cumulative host pools |

Each script generates all of its scenarios automatically. `--addition` has different meanings in the two scripts that accept it; see the method descriptions below.

## Requirements

- Python 3.8 or later.
- Python standard library only; no additional packages are required.
- Input networks and taxonomy tables supplied by the reader. The scripts do not download data or run host prediction tools.

Run commands from the directory containing the scripts. The examples use `python`; use the command for Python 3 in your environment. Display all options with, for example, `python expanding_host_ranges.py --help`.

## Input files

Use UTF-8, tab-separated text files. The examples below form a small synthetic dataset for trying all three scripts. Create a `data` directory and save each block under the indicated filename. These identifiers and taxonomic labels are illustrative, not research observations.

### Baseline network: `data/baseline.tsv`

```tsv
pOTU	vOTU
p1	v1
p2	v2
p3	v3
p4	v4
p5	v5
p6	v6
p7	v7
p8	v8
```

For `expanding_host_ranges.py` and `replacing_erroneous_links.py`, the header must contain the case-sensitive column names `pOTU` and `vOTU`. Additional columns are ignored. Duplicate edges are counted once; blank lines and rows with empty edge identifiers are skipped.

`inferring_hosts_for_unrecalled_viruses.py` also accepts a headerless file whose first two columns are `pOTU` and `vOTU`. Its named headers are case-insensitive, and empty edge identifiers are errors. A four-column baseline with embedded viral taxonomy is accepted, but embedded group labels are used only for diagnostics; the supplied current `vOTU`–`gvOTU` map is authoritative.

Only nodes appearing in an edge are represented. Isolated nodes cannot be specified in these edge lists.

### Host taxonomy: `data/taxonomy.tsv`

```tsv
pOTU	Order	Family	Genus
p1	O1	F1	G1
p2	O1	F1	G1
p3	O1	F1	G2
p4	O1	F1	G2
p5	O1	F2	G3
p6	O1	F2	G3
p7	O2	F3	G4
p8	O2	F3	G4
```

For the first two scripts, the header must contain `pOTU`, `Order`, `Family`, and `Genus` with exactly this capitalization. They restrict candidate hosts to pOTUs already present in the baseline, even if the taxonomy table contains additional hosts. A baseline host missing from the taxonomy is reported and cannot supply the missing taxonomic information for constrained sampling.

For `inferring_hosts_for_unrecalled_viruses.py`, supply the **complete OPG host taxonomy**, including hosts outside the baseline when available. All baseline pOTUs must be present. The named four-column format above is supported with case-insensitive headers. A headerless OPG table is also supported, with `pOTU`, order, family, and genus in columns 1, 6, 7, and 8, respectively (the expected source format has nine columns). The script rejects conflicting records and contradictory known genus–family, genus–order, or family–order assignments.

Empty values and the case-insensitive labels `na`, `n/a`, `none`, `null`, `nan`, `unclassified`, and `unknown` represent missing taxonomy. Missing labels do not establish a shared taxonomic group. Use a coherent taxonomy and consistent names across files. Unprefixed labels work with all three scripts; the new-virus inference script additionally strips `o__`, `f__`, and `g__` rank prefixes while retaining meaningful suffixes such as `_A`.

### Viral group map: `data/votu_gvotu.tsv`

Required only for `inferring_hosts_for_unrecalled_viruses.py`. Supply the complete current OVG mapping, including baseline viruses and candidate new viruses. `gvOTU` is the viral group used to share host-genus evidence.

```tsv
vOTU	gvOTU
v1	gv1
v2	gv1
v3	gv2
v4	gv2
v5	gv3
v6	gv3
v7	gv4
v8	gv4
v9	gv1
v10	gv2
```

The two-column header is optional. Conflicting assignments for the same vOTU are errors. Unknown group labels cannot link viruses through a shared group. Baseline viruses absent from the map cause an error unless `--allow-missing-baseline-groups` is supplied; that option retains their original edges while excluding them as group-based donors.

### Predicted host genera: `data/iphop.tsv`

Required only for `inferring_hosts_for_unrecalled_viruses.py`. Each row contains a vOTU and a Python dictionary literal whose string keys are predicted host genera.

```tsv
vOTU	host_genera
v9	{'G1': 0.95}
v10	{'G2': 0.90}
```

The header is optional. Dictionary values are ignored: this script uses genus membership, not prediction scores or a score threshold. Apply any required confidence filtering before supplying the file. Duplicate vOTU rows are merged by taking the union of genus keys. Predictions are filtered to known genera in the supplied OPG taxonomy and to viruses in the supplied OVG map before group-based propagation. Unmatched genera are reported in the audit outputs. An empty dictionary is allowed.

## Expanding the host ranges of existing viruses

`expanding_host_ranges.py` preserves every original edge and adds unique links between the existing baseline host and virus nodes. Taxonomic conditions are evaluated against **all original hosts of the target virus**, rather than the hosts added during the simulation.

| Mode | Eligible new host |
|---|---|
| `same_genus_diff_potu` | An unlinked pOTU whose known genus is represented among the original hosts |
| `same_family_diff_genus` | A host in a represented family, with a known genus absent from all original hosts |
| `same_order_diff_family` | A host in a represented order, with a known family absent from all original hosts |

At each step, an eligible virus is selected uniformly, then an eligible new host is sampled for that virus. This does not sample uniformly from the entire network's set of missing edges. The requested number of edges is added separately in each mode. If a mode lacks sufficient candidates, the script raises an error rather than relaxing its taxonomic constraint.

```shell
python expanding_host_ranges.py --fbaseline data/baseline.tsv --ftaxa data/taxonomy.tsv --addition 25 --outfile results/host_ranges --repeats 2 --seed-base 1000
```

`--addition 25` adds `floor(0.25 × baseline_edge_count + 0.5)` edges per network: two edges in the synthetic example. The command creates six networks, with ten edges each.

`--max-retry` defaults to `500`. It controls the number of rejection-sampling attempts before an exact search of remaining candidates. It does not relax the candidate definition.

## Replacing erroneous links

`replacing_erroneous_links.py` preserves the virus nodes, each virus's degree, and the total edge count. For selected viruses, it replaces all original hosts with distinct hosts that were not originally linked to that virus. Candidate hosts come from the baseline host set. Host degrees and the set of hosts still connected in the final network can change.

The first three modes use the taxonomic distances listed above, but each replacement is compared with the **individual original host being replaced**. A fourth mode, `cross_order`, requires the replacement host to have a known order different from that original host's order. These are simulation constraints; they do not establish whether a link is biologically correct or incorrect.

```shell
python replacing_erroneous_links.py --fbaseline data/baseline.tsv --ftaxa data/taxonomy.tsv --accuracy 75 --outfile results/replaced --repeats 2 --seed-base 1000
```

`--accuracy` is a percentage from `0` to `100`. The target number of accurate viruses is `floor(accuracy / 100 × baseline_virus_count + 0.5)`. Here, six of the eight viruses retain their original hosts, and two have all original hosts replaced. The command creates eight networks, with eight edges each.

For each replicate, the same viruses remain accurate across all four modes. A virus that cannot have all its hosts replaced in any one mode must remain accurate in every mode. Requests below this shared minimum accuracy are rejected. Replacement matching ensures unique hosts within each virus's final host set.

The reported metrics are defined as:

- **Accuracy:** the fraction of viruses retaining at least one original host.
- **Consistency:** the mean, over viruses, of the fraction of original hosts retained by each virus.

In this implementation, a virus retains either all or none of its original hosts. Consequently, realized consistency equals realized accuracy; they are not independent controls. These are virus-level summaries, not the fraction of original edges retained across the whole network. Requested and realized values are written as fractions in the metrics files.

`--max-retry` is retained for command-line compatibility and is unused by this replacement method.

## Inferring hosts for viruses absent from the baseline

`inferring_hosts_for_unrecalled_viruses.py` preserves the baseline and adds eligible viruses from the complete OVG map. A new virus's reference host genera are the union of:

1. Known host genera of baseline viruses in the same known gvOTU.
2. Its own filtered iPHoP genus predictions.
3. Filtered iPHoP genus predictions from other viruses in the same known gvOTU.

Only viruses absent from the baseline and having at least one same-genus candidate are eligible for the paired simulations. An unknown gvOTU does not enable group-based sharing, although direct predictions can still provide reference genera.

| Mode | Candidate host pool from the complete OPG taxonomy |
|---|---|
| `same_genus` | Hosts in the reference genera |
| `same_family` | Hosts in families represented by the reference genera, retaining known same-genus candidates whose family is missing |
| `same_order` | Hosts in orders represented by the reference genera or their known families, retaining narrower candidates whose broader ranks are missing |

These pools are cumulative: genus candidates are included in the family pool, and family candidates are included in the order pool. Unlike the host-range expansion modes, broader pools do **not** require a different genus or family. Sampled edge sets need not be nested, because each mode samples hosts separately.

```shell
python inferring_hosts_for_unrecalled_viruses.py --fbaseline data/baseline.tsv --ftaxa data/taxonomy.tsv --fvotu-gvotu data/votu_gvotu.tsv --fiphop data/iphop.tsv --addition 25 --outfile results/new_viruses --repeats 2 --seed-base 1000
```

Here, `--addition 25` adds `floor(0.25 × baseline_virus_count + 0.5)` **viruses**, not edges. The synthetic example adds both eligible new viruses (`v9` and `v10`), producing six networks with ten viruses and ten edges each.

Use `--n-new-votus N` instead of `--addition` to request an exact number. The two options are mutually exclusive. Omit both to use all eligible new viruses. A count exceeding the eligible library size is rejected.

Within a replicate, all three modes use the same selected viruses and the same assigned degree for each new virus. Target degree frequencies approximate the baseline virus-degree distribution using largest-remainder allocation. Degrees are assigned in descending order subject to same-genus capacity; infeasible quotas are reduced and the deviations are reported. Each mode then samples distinct hosts for each new virus from its candidate pool. Hosts outside the baseline may be added from the complete taxonomy. Connectance is measured using the actual connected host and virus nodes; it is not held fixed.

Add `--prepare-only` to write the candidate library and input audit without generating replicate networks. This reserves the chosen output prefix; use a new prefix for a later full run.

## Output files

`--outfile` is an output **prefix**, not a single filename. Parent directories are created as needed. All network files contain a `pOTU`/`vOTU` header and sorted, unique tab-separated edges.

### Existing-virus host-range expansion

For a prefix such as `results/host_ranges` and each of the three modes:

```text
results/host_ranges_<mode>/host_ranges_<mode>_rep001.txt
results/host_ranges_<mode>_bimat/real_add_edges/host_ranges_<mode>_rep001.txt
```

The first file is a network; the second is its metrics table. Metrics include requested and realized additions, edge counts, the number of viruses gaining hosts, and taxonomic summaries. Addition rates are stored as fractions. Cross-rank edge summaries can overlap and should not be summed as mutually exclusive categories.

### Link replacement

For a prefix such as `results/replaced` and each of the four modes:

```text
results/replaced_<mode>/replaced_<mode>_rep001.txt
results/replaced_<mode>_bimat/real_acc_con/replaced_<mode>_rep001.txt
```

The metrics table records requested and realized accuracy, the accurate-virus count, realized consistency, and counts of removed and replacement edges. In both scripts above, the `fout` column records the network filename without its directory. The `_bimat` directory names are output conventions; no additional software is required to generate these files.

These two scripts overwrite files with matching names when a prefix is reused. Use separate prefixes for independent runs. A failed run can leave outputs from earlier modes or replicates.

### New-virus host inference

For a prefix such as `results/new_viruses`:

| Output | Contents |
|---|---|
| `new_viruses_<mode>/new_viruses_<mode>_rep001.txt` | Baseline plus the inferred new-virus edges |
| `new_viruses_<mode>/new_viruses_<mode>_rep001_added_edges.tsv` | Added edges only |
| `new_viruses_paired/rep001_degree_quotas.tsv` | Per-virus target degree, assigned degree, and same-genus capacity |
| `new_viruses_paired/rep001_degree_distribution.tsv` | Baseline, target, and assigned degree distributions |
| `new_viruses_audit/input_summary.json` | Input counts, filtering diagnostics, and simulation rules |
| `new_viruses_audit/new_virus_library.tsv` | Candidate viruses, evidence sources, candidate capacities, and eligibility |
| `new_viruses_audit/baseline_virus_degree_distribution.tsv` | Baseline virus-degree distribution |
| `new_viruses_audit/iphop_filtered_for_OVG.txt` | Genus predictions retained after OPG and OVG filtering |
| `new_viruses_audit/iphop_unmatched_genera.tsv` | Predicted genera absent from the supplied taxonomy |
| `new_viruses_metrics.tsv` | One metrics row per mode and replicate |
| `new_viruses_run.json` | Parameters, run status, and output references |

Metrics include added node/edge counts, connectance, degree summaries, and capacity reductions. `TV` is total variation distance between normalized degree distributions; `KS` is their maximum cumulative-distribution difference. These are descriptive distances, not hypothesis-test results. Distances involving an empty new-virus distribution are left blank.

This script refuses to reuse any existing output associated with the prefix. Its records store input filenames and output references relative to the output parent directory. Retain your input versions separately when archiving an analysis.

## Reproducibility and interpretation

All scripts accept `--seed-base` (default `1000`) and `--repeats` (default `100`). The default generates 300 networks for either three-mode script and 400 networks for the four-mode replacement script, in addition to metrics and audit files.

Use the same code version, Python version, input files, and parameter values to reproduce sampled networks. Candidate order is stabilized before seeded random choices, so Python's hash randomization does not select a different starting order. Changing inputs, modes, sampling options, or the Python version can change outputs even when the base seed is unchanged. The new-virus script records separate seeds for virus selection, degree quotas, degree assignment, and host sampling.

The methods provide controlled perturbations of a reference network. They do not validate host predictions or establish that absent links are biologically impossible. Taxonomic completeness and candidate availability affect feasible perturbations and should be considered alongside the reported realized metrics. In the new-virus method, matching the baseline degree distribution is a modeling assumption; finite sample sizes and host capacity can produce deviations. Report the chosen scheme, input catalogue versions, parameter values, seeds, and realized perturbation measures with downstream results.
