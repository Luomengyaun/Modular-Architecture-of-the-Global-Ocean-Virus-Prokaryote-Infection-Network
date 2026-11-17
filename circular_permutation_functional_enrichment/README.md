# Circular permutation test for functional enrichment of virus–host matching regions

This script reproduces the functional enrichment analysis (panel D) for generalist viruses in the Ocean Viral–Prokaryotic Interaction Network (OVPIN). It compares the observed capture of functional genes in host-matching regions with a null distribution of **1,000 circular permutations**.

## Requirements

Python 3.10 or newer and NumPy. Run from the folder containing this script.

```sh
python -m pip install -r requirements.txt
```

## Inputs

| Input | Format |
| --- | --- |
| `pair_union_le10pct.tsv` | Filtered virus–host pairs. Required columns: `contig`, `host`, `merged_intervals`, `union_length_bp`, `contig_length_bp`, `union_coverage_pct`. Coordinates are 1-based and inclusive. Intervals are written as `start-end`, separated by semicolons. Each pair must cover at most 10% of the contig. |
| `retained_contig_votu.tsv` | Tab-separated `contig` and `vOTU` columns. |
| `matched_HQ_comp90_contigs.txt` | Selected high-quality contigs, one contig ID per line. |
| `retained_contig_votu_phrokka.txt` | Pharokka annotation table. Protein IDs are in column 1 and functional categories in column 5. |
| Prodigal protein FASTA files | Protein headers must have the form `>contig_1 # start # end # strand # ...`. Protein IDs must match the annotation table. Supply one or more directories with `--faa-dirs`, individual files with `--extra-faa`, or both. |

Input identifiers should use the same original contig IDs throughout. The repository contains the analysis code. Supply the corresponding study input files separately.

## Run the analysis

With the four default input filenames in the current folder and the FASTA files in `prodigal_faa/`:

```sh
python circular_permutation_functional_enrichment.py --faa-dirs prodigal_faa --workers 8
```

For other filenames, pass `--pairs`, `--retained`, `--hq-contigs` and `--annotations`. The defaults are `--n-random 1000`, `--seed 20261111` and `--batch-size 128`. Output files are written to `functional_enrichment_results/` by default.

## Permutation and inference

1. Merge qualifying host-matching intervals across hosts for each selected contig. Calculate the percentage of genes in each functional category captured by the resulting intervals, counting each gene once per selection. The denominator is the number of genes of that category on the contig. Evaluate direct CDS overlap and a second definition that includes up to three flanking genes on each side.
2. Keep gene coordinates fixed. In each replicate, draw one uniform relative shift per vOTU and convert it to an integer offset in base pairs for each member contig. Move the matching intervals together on each contig, wrapping at the sequence boundary. This preserves interval lengths and spacing on the circle.
3. For each category, calculate the mean capture percentage across eligible contigs. For the vOTU-level statistic, first average eligible contigs within each vOTU, then give each eligible vOTU equal weight. Units without genes of the category are excluded for that category.
4. For each statistic, calculate the upper and lower empirical tail probabilities as `(count + 1) / (permutations + 1)`. The two-sided P value is twice the smaller tail probability, capped at 1. Apply Benjamini–Hochberg correction jointly across the testable combinations of 10 categories, two region definitions, and two aggregation levels. The panel D figure shows direct overlap at the vOTU level and uses the adjusted P values from that shared family.

## Outputs

| File | Contents |
| --- | --- |
| `D_global_random_tests.tsv` | Results across both region definitions and levels, including observed and null capture, empirical P values, and the global BH adjustment. |
| `D_group_null_distributions.npz` | Permutation averages for each level, region definition and functional category. The zero-based category order follows `FUNCTIONS` in the script. |
| `D_function_boxplot.tsv` | Per-vOTU observed and permutation-mean values for direct overlap. `Random` and `OVIPN` are the null and observed series, respectively. |
| `D_function_significance.tsv` | Direct-overlap vOTU tests used in panel D. |
| `D_random_summary.tsv` | Input filenames, sample sizes, permutation count and seed. |
