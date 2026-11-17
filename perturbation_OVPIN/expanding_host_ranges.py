#!/usr/bin/env python3

import argparse
import math
import random
from collections import defaultdict
from pathlib import Path


# Three missing-edge scenarios used in this analysis:
#   1. same_genus_diff_potu: same genus, different pOTU
#   2. same_family_diff_genus: same family, different genus
#   3. same_order_diff_family: same order, different family
ALL_MODES = [
    "same_genus_diff_potu",
    "same_family_diff_genus",
    "same_order_diff_family",
]

MODE_LABELS = {
    "same_genus_diff_potu": "same genus, different pOTU",
    "same_family_diff_genus": "same family, different genus",
    "same_order_diff_family": "same order, different family",
}

MISSING_TAXA = {
    "",
    "na",
    "n/a",
    "none",
    "null",
    "nan",
    "unclassified",
    "unknown",
}


def round_half_up(value):
    """Round halves upward instead of using Python's bankers' rounding."""
    return int(math.floor(value + 0.5))


def normalize_taxon(value):
    value = value.strip()
    if value.lower() in MISSING_TAXA:
        return None
    return value


def split_line(line, use_tabs):
    """Split a tab-delimited file, with whitespace fallback for simple tables."""
    line = line.rstrip("\n\r")
    return line.split("\t") if use_tabs else line.split()


def read_baseline_network(filename):
    """Read a network containing columns named pOTU and vOTU."""
    positive_hosts = defaultdict(set)
    baseline_edges = set()

    with open(filename, "r", encoding="utf-8") as handle:
        first_line = handle.readline()
        if not first_line:
            raise ValueError("The baseline network file is empty.")

        use_tabs = "\t" in first_line
        header = split_line(first_line, use_tabs)
        if "pOTU" not in header or "vOTU" not in header:
            raise ValueError(
                "The baseline network must contain columns named 'pOTU' and 'vOTU'."
            )

        p_idx = header.index("pOTU")
        v_idx = header.index("vOTU")
        required_idx = max(p_idx, v_idx)

        for line_number, line in enumerate(handle, start=2):
            if not line.strip():
                continue

            fields = split_line(line, use_tabs)
            if len(fields) <= required_idx:
                raise ValueError(
                    f"Line {line_number} in {filename!r} has too few columns."
                )

            p = fields[p_idx].strip()
            v = fields[v_idx].strip()
            if not p or not v:
                continue

            edge = (p, v)
            if edge not in baseline_edges:
                baseline_edges.add(edge)
                positive_hosts[v].add(p)

    if not baseline_edges:
        raise ValueError("The baseline network contains no valid edges.")

    return positive_hosts, baseline_edges


def read_potu_taxonomy(filename, baseline_potus):
    """Read pOTU taxonomy. Order, Family and Genus are required."""
    required_columns = ["pOTU", "Order", "Family", "Genus"]
    taxonomy = {}

    with open(filename, "r", encoding="utf-8") as handle:
        first_line = handle.readline()
        if not first_line:
            raise ValueError("The taxonomy file is empty.")

        use_tabs = "\t" in first_line
        header = split_line(first_line, use_tabs)
        missing_columns = [c for c in required_columns if c not in header]
        if missing_columns:
            raise ValueError(
                "The taxonomy file is missing required columns: "
                + ", ".join(missing_columns)
            )

        indexes = {column: header.index(column) for column in required_columns}
        required_idx = max(indexes.values())

        for line_number, line in enumerate(handle, start=2):
            if not line.strip():
                continue

            fields = split_line(line, use_tabs)
            if len(fields) <= required_idx:
                raise ValueError(
                    f"Line {line_number} in {filename!r} has too few columns."
                )

            p = fields[indexes["pOTU"]].strip()
            if not p:
                continue

            record = {
                "Order": normalize_taxon(fields[indexes["Order"]]),
                "Family": normalize_taxon(fields[indexes["Family"]]),
                "Genus": normalize_taxon(fields[indexes["Genus"]]),
            }

            if p in taxonomy and taxonomy[p] != record:
                raise ValueError(
                    f"Conflicting taxonomy records were found for pOTU {p!r}."
                )
            taxonomy[p] = record

    missing_potus = sorted(set(baseline_potus) - set(taxonomy))
    if missing_potus:
        preview = ", ".join(missing_potus[:10])
        suffix = "" if len(missing_potus) <= 10 else ", ..."
        print(
            f"Warning: {len(missing_potus)} baseline pOTUs have no taxonomy record: "
            f"{preview}{suffix}"
        )

    return taxonomy


class TaxonomyIndex:
    """Indexes used to sample taxonomy-constrained missing edges."""

    def __init__(self, all_potus, taxonomy):
        self.all_potus = tuple(sorted(all_potus))
        self.taxonomy = taxonomy

        genus_groups = defaultdict(list)
        family_groups = defaultdict(list)
        order_groups = defaultdict(list)
        ordered_potus = []

        for p in self.all_potus:
            record = taxonomy.get(p)
            if record is None:
                continue

            genus = record["Genus"]
            family = record["Family"]
            order = record["Order"]

            if genus is not None:
                genus_groups[genus].append(p)

            # For the family mode, candidate genus must be known.
            if family is not None and genus is not None:
                family_groups[family].append(p)

            # For the order mode, candidate family must be known.
            if order is not None and family is not None:
                order_groups[order].append(p)

            if order is not None:
                ordered_potus.append(p)

        self.genus_groups = {
            key: tuple(sorted(value)) for key, value in genus_groups.items()
        }
        self.family_groups = {
            key: tuple(sorted(value)) for key, value in family_groups.items()
        }
        self.order_groups = {
            key: tuple(sorted(value)) for key, value in order_groups.items()
        }
        self.ordered_potus = tuple(sorted(ordered_potus))
        self.genus_sizes = {key: len(value) for key, value in self.genus_groups.items()}
        self.family_sizes = {key: len(value) for key, value in self.family_groups.items()}
        self.order_sizes = {key: len(value) for key, value in self.order_groups.items()}

    def host_taxa(self, hosts):
        genera = set()
        families = set()
        orders = set()

        for p in hosts:
            record = self.taxonomy.get(p)
            if record is None:
                continue
            if record["Genus"] is not None:
                genera.add(record["Genus"])
            if record["Family"] is not None:
                families.add(record["Family"])
            if record["Order"] is not None:
                orders.add(record["Order"])

        return genera, families, orders

    def relation_is_valid(self, mode, original_hosts, candidate):
        """
        Test whether candidate expands host range at the requested rank.

        Definitions are relative to all original hosts of the vOTU:
          same_genus_diff_potu:
              candidate genus is already represented among original hosts.
          same_family_diff_genus:
              candidate family is represented, but candidate genus is new.
          same_order_diff_family:
              candidate order is represented, but candidate family is new.
          cross_order:
              candidate order is new.
        """
        if candidate in original_hosts:
            return False

        record = self.taxonomy.get(candidate)
        if record is None:
            return False

        genera, families, orders = self.host_taxa(original_hosts)

        if mode == "same_genus_diff_potu":
            return record["Genus"] is not None and record["Genus"] in genera

        if mode == "same_family_diff_genus":
            return (
                record["Family"] is not None
                and record["Family"] in families
                and record["Genus"] is not None
                and record["Genus"] not in genera
            )

        if mode == "same_order_diff_family":
            return (
                record["Order"] is not None
                and record["Order"] in orders
                and record["Family"] is not None
                and record["Family"] not in families
            )

        if mode == "cross_order":
            return record["Order"] is not None and record["Order"] not in orders

        raise ValueError(f"Unknown mode: {mode}")

    def candidate_count(self, mode, original_hosts):
        """Count all unique pOTUs that can be added for one vOTU."""
        genera, families, orders = self.host_taxa(original_hosts)

        if mode == "same_genus_diff_potu":
            total = sum(self.genus_sizes.get(g, 0) for g in genera)
            existing_in_pool = sum(
                1
                for p in original_hosts
                if self.taxonomy.get(p, {}).get("Genus") in genera
            )
            return max(0, total - existing_in_pool)

        if mode == "same_family_diff_genus":
            total = sum(self.family_sizes.get(f, 0) for f in families)
            excluded = sum(self.genus_sizes.get(g, 0) for g in genera)
            return max(0, total - excluded)

        if mode == "same_order_diff_family":
            # Count within the actual order pools: family_sizes excludes
            # hosts with unknown genus, which are valid in the order mode.
            return sum(
                p not in original_hosts
                and self.taxonomy[p]["Family"] not in families
                for order in orders
                for p in self.order_groups.get(order, ())
            )

        if mode == "cross_order":
            excluded = sum(self.order_sizes.get(o, 0) for o in orders)
            return max(0, len(self.ordered_potus) - excluded)

        raise ValueError(f"Unknown mode: {mode}")

    @staticmethod
    def choose_from_groups(groups, rng):
        """Choose uniformly from the concatenation of several groups."""
        total = sum(len(group) for group in groups)
        if total == 0:
            return None

        index = rng.randrange(total)
        for group in groups:
            if index < len(group):
                return group[index]
            index -= len(group)

        raise RuntimeError("Internal weighted-group sampling error.")

    def choose_candidate(
        self,
        mode,
        original_hosts,
        assigned_hosts,
        rng,
        max_retry=500,
    ):
        """Select one valid candidate that is not already assigned."""
        genera, families, orders = self.host_taxa(original_hosts)

        if mode == "same_genus_diff_potu":
            groups = [self.genus_groups.get(g, ()) for g in sorted(genera)]
            source = None
        elif mode == "same_family_diff_genus":
            groups = [self.family_groups.get(f, ()) for f in sorted(families)]
            source = None
        elif mode == "same_order_diff_family":
            groups = [self.order_groups.get(o, ()) for o in sorted(orders)]
            source = None
        elif mode == "cross_order":
            source = self.ordered_potus
            groups = None
        else:
            raise ValueError(f"Unknown mode: {mode}")

        # Fast rejection sampling.
        for _ in range(max_retry):
            if groups is None:
                if not source:
                    return None
                candidate = rng.choice(source)
            else:
                candidate = self.choose_from_groups(groups, rng)
                if candidate is None:
                    return None

            if candidate in assigned_hosts:
                continue
            if self.relation_is_valid(mode, original_hosts, candidate):
                return candidate

        # Exact fallback. This is slower but guarantees correctness when
        # the valid fraction of a taxonomic group is very small.
        valid = []
        if groups is None:
            iterable = source
        else:
            iterable = (p for group in groups for p in group)

        seen = set()
        for candidate in iterable:
            if candidate in seen:
                continue
            seen.add(candidate)
            if candidate in assigned_hosts:
                continue
            if self.relation_is_valid(mode, original_hosts, candidate):
                valid.append(candidate)

        if not valid:
            return None
        return rng.choice(valid)


class DynamicVotuPool:
    """O(1) random choice and O(1) removal for eligible vOTUs."""

    def __init__(self, items):
        self.items = list(items)
        self.positions = {item: i for i, item in enumerate(self.items)}

    def __len__(self):
        return len(self.items)

    def choice(self, rng):
        if not self.items:
            return None
        return self.items[rng.randrange(len(self.items))]

    def remove(self, item):
        position = self.positions.pop(item, None)
        if position is None:
            return

        last = self.items.pop()
        if position < len(self.items):
            self.items[position] = last
            self.positions[last] = position


def count_cross_rank_votus(hosts_by_votu, taxonomy, rank):
    count = 0
    for hosts in hosts_by_votu.values():
        taxa = {
            taxonomy[p][rank]
            for p in hosts
            if p in taxonomy and taxonomy[p][rank] is not None
        }
        if len(taxa) >= 2:
            count += 1
    return count


def classify_added_edges(added_edges, original_hosts, taxonomy):
    """Summarize the taxonomic distance of newly added edges."""
    same_genus = 0
    cross_genus = 0
    cross_family = 0
    cross_order = 0
    unknown = 0

    for p, v in added_edges:
        candidate = taxonomy.get(p)
        if candidate is None:
            unknown += 1
            continue

        original_records = [
            taxonomy[h]
            for h in original_hosts[v]
            if h in taxonomy
        ]

        genera = {r["Genus"] for r in original_records if r["Genus"] is not None}
        families = {r["Family"] for r in original_records if r["Family"] is not None}
        orders = {r["Order"] for r in original_records if r["Order"] is not None}

        if candidate["Genus"] is not None and candidate["Genus"] in genera:
            same_genus += 1
        else:
            cross_genus += 1

        if candidate["Family"] is not None and candidate["Family"] not in families:
            cross_family += 1

        if candidate["Order"] is not None and candidate["Order"] not in orders:
            cross_order += 1

    return {
        "same_genus_added_edges": same_genus,
        "cross_genus_added_edges": cross_genus,
        "cross_family_added_edges": cross_family,
        "cross_order_added_edges": cross_order,
        "unknown_taxonomy_added_edges": unknown,
    }


def add_missing_edges_by_taxonomy(
    fbaseline,
    ftaxa,
    fout,
    fm,
    addition=0.20,
    mode="same_genus_diff_potu",
    seed=123,
    max_retry=500,
    prepared=None,
):
    """
    Add missing virus-host edges to an existing OVPIN.

    Parameters
    ----------
    addition:
        Number of added edges as a fraction of the original edge count.
        For example, 0.20 adds round(0.20 * baseline_edges) new edges.

    mode:
        same_genus_diff_potu
        same_family_diff_genus
        same_order_diff_family

    The function preserves all original edges and nodes. It only adds unique
    pOTU-vOTU edges. A vOTU is selected uniformly among vOTUs that still have
    at least one eligible candidate host under the requested mode.
    """
    if addition < 0:
        raise ValueError("addition must be greater than or equal to zero.")
    if mode not in ALL_MODES:
        raise ValueError(
            f"Unknown mode {mode!r}. Expected one of: {', '.join(ALL_MODES)}"
        )
    if max_retry < 1:
        raise ValueError("max_retry must be at least 1.")

    rng = random.Random(seed)

    if prepared is None:
        positive_hosts, baseline_edges = read_baseline_network(fbaseline)
        all_votus = sorted(positive_hosts)
        all_potus = sorted({p for p, _ in baseline_edges})
        taxonomy = read_potu_taxonomy(ftaxa, all_potus)
        index = TaxonomyIndex(all_potus, taxonomy)
        capacity_by_mode = {}
    else:
        positive_hosts = prepared["positive_hosts"]
        baseline_edges = prepared["baseline_edges"]
        all_votus = prepared["all_votus"]
        all_potus = prepared["all_potus"]
        taxonomy = prepared["taxonomy"]
        index = prepared["index"]
        capacity_by_mode = prepared["capacity_by_mode"]

    n_baseline_edges = len(baseline_edges)
    requested_added_edges = round_half_up(addition * n_baseline_edges)

    if mode not in capacity_by_mode:
        capacity_by_mode[mode] = {
            v: index.candidate_count(mode, positive_hosts[v])
            for v in all_votus
        }
    capacity = capacity_by_mode[mode]
    maximum_added_edges = sum(capacity.values())

    if requested_added_edges > maximum_added_edges:
        maximum_rate = maximum_added_edges / n_baseline_edges
        raise ValueError(
            f"Requested addition is impossible for mode {mode!r}. "
            f"Requested {requested_added_edges} new edges, but only "
            f"{maximum_added_edges} eligible non-edges exist. "
            f"Maximum addition relative to baseline: {maximum_rate:.6f}."
        )

    assigned_hosts = {
        v: set(hosts)
        for v, hosts in positive_hosts.items()
    }
    added_count_by_votu = defaultdict(int)
    eligible = DynamicVotuPool(
        [v for v in all_votus if capacity[v] > 0]
    )

    added_edges = set()

    while len(added_edges) < requested_added_edges:
        if len(eligible) == 0:
            raise RuntimeError(
                "No eligible vOTUs remain before the requested number of "
                "new edges was reached."
            )

        v = eligible.choice(rng)
        candidate = index.choose_candidate(
            mode=mode,
            original_hosts=positive_hosts[v],
            assigned_hosts=assigned_hosts[v],
            rng=rng,
            max_retry=max_retry,
        )

        if candidate is None:
            # Exact candidate selection found no remaining host, so this vOTU
            # is exhausted even if the count estimate was conservative.
            eligible.remove(v)
            continue

        edge = (candidate, v)
        if edge in baseline_edges or edge in added_edges:
            raise RuntimeError("Internal error: attempted to add a duplicate edge.")

        added_edges.add(edge)
        assigned_hosts[v].add(candidate)
        added_count_by_votu[v] += 1

        if added_count_by_votu[v] >= capacity[v]:
            eligible.remove(v)

    final_edges = baseline_edges | added_edges

    if len(final_edges) != n_baseline_edges + requested_added_edges:
        raise RuntimeError("Internal error: final edge count is incorrect.")

    output_parent = Path(fout).resolve().parent
    output_parent.mkdir(parents=True, exist_ok=True)

    with open(fout, "w", encoding="utf-8") as output:
        output.write("pOTU\tvOTU\n")
        for p, v in sorted(final_edges):
            output.write(f"{p}\t{v}\n")

    realized_addition = len(added_edges) / n_baseline_edges
    v_with_added_hosts = sum(count > 0 for count in added_count_by_votu.values())
    newly_multi_host_votus = sum(
        len(positive_hosts[v]) == 1 and len(assigned_hosts[v]) >= 2
        for v in all_votus
    )

    baseline_cross_genus = count_cross_rank_votus(
        positive_hosts, taxonomy, "Genus"
    )
    final_cross_genus = count_cross_rank_votus(
        assigned_hosts, taxonomy, "Genus"
    )
    baseline_cross_family = count_cross_rank_votus(
        positive_hosts, taxonomy, "Family"
    )
    final_cross_family = count_cross_rank_votus(
        assigned_hosts, taxonomy, "Family"
    )
    baseline_cross_order = count_cross_rank_votus(
        positive_hosts, taxonomy, "Order"
    )
    final_cross_order = count_cross_rank_votus(
        assigned_hosts, taxonomy, "Order"
    )

    taxonomic_summary = classify_added_edges(
        added_edges, positive_hosts, taxonomy
    )

    metrics_parent = Path(fm).resolve().parent
    metrics_parent.mkdir(parents=True, exist_ok=True)

    with open(fm, "w", encoding="utf-8") as metrics_output:
        metrics_output.write(
            "fout\tmode\trequested_addition\trealized_addition\t"
            "baseline_edges\trequested_added_edges\trealized_added_edges\t"
            "final_edges\tvOTUs_with_added_hosts\tnewly_multi_host_vOTUs\t"
            "baseline_cross_genus_vOTUs\tfinal_cross_genus_vOTUs\t"
            "baseline_cross_family_vOTUs\tfinal_cross_family_vOTUs\t"
            "baseline_cross_order_vOTUs\tfinal_cross_order_vOTUs\t"
            "same_genus_added_edges\tcross_genus_added_edges\t"
            "cross_family_added_edges\tcross_order_added_edges\t"
            "unknown_taxonomy_added_edges\n"
        )
        metrics_output.write(
            f"{Path(fout).name}\t{mode}\t{addition:.6f}\t{realized_addition:.6f}\t"
            f"{n_baseline_edges}\t{requested_added_edges}\t{len(added_edges)}\t"
            f"{len(final_edges)}\t{v_with_added_hosts}\t"
            f"{newly_multi_host_votus}\t"
            f"{baseline_cross_genus}\t{final_cross_genus}\t"
            f"{baseline_cross_family}\t{final_cross_family}\t"
            f"{baseline_cross_order}\t{final_cross_order}\t"
            f"{taxonomic_summary['same_genus_added_edges']}\t"
            f"{taxonomic_summary['cross_genus_added_edges']}\t"
            f"{taxonomic_summary['cross_family_added_edges']}\t"
            f"{taxonomic_summary['cross_order_added_edges']}\t"
            f"{taxonomic_summary['unknown_taxonomy_added_edges']}\n"
        )

    print("Missing-edge network generated.")
    print(f"Mode: {mode} ({MODE_LABELS[mode]})")
    print(f"Input file: {fbaseline}")
    print(f"Output file: {fout}")
    print(f"Metrics file: {fm}")
    print(f"Seed: {seed}")
    print(f"Baseline edges: {n_baseline_edges}")
    print(f"Requested added edges: {requested_added_edges}")
    print(f"Realized added edges: {len(added_edges)}")
    print(f"Realized addition: {realized_addition:.6f}")
    print(f"Final edges: {len(final_edges)}")
    print(f"vOTUs receiving new hosts: {v_with_added_hosts}")
    print(f"Newly multi-host vOTUs: {newly_multi_host_votus}")


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Generate repeated OVPIN networks by adding taxonomy-constrained "
            "missing virus-host edges."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--addition",
        type=float,
        required=True,
        help=(
            "Number of added edges as a percentage of the original edge count. "
            "For example, 20 adds 20%% as many new edges as the baseline has."
        ),
    )
    parser.add_argument(
        "--fbaseline",
        required=True,
        help="Baseline OVPIN file containing pOTU and vOTU columns.",
    )
    parser.add_argument(
        "--ftaxa",
        required=True,
        help=(
            "pOTU taxonomy file containing pOTU, Order, Family and Genus columns."
        ),
    )
    parser.add_argument(
        "--outfile",
        required=True,
        help=(
            "Output prefix. Three folders are generated: "
            "same-genus/different-pOTU, same-family/different-genus, "
            "and same-order/different-family."
        ),
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=100,
        help="Number of networks generated for each addition mode.",
    )
    parser.add_argument(
        "--seed-base",
        type=int,
        default=1000,
        help="Base random seed.",
    )
    parser.add_argument(
        "--max-retry",
        type=int,
        default=500,
        help="Fast candidate-sampling attempts before exact fallback.",
    )

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    if args.addition < 0:
        parser.error("--addition must be greater than or equal to 0.")
    if args.repeats < 1:
        parser.error("--repeats must be at least 1.")
    if args.max_retry < 1:
        parser.error("--max-retry must be at least 1.")
    if not Path(args.fbaseline).is_file():
        parser.error(f"Baseline input file not found: {args.fbaseline!r}.")
    if not Path(args.ftaxa).is_file():
        parser.error(f"Taxonomy input file not found: {args.ftaxa!r}.")

    addition = args.addition / 100.0

    # Read the large input files only once, then reuse them for every mode
    # and replicate.
    positive_hosts, baseline_edges = read_baseline_network(args.fbaseline)
    all_votus = sorted(positive_hosts)
    all_potus = sorted({p for p, _ in baseline_edges})
    taxonomy = read_potu_taxonomy(args.ftaxa, all_potus)
    index = TaxonomyIndex(all_potus, taxonomy)
    prepared = {
        "positive_hosts": positive_hosts,
        "baseline_edges": baseline_edges,
        "all_votus": all_votus,
        "all_potus": all_potus,
        "taxonomy": taxonomy,
        "index": index,
        "capacity_by_mode": {},
    }

    output_prefix = Path(args.outfile)
    parent = output_prefix.parent
    stem = output_prefix.name

    for mode_index, mode in enumerate(ALL_MODES):
        mode_dir = parent / f"{stem}_{mode}"
        mode_dir.mkdir(parents=True, exist_ok=True)
        file_prefix = mode_dir.name

        metrics_dir = parent / f"{stem}_{mode}_bimat" / "real_add_edges"
        metrics_dir.mkdir(parents=True, exist_ok=True)

        for i in range(1, args.repeats + 1):
            network_file = mode_dir / f"{file_prefix}_rep{i:03d}.txt"
            metrics_file = metrics_dir / f"{file_prefix}_rep{i:03d}.txt"
            seed = args.seed_base + mode_index * 1000000 + i

            add_missing_edges_by_taxonomy(
                fbaseline=args.fbaseline,
                ftaxa=args.ftaxa,
                fout=str(network_file),
                fm=str(metrics_file),
                addition=addition,
                mode=mode,
                seed=seed,
                max_retry=args.max_retry,
                prepared=prepared,
            )


if __name__ == "__main__":
    main()
