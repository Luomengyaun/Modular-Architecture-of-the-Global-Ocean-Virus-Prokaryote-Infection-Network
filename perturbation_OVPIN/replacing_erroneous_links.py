#!/usr/bin/env python3

import argparse
import math
import random
from collections import defaultdict
from pathlib import Path


MODE_ORDER = [
    "same_genus_diff_potu",
    "same_family_diff_genus",
    "same_order_diff_family",
    "cross_order",
]

MODE_LABELS = {
    "same_genus_diff_potu": "same genus, different pOTU",
    "same_family_diff_genus": "same family, different genus",
    "same_order_diff_family": "same order, different family",
    "cross_order": "different order",
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


def read_baseline_network(filename):
    positive_hosts = defaultdict(set)
    baseline_edges = set()

    with open(filename, "r", encoding="utf-8") as handle:
        first_line = handle.readline()
        if not first_line:
            raise ValueError("The baseline network file is empty.")

        header = first_line.rstrip("\n").split("\t")
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

            fields = line.rstrip("\n").split("\t")
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
    required_columns = ["pOTU", "Order", "Family", "Genus"]
    taxonomy = {}

    with open(filename, "r", encoding="utf-8") as handle:
        first_line = handle.readline()
        if not first_line:
            raise ValueError("The taxonomy file is empty.")

        header = first_line.rstrip("\n").split("\t")
        missing_columns = [column for column in required_columns if column not in header]
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

            fields = line.rstrip("\n").split("\t")
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
            f"Warning: {len(missing_potus)} baseline pOTUs have no taxonomy record "
            f"and therefore cannot be used as constrained replacement hosts: "
            f"{preview}{suffix}"
        )

    return taxonomy


class CandidateProvider:
    """Provide taxonomy-constrained replacement hosts for one perturbation mode."""

    def __init__(self, mode, all_potus, taxonomy):
        if mode not in MODE_ORDER:
            raise ValueError(f"Unknown perturbation mode: {mode}")

        self.mode = mode
        self.taxonomy = taxonomy
        self.all_potus = tuple(sorted(all_potus))

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
            if family is not None:
                family_groups[family].append(p)
            if order is not None:
                order_groups[order].append(p)
                ordered_potus.append(p)

        self.genus_groups = {key: tuple(value) for key, value in genus_groups.items()}
        self.family_groups = {key: tuple(value) for key, value in family_groups.items()}
        self.order_groups = {key: tuple(value) for key, value in order_groups.items()}
        self.ordered_potus = tuple(ordered_potus)

    def base_pool(self, original_host):
        record = self.taxonomy.get(original_host)
        if record is None:
            return ()

        if self.mode == "same_genus_diff_potu":
            genus = record["Genus"]
            if genus is None:
                return ()
            return self.genus_groups.get(genus, ())

        if self.mode == "same_family_diff_genus":
            family = record["Family"]
            genus = record["Genus"]
            if family is None or genus is None:
                return ()
            return self.family_groups.get(family, ())

        if self.mode == "same_order_diff_family":
            order = record["Order"]
            family = record["Family"]
            if order is None or family is None:
                return ()
            return self.order_groups.get(order, ())

        order = record["Order"]
        if order is None:
            return ()
        return self.ordered_potus

    def relation_is_valid(self, original_host, candidate_host):
        original = self.taxonomy.get(original_host)
        candidate = self.taxonomy.get(candidate_host)
        if original is None or candidate is None:
            return False

        if self.mode == "same_genus_diff_potu":
            return (
                original["Genus"] is not None
                and candidate["Genus"] == original["Genus"]
                and candidate_host != original_host
            )

        if self.mode == "same_family_diff_genus":
            return (
                original["Family"] is not None
                and original["Genus"] is not None
                and candidate["Family"] == original["Family"]
                and candidate["Genus"] is not None
                and candidate["Genus"] != original["Genus"]
            )

        if self.mode == "same_order_diff_family":
            return (
                original["Order"] is not None
                and original["Family"] is not None
                and candidate["Order"] == original["Order"]
                and candidate["Family"] is not None
                and candidate["Family"] != original["Family"]
            )

        return (
            original["Order"] is not None
            and candidate["Order"] is not None
            and candidate["Order"] != original["Order"]
        )

    @staticmethod
    def _rotated_indexes(length, rng):
        if length == 0:
            return

        if rng is None:
            for index in range(length):
                yield index
            return

        start = rng.randrange(length)
        direction = 1 if rng.random() < 0.5 else -1
        for offset in range(length):
            yield (start + direction * offset) % length

    def iter_candidates(self, original_host, positive_hosts, rng=None):
        pool = self.base_pool(original_host)

        for index in self._rotated_indexes(len(pool), rng):
            candidate = pool[index]
            if candidate in positive_hosts:
                continue
            if self.relation_is_valid(original_host, candidate):
                yield candidate

    def candidate_pool_size_upper_bound(self, original_host):
        return len(self.base_pool(original_host))


def find_replacement_matching(positive_hosts, provider, target_size=None, rng=None):
    """
    Match original positive hosts to distinct constrained non-positive hosts.

    The returned dictionary maps each replaced original host to its replacement.
    When target_size is None, a maximum matching is returned.
    """
    positives = set(positive_hosts)
    # Stabilize the starting order before any seeded shuffle.
    originals = sorted(positives)

    if rng is None:
        originals.sort(
            key=lambda p: (
                provider.candidate_pool_size_upper_bound(p),
                p,
            )
        )
    else:
        rng.shuffle(originals)

    candidate_to_original = {}

    def augment(original_host, seen_candidates):
        for candidate in provider.iter_candidates(original_host, positives, rng):
            if candidate in seen_candidates:
                continue

            seen_candidates.add(candidate)
            previous_original = candidate_to_original.get(candidate)

            if previous_original is None or augment(previous_original, seen_candidates):
                candidate_to_original[candidate] = original_host
                return True

        return False

    matched_count = 0
    for original_host in originals:
        if augment(original_host, set()):
            matched_count += 1
            if target_size is not None and matched_count >= target_size:
                break

    original_to_candidate = {
        original: candidate
        for candidate, original in candidate_to_original.items()
    }

    if target_size is not None and len(original_to_candidate) > target_size:
        # This should not normally occur, but keep exactly the requested number.
        items = sorted(original_to_candidate.items())
        if rng is not None:
            rng.shuffle(items)
        original_to_candidate = dict(items[:target_size])

    return original_to_candidate


def calculate_mode_min_retained(positive_hosts, provider):
    min_retained = {}

    for v, hosts in positive_hosts.items():
        matching = find_replacement_matching(hosts, provider, target_size=None, rng=None)
        max_replaceable = len(matching)
        min_retained[v] = len(hosts) - max_replaceable

    return min_retained


def calibrate_retained_counts_accuracy_only(
    all_votus,
    degree,
    min_retained,
    accuracy,
    rng,
):
    """Choose which vOTUs remain accurate when only Accuracy is controlled.

    Accurate vOTUs retain all original positive hosts.
    Inaccurate vOTUs retain none of their original positive hosts.
    Therefore Consistency is not an input constraint; it is a realized output metric.
    """
    n_votus = len(all_votus)
    target_accurate_count = round_half_up(accuracy * n_votus)
    realized_target_accuracy = target_accurate_count / n_votus

    # Any vOTU that cannot have all original hosts replaced in at least one mode
    # must remain accurate so that all four perturbation modes can share the same
    # accurate/inaccurate vOTU assignment.
    forced_accurate = {v for v in all_votus if min_retained[v] > 0}

    if target_accurate_count < len(forced_accurate):
        minimum_accuracy = len(forced_accurate) / n_votus
        raise ValueError(
            "Requested accuracy is impossible under the four shared taxonomy "
            "constraints. "
            f"At least {len(forced_accurate)} of {n_votus} vOTUs must retain "
            "their original hosts. "
            f"Minimum achievable accuracy: {minimum_accuracy:.6f}."
        )

    optional_votus = sorted(v for v in all_votus if v not in forced_accurate)
    n_optional_accurate = target_accurate_count - len(forced_accurate)

    rng.shuffle(optional_votus)
    selected_optional = set(optional_votus[:n_optional_accurate])
    accurate_votus = set(forced_accurate) | selected_optional

    retained_count = {
        v: degree[v] if v in accurate_votus else 0
        for v in all_votus
    }

    final_accurate_count = sum(retained_count[v] >= 1 for v in all_votus)
    final_accuracy = final_accurate_count / n_votus
    final_consistency = sum(
        retained_count[v] / degree[v]
        for v in all_votus
    ) / n_votus

    if final_accurate_count != target_accurate_count:
        raise RuntimeError(
            "Internal calibration error: the accurate-vOTU count does not "
            "match the target."
        )

    return retained_count, final_accuracy, final_consistency, final_accurate_count


def construct_perturbed_network(
    positive_hosts,
    baseline_edges,
    all_votus,
    degree,
    retained_count,
    provider,
    rng,
):
    assigned_hosts = {}

    for v in all_votus:
        positives = positive_hosts[v]
        n_retained = retained_count[v]
        n_incorrect = degree[v] - n_retained

        if n_incorrect == 0:
            assigned_hosts[v] = set(positives)
            continue

        matching = find_replacement_matching(
            positives,
            provider,
            target_size=n_incorrect,
            rng=rng,
        )

        if len(matching) != n_incorrect:
            raise RuntimeError(
                f"Could not generate {n_incorrect} constrained replacement hosts "
                f"for vOTU {v!r} in mode {provider.mode!r}."
            )

        removed_original_hosts = set(matching)
        retained_original_hosts = set(positives) - removed_original_hosts
        replacement_hosts = set(matching.values())
        assigned = retained_original_hosts | replacement_hosts

        if len(assigned) != degree[v]:
            raise RuntimeError(
                f"Internal error: degree of vOTU {v!r} changed in mode "
                f"{provider.mode!r}."
            )

        assigned_hosts[v] = assigned

    final_edges = {
        (p, v)
        for v, hosts in assigned_hosts.items()
        for p in hosts
    }

    if len(final_edges) != len(baseline_edges):
        raise RuntimeError(
            f"Internal error: total edge count changed in mode {provider.mode!r}."
        )

    return assigned_hosts, final_edges


def calculate_realized_metrics(
    assigned_hosts,
    positive_hosts,
    degree,
    all_votus,
):
    retained_positive_number = {
        v: len(assigned_hosts[v] & positive_hosts[v])
        for v in all_votus
    }

    final_accurate_count = sum(
        retained_positive_number[v] >= 1
        for v in all_votus
    )
    final_accuracy = final_accurate_count / len(all_votus)
    final_consistency = sum(
        retained_positive_number[v] / degree[v]
        for v in all_votus
    ) / len(all_votus)

    return final_accuracy, final_consistency, final_accurate_count


def write_network(filename, final_edges):
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", encoding="utf-8") as output:
        output.write("pOTU\tvOTU\n")
        for p, v in sorted(final_edges):
            output.write(f"{p}\t{v}\n")


def write_metrics_file(
    path,
    replicate,
    fout,
    mode,
    requested_accuracy,
    realized_accuracy,
    accurate_count,
    n_votus,
    realized_consistency,
    removed_count,
    added_count,
):
    """Write realized metrics, recording the network filename without its directory."""
    path.parent.mkdir(parents=True, exist_ok=True)

    with open(path, "w", encoding="utf-8") as output:
        output.write(
            "replicate\tfout\tmode\trequested_accuracy\trealized_accuracy\t"
            "accurate_vOTUs\trealized_consistency\t"
            "removed_original_edges\t"
            "added_replacement_edges\n"
        )
        output.write(
            f"{replicate}\t{Path(fout).name}\t{mode}\t"
            f"{requested_accuracy:.6f}\t{realized_accuracy:.6f}\t"
            f"{accurate_count}/{n_votus}\t"
            f"{realized_consistency:.6f}\t"
            f"{removed_count}\t{added_count}\n"
        )


def build_parser():
    parser = argparse.ArgumentParser(
        description=(
            "Generate four sets of accuracy-only taxonomy-constrained perturbed OVPIN networks: "
            "same genus/different pOTU, same family/different genus, "
            "same order/different family, and cross-order."
        ),
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    parser.add_argument(
        "--accuracy",
        type=float,
        required=True,
        help="Target final accuracy as a percentage, for example 80.",
    )
    parser.add_argument(
        "--fbaseline",
        required=True,
        help="Baseline network file containing pOTU and vOTU columns.",
    )
    parser.add_argument(
        "--ftaxa",
        required=True,
        help=(
            "pOTU taxonomy file containing pOTU, Order, Family, and Genus columns."
        ),
    )
    parser.add_argument(
        "--outfile",
        required=True,
        help=(
            "Base output prefix. Four directories are created by appending the "
            "four perturbation mode names."
        ),
    )
    parser.add_argument(
        "--repeats",
        type=int,
        default=100,
        help="Number of replicate networks generated for each perturbation mode.",
    )
    parser.add_argument(
        "--seed-base",
        type=int,
        default=1000,
        help="Base random seed; replicate i uses seed-base + i.",
    )
    parser.add_argument(
        "--max-retry",
        type=int,
        default=100,
        help="Retained for command-line compatibility; unused in accuracy-only mode.",
    )

    return parser


def main():
    parser = build_parser()
    args = parser.parse_args()

    if not 0 <= args.accuracy <= 100:
        parser.error("--accuracy must be between 0 and 100.")
    if args.repeats < 1:
        parser.error("--repeats must be at least 1.")
    if args.max_retry < 1:
        parser.error("--max-retry must be at least 1.")
    if not Path(args.fbaseline).is_file():
        parser.error(f"Baseline input file not found: {args.fbaseline!r}.")
    if not Path(args.ftaxa).is_file():
        parser.error(f"Taxonomy input file not found: {args.ftaxa!r}.")

    accuracy = args.accuracy / 100.0

    positive_hosts, baseline_edges = read_baseline_network(args.fbaseline)
    all_votus = sorted(positive_hosts)
    all_potus = sorted({p for p, _ in baseline_edges})
    degree = {v: len(positive_hosts[v]) for v in all_votus}

    taxonomy = read_potu_taxonomy(args.ftaxa, all_potus)

    providers = {
        mode: CandidateProvider(mode, all_potus, taxonomy)
        for mode in MODE_ORDER
    }

    print("Precomputing the maximum replaceable-host capacity for each mode...")
    mode_min_retained = {}
    for mode in MODE_ORDER:
        mode_min_retained[mode] = calculate_mode_min_retained(
            positive_hosts,
            providers[mode],
        )
        forced = sum(value > 0 for value in mode_min_retained[mode].values())
        print(
            f"  {mode}: {forced}/{len(all_votus)} vOTUs cannot be made fully "
            "inaccurate under this constraint."
        )

    # Use a shared per-vOTU lower bound so all four modes have the same
    # retained-positive count for every vOTU in a given replicate.
    shared_min_retained = {
        v: max(mode_min_retained[mode][v] for mode in MODE_ORDER)
        for v in all_votus
    }

    output_dirs = {
        mode: Path(f"{args.outfile}_{mode}")
        for mode in MODE_ORDER
    }
    metrics_dirs = {}

    for mode, directory in output_dirs.items():
        directory.mkdir(parents=True, exist_ok=True)

        # Match the directory layout used by the original random-perturbation script:
        # <mode_output>_bimat/real_acc_con/<mode_output>_repXXX.txt
        metrics_directory = Path(f"{directory}_bimat") / "real_acc_con"
        metrics_directory.mkdir(parents=True, exist_ok=True)
        metrics_dirs[mode] = metrics_directory

    print()
    print(f"Baseline vOTUs: {len(all_votus)}")
    print(f"Baseline pOTUs: {len(all_potus)}")
    print(f"Baseline edges: {len(baseline_edges)}")
    print(f"Requested final accuracy: {accuracy:.6f}")
    print("Mode: accuracy-only; consistency is measured, not constrained.")
    print(f"Replicates per mode: {args.repeats}")
    print()

    for replicate in range(1, args.repeats + 1):
        seed = args.seed_base + replicate
        calibration_rng = random.Random(seed)

        retained_count, calibrated_accuracy, calibrated_consistency, calibrated_count = (
            calibrate_retained_counts_accuracy_only(
                all_votus=all_votus,
                degree=degree,
                min_retained=shared_min_retained,
                accuracy=accuracy,
                rng=calibration_rng,
            )
        )

        for mode_index, mode in enumerate(MODE_ORDER):
            mode_seed = seed * 100 + mode_index
            mode_rng = random.Random(mode_seed)

            assigned_hosts, final_edges = construct_perturbed_network(
                positive_hosts=positive_hosts,
                baseline_edges=baseline_edges,
                all_votus=all_votus,
                degree=degree,
                retained_count=retained_count,
                provider=providers[mode],
                rng=mode_rng,
            )

            final_accuracy, final_consistency, final_accurate_count = (
                calculate_realized_metrics(
                    assigned_hosts=assigned_hosts,
                    positive_hosts=positive_hosts,
                    degree=degree,
                    all_votus=all_votus,
                )
            )

            if not math.isclose(
                final_accuracy,
                calibrated_accuracy,
                rel_tol=0,
                abs_tol=1e-12,
            ):
                raise RuntimeError(
                    f"Mode {mode!r} did not preserve the calibrated accuracy."
                )

            if not math.isclose(
                final_consistency,
                calibrated_consistency,
                rel_tol=0,
                abs_tol=1e-12,
            ):
                raise RuntimeError(
                    f"Mode {mode!r} did not preserve the calibrated consistency."
                )

            directory = output_dirs[mode]
            prefix = directory.name
            fout = directory / f"{prefix}_rep{replicate:03d}.txt"
            write_network(fout, final_edges)

            removed_edges = baseline_edges - final_edges
            added_edges = final_edges - baseline_edges

            fm = (
                metrics_dirs[mode]
                / f"{prefix}_rep{replicate:03d}.txt"
            )

            write_metrics_file(
                path=fm,
                replicate=replicate,
                fout=str(fout),
                mode=mode,
                requested_accuracy=accuracy,
                realized_accuracy=final_accuracy,
                accurate_count=final_accurate_count,
                n_votus=len(all_votus),
                realized_consistency=final_consistency,
                removed_count=len(removed_edges),
                added_count=len(added_edges),
            )

        print(
            f"Completed replicate {replicate:03d}/{args.repeats}: "
            f"accuracy={calibrated_accuracy:.6f}, "
            f"consistency={calibrated_consistency:.6f}"
        )

    print()
    print("Finished. Four network directories were generated:")
    for mode in MODE_ORDER:
        print(f"  {output_dirs[mode]}  ({MODE_LABELS[mode]})")


if __name__ == "__main__":
    main()
