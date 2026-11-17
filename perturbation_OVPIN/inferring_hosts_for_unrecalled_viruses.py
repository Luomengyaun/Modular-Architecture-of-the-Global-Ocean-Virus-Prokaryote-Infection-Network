#!/usr/bin/env python3
"""Expand OVPIN with new OVG vOTUs and OPG hosts in three paired scenarios.

Python >= 3.8; standard library only. See README.md for the
scientific assumptions, input formats, output files and command examples.
"""

import argparse
import ast
import bisect
import csv
import hashlib
import heapq
import json
import math
import random
import sys
from collections import Counter, defaultdict, namedtuple
from pathlib import Path


MODES = ("same_genus", "same_family", "same_order")
MISSING = {"", "na", "n/a", "none", "null", "nan", "unknown", "unclassified"}
Taxon = namedtuple("Taxon", "order family genus")
Candidate = namedtuple(
    "Candidate", "votu group source_bc source_direct source_group genera"
)

def normalize_taxon(value, rank=None):
    """Remove a rank prefix, but preserve meaningful GTDB suffixes such as _A."""
    value = value.strip()
    if rank and value.startswith(rank + "__"):
        value = value[3:].strip()
    return None if value.lower() in MISSING else value


def data_lines(filename):
    with open(filename, "r", encoding="utf-8-sig") as handle:
        for number, line in enumerate(handle, 1):
            if line.strip():
                yield number, line.rstrip("\r\n").split("\t")


def fail_row(filename, number, message):
    raise ValueError("{}:{}: {}".format(filename, number, message))


def read_baseline(filename):
    """Use only the edge columns as truth; embedded taxonomy is diagnostic."""
    edges = set()
    hosts = defaultdict(set)
    embedded_groups = defaultdict(set)
    stats = Counter()
    indexes = (0, 1)
    first = True
    for number, fields in data_lines(filename):
        header = [s.strip().lower() for s in fields]
        if first and "potu" in header and "votu" in header:
            indexes = (header.index("potu"), header.index("votu"))
            first = False
            continue
        first = False
        if len(fields) <= max(indexes):
            fail_row(filename, number, "expected pOTU and vOTU columns")
        p, v = (fields[i].strip() for i in indexes)
        if not p or not v:
            fail_row(filename, number, "empty pOTU or vOTU")
        stats["data_rows"] += 1
        if (p, v) in edges:
            stats["duplicate_edges"] += 1
        edges.add((p, v))
        hosts[v].add(p)
        # A four-column baseline may contain an old [fvOTU, gvOTU, vOTU]
        # list in column 4. Never use it as a fallback for the current map.
        if indexes == (0, 1) and len(fields) >= 4:
            try:
                value = ast.literal_eval(fields[3])
                if isinstance(value, (list, tuple)) and len(value) >= 2:
                    if isinstance(value[1], str):
                        group = normalize_taxon(value[1])
                        if group:
                            embedded_groups[v].add(group)
                else:
                    stats["unreadable_embedded_virus_taxonomy_rows"] += 1
            except (ValueError, SyntaxError, TypeError, RecursionError):
                stats["unreadable_embedded_virus_taxonomy_rows"] += 1
    if not edges:
        raise ValueError("The baseline network has no valid edges.")
    return dict(hosts), edges, dict(embedded_groups), dict(stats)


def read_opg_taxonomy(filename):
    """Read the full OPG pOTU table, not just the baseline's host taxonomy."""
    records = {}
    stats = Counter()
    indexes = {"potu": 0, "order": 5, "family": 6, "genus": 7}
    first = True
    for number, fields in data_lines(filename):
        header = [s.strip().lower() for s in fields]
        if first and "potu" in header:
            missing = set(indexes) - set(header)
            if missing:
                fail_row(filename, number, "missing headers: " + ", ".join(sorted(missing)))
            indexes = {name: header.index(name) for name in indexes}
            first = False
            continue
        first = False
        if len(fields) <= max(indexes.values()):
            fail_row(filename, number, "expected the OPG 9-column table or named columns")
        p = fields[indexes["potu"]].strip()
        if not p:
            fail_row(filename, number, "empty pOTU")
        record = Taxon(
            normalize_taxon(fields[indexes["order"]], "o"),
            normalize_taxon(fields[indexes["family"]], "f"),
            normalize_taxon(fields[indexes["genus"]], "g"),
        )
        stats["data_rows"] += 1
        if p in records:
            if records[p] != record:
                fail_row(filename, number, "conflicting taxonomy for pOTU " + p)
            stats["duplicate_identical_potu_rows"] += 1
        records[p] = record
    if not records:
        raise ValueError("The OPG taxonomy file has no pOTUs.")
    for rank in ("genus", "family", "order"):
        stats["potus_with_unknown_" + rank] = sum(
            getattr(r, rank) is None for r in records.values()
        )
    return records, dict(stats)


def read_votu_groups(filename):
    result = {}
    stats = Counter()
    first = True
    for number, fields in data_lines(filename):
        if first and [s.strip().lower() for s in fields[:2]] == ["votu", "gvotu"]:
            first = False
            continue
        first = False
        if len(fields) < 2:
            fail_row(filename, number, "expected vOTU and gvOTU")
        v = fields[0].strip()
        group = normalize_taxon(fields[1])
        if not v:
            fail_row(filename, number, "empty vOTU")
        if group is not None:
            group = sys.intern(group)
        stats["data_rows"] += 1
        if v in result:
            if result[v] != group:
                fail_row(filename, number, "conflicting gvOTU assignments for " + v)
            stats["duplicate_identical_assignments"] += 1
        result[v] = group
    if not result:
        raise ValueError("The vOTU-gvOTU map is empty.")
    stats["votus_with_unknown_gvotu"] = sum(g is None for g in result.values())
    return result, dict(stats)


def read_iphop(filename, known_genera):
    """Safely parse Python dict literals; filter BEFORE donor propagation."""
    original = defaultdict(set)
    row_counts = Counter()
    stats = Counter()
    first = True
    for number, fields in data_lines(filename):
        if first and fields[0].strip().lower() == "votu":
            first = False
            continue
        first = False
        if len(fields) != 2:
            fail_row(filename, number, "expected vOTU<TAB>{host_genus: value, ...}")
        v = fields[0].strip()
        if not v:
            fail_row(filename, number, "empty vOTU")
        try:
            value = ast.literal_eval(fields[1])
        except (ValueError, SyntaxError, TypeError, RecursionError) as exc:
            fail_row(filename, number, "invalid dictionary literal: " + str(exc))
        if not isinstance(value, dict) or not all(isinstance(g, str) for g in value):
            fail_row(filename, number, "host predictions must be a dictionary with string keys")
        stats["data_rows"] += 1
        row_counts[v] += 1
        original[v]  # Preserve records with an empty prediction dictionary.
        for g in value:
            g = normalize_taxon(g, "g")
            if g is not None:
                original[v].add(g)
            else:
                stats["unknown_genus_keys"] += 1
    retained = {}
    missing = Counter()
    for v, genera in original.items():
        keep = genera & known_genera
        for g in genera - known_genera:
            missing[g] += 1
        if keep:
            retained[v] = keep
            if keep != genera:
                stats["partly_filtered_votus"] += 1
        else:
            stats["votus_without_matched_genus"] += 1
    stats["input_unique_votus"] = len(original)
    stats["retained_votus_before_OVG_filter"] = len(retained)
    stats["duplicate_votu_rows_merged"] = sum(n - 1 for n in row_counts.values())
    stats["unmatched_unique_genera"] = len(missing)
    stats["input_votu_genus_pairs"] = sum(len(gs) for gs in original.values())
    stats["retained_votu_genus_pairs_before_OVG_filter"] = sum(map(len, retained.values()))
    return retained, missing, dict(stats)


class HostPool:
    """A union of shared, disjoint host groups; no expanded host list per virus."""

    __slots__ = ("groups", "ends", "size")

    def __init__(self, groups):
        self.groups = tuple(group for group in groups if group)
        total = 0
        ends = []
        for group in self.groups:
            total += len(group)
            ends.append(total)
        self.ends = tuple(ends)
        self.size = total

    def sample(self, degree, rng):
        if not 0 <= degree <= self.size:
            raise ValueError("Requested degree exceeds candidate-host capacity.")
        hosts = []
        for position in rng.sample(range(self.size), degree):
            group_index = bisect.bisect_right(self.ends, position)
            previous = self.ends[group_index - 1] if group_index else 0
            hosts.append(self.groups[group_index][position - previous])
        return hosts


class HostIndex:
    """OPG indexes with cumulative ranges, including lower-rank unknown hosts."""

    def __init__(self, taxonomy):
        self.taxonomy = taxonomy
        groups = {name: defaultdict(list) for name in (
            "genus", "family", "order", "genus_missing_family",
            "family_missing_order", "genus_missing_family_order",
        )}
        genus_families = defaultdict(set)
        genus_orders = defaultdict(set)
        family_orders = defaultdict(set)
        for p, t in taxonomy.items():
            for rank in ("genus", "family", "order"):
                value = getattr(t, rank)
                if value is not None:
                    groups[rank][value].append(p)
            if t.genus:
                if t.family:
                    genus_families[t.genus].add(t.family)
                else:
                    groups["genus_missing_family"][t.genus].append(p)
                if t.order:
                    genus_orders[t.genus].add(t.order)
                elif not t.family:
                    groups["genus_missing_family_order"][t.genus].append(p)
            if t.family:
                if t.order:
                    family_orders[t.family].add(t.order)
                else:
                    groups["family_missing_order"][t.family].append(p)
        for label, mapping in (
            ("genus-to-family", genus_families),
            ("genus-to-order", genus_orders),
            ("family-to-order", family_orders),
        ):
            conflicts = [key for key, values in mapping.items() if len(values) > 1]
            if conflicts:
                raise ValueError("Conflicting OPG {} taxonomy for: {}".format(
                    label, ", ".join(sorted(conflicts)[:10])
                ))
        self.groups = {
            name: {key: tuple(sorted(ps)) for key, ps in values.items()}
            for name, values in groups.items()
        }
        self.genus_families = genus_families
        self.genus_orders = genus_orders
        self.family_orders = family_orders
        self.cache = {}

    @property
    def known_genera(self):
        return set(self.groups["genus"])

    def reference_ranks(self, genera):
        families = set()
        orders = set()
        for g in genera:
            families.update(self.genus_families.get(g, ()))
            orders.update(self.genus_orders.get(g, ()))
        for f in families:
            orders.update(self.family_orders.get(f, ()))
        return families, orders

    def pool(self, genera, mode):
        key = (genera, mode)
        if key in self.cache:
            return self.cache[key]
        families, orders = self.reference_ranks(genera)
        if mode == "same_genus":
            parts = [self.groups["genus"][g] for g in genera]
        elif mode == "same_family":
            parts = [self.groups["family"][f] for f in sorted(families)]
            parts += [self.groups["genus_missing_family"].get(g, ()) for g in genera]
        elif mode == "same_order":
            parts = [self.groups["order"][o] for o in sorted(orders)]
            parts += [self.groups["family_missing_order"].get(f, ()) for f in sorted(families)]
            parts += [self.groups["genus_missing_family_order"].get(g, ()) for g in genera]
        else:
            raise ValueError("Unknown mode: " + mode)
        result = HostPool(parts)
        self.cache[key] = result
        return result


def build_library(baseline_hosts, taxonomy, votu_groups, iphop, index):
    baseline_votus = set(baseline_hosts)
    bc_groups = set()
    bc_genera = defaultdict(set)
    for v, hosts in baseline_hosts.items():
        group = votu_groups.get(v)
        if group is None:
            continue
        bc_groups.add(group)
        for p in hosts:
            genus = taxonomy[p].genus
            if genus:
                bc_genera[group].add(genus)
    group_donors = Counter()
    group_iphop_genera = defaultdict(set)
    for v, genera in iphop.items():
        group = votu_groups.get(v)
        if group is not None:
            group_donors[group] += 1
            group_iphop_genera[group].update(genera)

    signatures = {}
    records = []
    stats = Counter()
    for v, group in votu_groups.items():
        if v in baseline_votus:
            continue
        source_bc = group is not None and group in bc_groups
        source_direct = v in iphop
        # A virus is not its own "other vOTU" donor for source 3.
        other_donors = group_donors.get(group, 0) - int(source_direct and group is not None)
        source_group = group is not None and other_donors > 0
        if not (source_bc or source_direct or source_group):
            continue
        genera = set(bc_genera.get(group, ()))
        genera.update(iphop.get(v, ()))
        genera.update(group_iphop_genera.get(group, ()))
        signature = tuple(sorted(genera))
        signature = signatures.setdefault(signature, signature)
        record = Candidate(v, group, source_bc, source_direct, source_group, signature)
        records.append(record)
        stats["source_bc_gvotu"] += int(source_bc)
        stats["source_iphop_direct"] += int(source_direct)
        stats["source_iphop_gvotu"] += int(source_group)
        stats["multiple_evidence_sources"] += int(sum((source_bc, source_direct, source_group)) > 1)
    records.sort(key=lambda r: r.votu)
    eligible = [r for r in records if index.pool(r.genera, "same_genus").size > 0]
    stats["complete_new_virus_library"] = len(records)
    stats["eligible_for_paired_three_modes"] = len(eligible)
    stats["no_same_genus_candidate"] = len(records) - len(eligible)
    return records, eligible, dict(stats)


def target_degrees(baseline_hist, count, rng):
    """Largest-remainder proportional allocation, with seeded random tie breaks."""
    if count == 0:
        return []
    total = sum(baseline_hist.values())
    allocated = {}
    remainders = []
    for degree, frequency in sorted(baseline_hist.items()):
        quotient, remainder = divmod(count * frequency, total)
        allocated[degree] = quotient
        remainders.append((remainder, rng.random(), degree))
    remaining = count - sum(allocated.values())
    for _, _, degree in sorted(remainders, reverse=True)[:remaining]:
        allocated[degree] += 1
    return [degree for degree in sorted(allocated, reverse=True)
            for _ in range(allocated[degree])]


def assign_degrees(selected, targets, index, rng):
    """Match high degrees first, uniformly across capacity-eligible viruses.

    If no remaining virus can take a target degree, use the largest remaining
    capacity. This alters the quota only after reassignment cannot avoid it.
    """
    if len(selected) != len(targets):
        raise ValueError("Degree count differs from the selected-virus count.")
    waiting = sorted(
        ((index.pool(r.genera, "same_genus").size, r) for r in selected),
        key=lambda item: (-item[0], item[1].votu),
    )
    active = []
    position = 0
    assigned = {}
    for target in sorted(targets, reverse=True):
        while position < len(waiting) and waiting[position][0] >= target:
            active.append(waiting[position])
            position += 1
        if not active:
            if position >= len(waiting):
                raise RuntimeError("No viruses remain for degree assignment.")
            maximum = waiting[position][0]
            while position < len(waiting) and waiting[position][0] == maximum:
                active.append(waiting[position])
                position += 1
        pick = rng.randrange(len(active))
        capacity, record = active[pick]
        active[pick] = active[-1]
        active.pop()
        degree = min(target, capacity)
        if degree < 1:
            raise RuntimeError("A zero-capacity virus entered the selected set.")
        assigned[record.votu] = (target, degree, capacity)
    return assigned


def seed_for(base, repetition, purpose):
    value = "{}|{}|{}".format(base, repetition, purpose).encode("utf-8")
    return int.from_bytes(hashlib.sha256(value).digest()[:8], "big")


def mean_degree(hist):
    n = sum(hist.values())
    return sum(k * frequency for k, frequency in hist.items()) / n if n else 0.0


def distribution_distance(first, second):
    """Total variation and maximum CDF difference, not significance tests."""
    n1, n2 = sum(first.values()), sum(second.values())
    if not n1 or not n2:
        return "", ""
    tv = 0.0
    cdf1 = cdf2 = ks = 0.0
    for degree in sorted(set(first) | set(second)):
        p, q = first.get(degree, 0) / n1, second.get(degree, 0) / n2
        tv += abs(p - q)
        cdf1 += p
        cdf2 += q
        ks = max(ks, abs(cdf1 - cdf2))
    return tv / 2, ks


def write_json(path, value):
    with path.open("w", encoding="utf-8") as handle:
        json.dump(value, handle, ensure_ascii=False, indent=2)
        handle.write("\n")


def write_distribution(path, baseline_hist, target_hist=None, assigned_hist=None):
    target_hist = target_hist or Counter()
    assigned_hist = assigned_hist or Counter()
    histograms = (baseline_hist, target_hist, assigned_hist)
    totals = [sum(h.values()) for h in histograms]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("degree", "baseline_count", "baseline_fraction",
                         "target_count", "target_fraction", "assigned_count", "assigned_fraction"))
        for degree in sorted(set().union(*(set(h) for h in histograms))):
            row = [degree]
            for hist, total in zip(histograms, totals):
                number = hist.get(degree, 0)
                row.extend((number, number / total if total else 0.0))
            writer.writerow(row)


def write_audit(audit_dir, records, index, baseline_hist, iphop, unmatched_genera, summary):
    audit_dir.mkdir(parents=True)
    write_json(audit_dir / "input_summary.json", summary)
    write_distribution(audit_dir / "baseline_virus_degree_distribution.tsv", baseline_hist)
    with (audit_dir / "new_virus_library.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("vOTU", "gvOTU", "source_bc_gvotu", "source_iphop_direct",
                         "source_iphop_gvotu", "reference_genera", "same_genus_capacity",
                         "same_family_capacity", "same_order_capacity", "eligible", "status"))
        for r in records:
            capacities = [index.pool(r.genera, mode).size for mode in MODES]
            if not capacities[0] <= capacities[1] <= capacities[2]:
                raise RuntimeError("Candidate pools are not cumulative for " + r.votu)
            eligible = capacities[0] > 0
            writer.writerow((r.votu, r.group or "", int(r.source_bc), int(r.source_direct),
                             int(r.source_group), ";".join(r.genera), *capacities,
                             int(eligible), "eligible" if eligible else "no_known_reference_genus"))
    with (audit_dir / "iphop_filtered_for_OVG.txt").open("w", encoding="utf-8") as handle:
        for v, genera in sorted(iphop.items()):
            handle.write(v + "\t" + repr({g: "" for g in sorted(genera)}) + "\n")
    with (audit_dir / "iphop_unmatched_genera.tsv").open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        writer.writerow(("genus", "input_vOTU_count"))
        writer.writerows(sorted(unmatched_genera.items(), key=lambda item: (-item[1], item[0])))


def choose_count(args, baseline_v_count, eligible_count):
    if args.n_new_votus is not None:
        count = args.n_new_votus
    elif args.addition is not None:
        count = int(math.floor(args.addition * baseline_v_count / 100.0 + 0.5))
    else:
        count = eligible_count
    if count > eligible_count:
        raise ValueError(
            "Requested {} new vOTUs, but only {} have same-genus candidates. "
            "Lower --addition/--n-new-votus or use all eligible viruses by omitting both."
            .format(count, eligible_count)
        )
    return count


def run(args):
    for label, filename in (("baseline", args.fbaseline), ("OPG taxonomy", args.ftaxa),
                            ("vOTU-gvOTU map", args.fvotu_gvotu), ("iPHoP", args.fiphop)):
        if not Path(filename).is_file():
            raise ValueError("{} file not found: {}".format(label, filename))
    prefix = Path(args.outfile)
    stem = prefix.name
    if not stem or stem in (".", ".."):
        raise ValueError("--outfile must include an output prefix name.")
    parent = prefix.parent
    audit_dir = parent / (stem + "_audit")
    paired_dir = parent / (stem + "_paired")
    mode_dirs = {mode: parent / (stem + "_" + mode) for mode in MODES}
    metrics_path = parent / (stem + "_metrics.tsv")
    manifest_path = parent / (stem + "_run.json")
    for path in [audit_dir, paired_dir, metrics_path, manifest_path] + list(mode_dirs.values()):
        if path.exists():
            raise ValueError("Output already exists; choose a new --outfile prefix: " + str(path))

    print("Reading baseline and complete OPG taxonomy...", flush=True)
    baseline_hosts, baseline_edges, embedded, baseline_stats = read_baseline(args.fbaseline)
    taxonomy, taxonomy_stats = read_opg_taxonomy(args.ftaxa)
    baseline_votus = set(baseline_hosts)
    baseline_potus = {p for p, v in baseline_edges}
    absent_hosts = baseline_potus - set(taxonomy)
    if absent_hosts:
        raise ValueError("{} baseline pOTUs are absent from OPG, for example: {}".format(
            len(absent_hosts), ", ".join(sorted(absent_hosts)[:10])
        ))
    index = HostIndex(taxonomy)
    print("Reading current vOTU-gvOTU map and filtering iPHoP genera...", flush=True)
    votu_groups, group_stats = read_votu_groups(args.fvotu_gvotu)
    absent_votus = baseline_votus - set(votu_groups)
    if absent_votus and not args.allow_missing_baseline_groups:
        raise ValueError(
            "{} baseline vOTUs are absent from the current vOTU-gvOTU map, for example: {}. "
            "Use --allow-missing-baseline-groups only to explicitly retain these original "
            "nodes without using them as gvOTU donors."
            .format(len(absent_votus), ", ".join(sorted(absent_votus)[:10]))
        )
    embedded_mismatches = sum(
        v in votu_groups and (len(groups) != 1 or votu_groups[v] not in groups)
        for v, groups in embedded.items()
    )
    iphop_all, unmatched_genera, iphop_stats = read_iphop(args.fiphop, index.known_genera)
    iphop = {v: gs for v, gs in iphop_all.items() if v in votu_groups}
    iphop_stats["retained_votus_outside_OVG_map_excluded"] = len(iphop_all) - len(iphop)
    iphop_stats["retained_reference_votus_in_OVG"] = len(iphop)
    print("Building the new-virus library and nested host candidate pools...", flush=True)
    records, eligible, library_stats = build_library(
        baseline_hosts, taxonomy, votu_groups, iphop, index
    )
    n_selected = choose_count(args, len(baseline_votus), len(eligible))
    baseline_hist = Counter(len(ps) for ps in baseline_hosts.values())
    base_mean = mean_degree(baseline_hist)
    base_conn = len(baseline_edges) / (len(baseline_votus) * len(baseline_potus))
    summary = {
        "inputs": {"baseline": Path(args.fbaseline).name,
                   "OPG_taxonomy": Path(args.ftaxa).name,
                   "votu_gvotu": Path(args.fvotu_gvotu).name,
                   "iphop": Path(args.fiphop).name},
        "baseline": dict(baseline_stats, edges=len(baseline_edges), votus=len(baseline_votus),
                         potus=len(baseline_potus), mean_virus_degree=base_mean, connectance=base_conn),
        "OPG_taxonomy": dict(taxonomy_stats, unique_potus=len(taxonomy),
                             known_genera=len(index.known_genera)),
        "votu_gvotu": dict(group_stats, unique_OVG_votus=len(votu_groups),
                           baseline_votus_missing_map=len(absent_votus),
                           baseline_votus_with_unknown_group=sum(votu_groups.get(v) is None for v in baseline_votus if v in votu_groups),
                           embedded_baseline_group_mismatches=embedded_mismatches),
        "iphop": iphop_stats,
        "new_virus_library": library_stats,
        "new_votus_per_repetition": n_selected,
        "reference_rule": "Union of BC same-gvOTU host genera, direct filtered iPHoP genera, and same-gvOTU filtered iPHoP genera; OPG taxonomy is authoritative.",
        "candidate_rule": "Cumulative genus <= family <= order; unknown taxa are never shared labels; known narrower candidates are retained.",
        "degree_rule": "Largest-remainder baseline degree frequencies, then descending capacity-constrained randomized assignment; reduce only infeasible quotas and report deviations.",
        "connectance_rule": "Measured using actual connected nodes; not constrained.",
    }
    parent.mkdir(parents=True, exist_ok=True)
    recorded_args = dict(vars(args))
    for name in ("fbaseline", "ftaxa", "fvotu_gvotu", "fiphop", "outfile"):
        recorded_args[name] = Path(recorded_args[name]).name
    manifest = {"status": "preparing", "arguments": recorded_args, "modes": MODES,
                "summary": (Path(audit_dir.name) / "input_summary.json").as_posix()}
    write_json(manifest_path, manifest)
    write_audit(audit_dir, records, index, baseline_hist, iphop, unmatched_genera, summary)
    print("New-virus library: {}; eligible: {}; selected per repetition: {}.".format(
        len(records), len(eligible), n_selected
    ), flush=True)
    if absent_votus:
        print("Warning: {} baseline viruses cannot donate through the current group map.".format(len(absent_votus)), flush=True)
    if embedded_mismatches:
        print("Note: {} baseline embedded group labels differ; the supplied current map is used.".format(embedded_mismatches), flush=True)
    if args.prepare_only:
        manifest["status"] = "prepared_only"
        write_json(manifest_path, manifest)
        print("Preparation complete. Audit: " + str(audit_dir), flush=True)
        return

    # Release the full OVG map and propagation intermediates before replicates.
    del votu_groups, iphop_all, iphop, embedded, baseline_hosts, records
    paired_dir.mkdir()
    for directory in mode_dirs.values():
        directory.mkdir()
    sorted_baseline = sorted(baseline_edges)
    metrics_writer = None
    manifest["status"] = "running"
    write_json(manifest_path, manifest)
    with metrics_path.open("w", encoding="utf-8", newline="") as metrics_file:
        for repetition in range(1, args.repeats + 1):
            selection_seed = seed_for(args.seed_base, repetition, "selection")
            quota_seed = seed_for(args.seed_base, repetition, "quota")
            assignment_seed = seed_for(args.seed_base, repetition, "assignment")
            selection_rng = random.Random(selection_seed)
            selected = (list(eligible) if n_selected == len(eligible)
                        else selection_rng.sample(eligible, n_selected))
            selected.sort(key=lambda r: r.votu)
            targets = target_degrees(baseline_hist, n_selected, random.Random(quota_seed))
            assignments = assign_degrees(selected, targets, index, random.Random(assignment_seed))
            target_hist = Counter(targets)
            assigned_hist = Counter(value[1] for value in assignments.values())
            added_count = sum(k * n for k, n in assigned_hist.items())
            reduced = sum(target > degree for target, degree, capacity in assignments.values())
            lost = sum(target - degree for target, degree, capacity in assignments.values())
            target_tv, target_ks = distribution_distance(baseline_hist, target_hist)
            actual_tv, actual_ks = distribution_distance(baseline_hist, assigned_hist)
            capacity_tv, capacity_ks = distribution_distance(target_hist, assigned_hist)
            rep = "rep{:03d}".format(repetition)
            with (paired_dir / (rep + "_degree_quotas.tsv")).open("w", encoding="utf-8", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t")
                writer.writerow(("vOTU", "target_degree", "assigned_degree", "same_genus_capacity"))
                for r in selected:
                    writer.writerow((r.votu, *assignments[r.votu]))
            write_distribution(paired_dir / (rep + "_degree_distribution.tsv"),
                               baseline_hist, target_hist, assigned_hist)

            for mode in MODES:
                host_seed = seed_for(args.seed_base, repetition, "hosts:" + mode)
                rng = random.Random(host_seed)
                added = []
                added_hosts = set()
                for r in selected:
                    degree = assignments[r.votu][1]
                    hosts = index.pool(r.genera, mode).sample(degree, rng)
                    added.extend((p, r.votu) for p in hosts)
                    added_hosts.update(hosts)
                if len(added) != added_count:
                    raise RuntimeError("Added edge count differs from the paired degree quota.")
                added.sort()
                file_stem = "{}_{}_{}".format(stem, mode, rep)
                network_path = mode_dirs[mode] / (file_stem + ".txt")
                added_path = mode_dirs[mode] / (file_stem + "_added_edges.tsv")
                with network_path.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t")
                    writer.writerow(("pOTU", "vOTU"))
                    writer.writerows(heapq.merge(sorted_baseline, added))
                with added_path.open("w", encoding="utf-8", newline="") as handle:
                    writer = csv.writer(handle, delimiter="\t")
                    writer.writerow(("pOTU", "vOTU"))
                    writer.writerows(added)
                final_v = len(baseline_votus) + n_selected
                new_p = len(added_hosts - baseline_potus)
                final_p = len(baseline_potus) + new_p
                final_e = len(baseline_edges) + added_count
                actual_conn = final_e / (final_v * final_p)
                row = {
                    "network_file": network_path.relative_to(parent).as_posix(),
                    "mode": mode, "repetition": repetition,
                    "seed_base": args.seed_base, "selection_seed": selection_seed,
                    "quota_seed": quota_seed, "assignment_seed": assignment_seed, "host_seed": host_seed,
                    "baseline_votus": len(baseline_votus), "baseline_potus": len(baseline_potus),
                    "baseline_edges": len(baseline_edges), "new_virus_library_size": library_stats["complete_new_virus_library"],
                    "eligible_new_votus": len(eligible), "selected_new_votus": n_selected,
                    "new_votus_percent_of_baseline": 100.0 * n_selected / len(baseline_votus),
                    "added_edges": added_count, "new_potus": new_p,
                    "final_votus": final_v, "final_potus": final_p, "final_edges": final_e,
                    "baseline_connectance": base_conn, "final_connectance": actual_conn,
                    "connectance_ratio_to_baseline": actual_conn / base_conn,
                    "baseline_virus_mean_degree": base_mean,
                    "target_new_virus_mean_degree": mean_degree(target_hist),
                    "assigned_new_virus_mean_degree": mean_degree(assigned_hist),
                    "final_all_virus_mean_degree": final_e / final_v,
                    "baseline_virus_max_degree": max(baseline_hist),
                    "target_new_virus_max_degree": max(target_hist, default=0),
                    "assigned_new_virus_max_degree": max(assigned_hist, default=0),
                    "baseline_single_host_fraction": baseline_hist.get(1, 0) / len(baseline_votus),
                    "assigned_single_host_fraction": assigned_hist.get(1, 0) / n_selected if n_selected else 0.0,
                    "capacity_reduced_votus": reduced, "capacity_lost_edges": lost,
                    "target_vs_baseline_TV": target_tv, "target_vs_baseline_KS": target_ks,
                    "assigned_vs_baseline_TV": actual_tv, "assigned_vs_baseline_KS": actual_ks,
                    "assigned_vs_target_TV": capacity_tv, "assigned_vs_target_KS": capacity_ks,
                }
                if metrics_writer is None:
                    metrics_writer = csv.DictWriter(metrics_file, fieldnames=list(row), delimiter="\t")
                    metrics_writer.writeheader()
                metrics_writer.writerow(row)
                metrics_file.flush()
            print("{}: {} new viruses, {} new edges per mode, {} capacity-reduced quotas.".format(
                rep, n_selected, added_count, reduced
            ), flush=True)
    manifest["status"] = "complete"
    manifest["metrics"] = metrics_path.name
    manifest["repetitions_completed"] = args.repeats
    write_json(manifest_path, manifest)
    print("Finished. Metrics: " + str(metrics_path), flush=True)


def build_parser():
    parser = argparse.ArgumentParser(
        description="Expand OVPIN with new OVG viruses and OPG hosts in paired cumulative genus/family/order scenarios.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--fbaseline", required=True,
                        help="Original pOTU-vOTU edges; optional header; the four-column taxonomy file is accepted.")
    parser.add_argument("--fvotu-gvotu", required=True,
                        help="Current complete OVG vOTU-gvOTU map; optional vOTU/gvOTU header.")
    parser.add_argument("--fiphop", required=True,
                        help="vOTU and Python dictionary of predicted genera; raw or already OPG-filtered input.")
    parser.add_argument("--ftaxa", required=True,
                        help="Complete OPG taxonomy: 9 columns without header, or named pOTU/Order/Family/Genus columns.")
    parser.add_argument("--outfile", required=True,
                        help="NEW output prefix, e.g. results/expanded. Existing outputs are not overwritten.")
    selection = parser.add_mutually_exclusive_group()
    selection.add_argument("--addition", type=float,
                           help="New vOTUs as a percentage of ORIGINAL VIRUS count (not edge count). Omit both selection options to use all eligible new viruses.")
    selection.add_argument("--n-new-votus", type=int, help="Exact number of eligible new viruses per repetition.")
    parser.add_argument("--repeats", type=int, default=100, help="Number of paired three-mode repetitions.")
    parser.add_argument("--seed-base", type=int, default=1000)
    parser.add_argument("--prepare-only", action="store_true",
                        help="Write candidate library and input audit without generating networks.")
    parser.add_argument("--allow-missing-baseline-groups", action="store_true",
                        help="Explicitly allow baseline vOTUs absent from the new map; keep their edges but do not use them as group donors.")
    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.repeats < 1:
        parser.error("--repeats must be at least 1")
    if args.addition is not None and (not math.isfinite(args.addition) or args.addition < 0):
        parser.error("--addition must be finite and non-negative")
    if args.n_new_votus is not None and args.n_new_votus < 0:
        parser.error("--n-new-votus must be non-negative")
    try:
        run(args)
    except (OSError, ValueError) as exc:
        parser.exit(2, "Error: {}\n".format(exc))


if __name__ == "__main__":
    main()
