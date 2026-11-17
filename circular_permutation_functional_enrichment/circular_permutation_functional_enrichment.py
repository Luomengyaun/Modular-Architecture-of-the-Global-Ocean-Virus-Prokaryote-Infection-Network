#!/usr/bin/env python3
"""Test functional enrichment of virus–host matching regions by circular permutation."""

import argparse
import csv
import hashlib
import itertools
import math
import multiprocessing as mp
import os
import re
import time
from collections import defaultdict, deque
from concurrent.futures import ProcessPoolExecutor
from contextlib import contextmanager
from pathlib import Path

for variable in ("OPENBLAS_NUM_THREADS", "OMP_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(variable, "1")

import numpy as np


FUNCTIONS = {
    "head and packaging": "Head & packaging",
    "connector": "Connector",
    "tail": "Tail",
    "DNA, RNA and nucleotide metabolism": "Nucleotide metabolism",
    "integration and excision": "Integration & excision",
    "lysis": "Lysis",
    "transcription regulation": "Transcriptional regulation",
    "moron, auxiliary metabolic gene and host takeover": "AMG",
    "other": "Other",
    "unknown function": "Unknown function",
}
MODES = ("overlap_only", "flank3")
LEVELS = ("contig", "vOTU")
FUNCTION_KEYS = tuple(FUNCTIONS)
CATEGORY_INDEX = {name.casefold(): index for index, name in enumerate(FUNCTION_KEYS)}
INTERVAL = re.compile(r"(\d+)\s*-\s*(\d+)")


@contextmanager
def tsv_rows(path, required):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        yield checked_reader(handle, required, str(path))


def checked_reader(handle, required, source):
    rows = csv.DictReader(handle, delimiter="\t")
    columns = rows.fieldnames or []
    if len(set(columns)) != len(columns) or not set(required).issubset(columns):
        raise ValueError(f"{source}: required columns {required}, found {columns}")
    return rows


def write_tsv(path, header, rows):
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def merge_intervals(intervals):
    merged = []
    for start, end in sorted(intervals):
        if merged and start <= merged[-1][1] + 1:
            merged[-1] = (merged[-1][0], max(end, merged[-1][1]))
        else:
            merged.append((start, end))
    return tuple(merged)


def parse_intervals(value):
    intervals = []
    for token in re.split(r"[;,]", value):
        hit = INTERVAL.fullmatch(token.strip())
        if hit is None:
            raise ValueError(f"invalid nucleotide interval {token!r}")
        start, end = sorted(map(int, hit.groups()))
        if start < 1:
            raise ValueError(f"invalid 1-based coordinate {token!r}")
        intervals.append((start, end))
    return intervals


def read_retained(path):
    mapping = {}
    with tsv_rows(path, ("contig", "vOTU")) as rows:
        for row in rows:
            contig, votu = row["contig"].strip(), row["vOTU"].strip()
            if not contig or not votu or votu == "NA":
                raise ValueError(f"{path}:{rows.line_num}: missing contig/vOTU")
            if mapping.setdefault(contig, votu) != votu:
                raise ValueError(f"{path}: conflicting vOTU assignments for {contig}")
    if not mapping:
        raise ValueError(f"{path}: no retained contigs")
    return mapping


def read_qualified_pairs(path, mapping):
    pairs = defaultdict(dict)
    lengths = {}
    required = ("contig", "host", "merged_intervals", "union_length_bp",
                "contig_length_bp", "union_coverage_pct")
    with tsv_rows(path, required) as rows:
        for row in rows:
            try:
                contig, host = row["contig"].strip(), row["host"].strip()
                if contig not in mapping or not host:
                    raise ValueError(f"missing vOTU mapping or host for {contig!r}/{host!r}")
                if host in pairs[contig]:
                    raise ValueError(f"duplicate qualified pair {contig}/{host}")
                length = int(row["contig_length_bp"])
                if length <= 0 or lengths.setdefault(contig, length) != length:
                    raise ValueError(f"invalid or inconsistent contig length for {contig}")
                intervals = merge_intervals(parse_intervals(row["merged_intervals"]))
                bp = sum(end - start + 1 for start, end in intervals)
                reported_pct = number(row["union_coverage_pct"], 0, 10)
                if (intervals[-1][1] > length or bp != int(row["union_length_bp"])
                        or 10 * bp > length or abs(100 * bp / length - reported_pct) > 1e-6):
                    raise ValueError(f"interval bounds or coverage disagree for {contig}/{host}")
                pairs[contig][host] = intervals
            except ValueError as error:
                raise ValueError(f"{path}:{rows.line_num}: {error}") from error
    if not pairs:
        raise ValueError(f"{path}: no qualified virus–host pairs")
    return pairs, lengths


def build_contig_unions(pairs):
    return {contig: merge_intervals(itertools.chain.from_iterable(hosts.values()))
            for contig, hosts in pairs.items()}


def read_hq_intervals(target_file, contig_unions, lengths, mapping):
    data = {}
    with target_file.open(encoding="utf-8-sig") as source:
        for line in source:
            fields = line.split()
            if not fields or fields[0].casefold() in {"contig", "contig_id", "viral_contig", "qseqid"}:
                continue
            contig = fields[0]
            if contig in data:
                continue
            if contig not in contig_unions or contig not in mapping:
                raise ValueError(f"HQ contig absent from qualified pairs or vOTU mapping: {contig}")
            length = lengths[contig]
            data[contig] = mapping[contig], length, contig_unions[contig]
    if not data:
        raise ValueError("No HQ position intervals")
    return data


def number(value, lo, hi):
    x = float(value)
    if not math.isfinite(x) or not lo <= x <= hi:
        raise ValueError(f"Number outside [{lo}, {hi}]: {value}")
    return x


def target_contig(protein, hq, location):
    contig, _, ordinal = protein.rpartition("_")
    if contig not in hq:
        return None
    if not ordinal.isascii() or not ordinal.isdecimal():
        raise ValueError(f"{location}: invalid protein ID {protein!r}, expected contig_N")
    return contig


def read_hq_annotations(path, hq):
    """Read protein ID from column 1 and Pharokka category from column 5."""
    annotations = {}
    headers = {"id", "protein", "protein_id", "gene", "gene_id", "query", "query_id"}
    with path.open(encoding="utf-8-sig") as source:
        for line_no, line in enumerate(source, 1):
            if not line.strip():
                continue
            protein = line.split(None, 1)[0]
            if protein.casefold() in headers:
                continue
            location = f"{path}:{line_no}"
            if target_contig(protein, hq, location) is None:
                continue
            fields = line.rstrip("\r\n").split("\t")
            if len(fields) == 1:
                fields = re.split(r"\s{2,}", line.strip())
            if len(fields) >= 5:
                label = " ".join(fields[4].casefold().split())
            else:
                # With single-space formatting, the known category must be a full suffix.
                parts = line.split(None, 3)
                tail = " ".join(parts[3].casefold().split()) if len(parts) == 4 else ""
                label = next((key for key in sorted(CATEGORY_INDEX, key=len, reverse=True)
                              if tail.endswith(" " + key)), None)
                if label is None:
                    raise ValueError(f"{location}: use 5-column tab-separated Pharokka annotations")
            previous = annotations.setdefault(protein, label)
            if previous != label:
                raise ValueError(f"{location}: conflicting functional categories for {protein}")
    if not annotations:
        raise ValueError(f"{path}: no Pharokka protein IDs match the HQ contigs")
    return annotations


def find_faa_files(dirs, extra):
    paths = []
    for directory in dirs:
        found = sorted(path for path in directory.rglob("*.faa") if path.is_file())
        if not found:
            raise FileNotFoundError(f"No .faa files found in {directory}")
        paths.extend(found)
    for path in extra:
        if not path.is_file():
            raise FileNotFoundError(f"Missing Prodigal FAA: {path}")
        paths.append(path)
    return list(dict.fromkeys(path.resolve() for path in paths))


def read_hq_genes(hq, annotations, faa_files):
    genes = defaultdict(list)
    coords = {}
    for i, path in enumerate(faa_files, 1):
        print(f"C: reading Prodigal FAA {i}/{len(faa_files)}: {path}", flush=True)
        with path.open("rb") as source:
            for line_no, raw in enumerate(source, 1):
                if not raw.startswith(b">"):
                    continue
                fields = raw[1:].decode("utf-8").split("#", 4)
                tokens = fields[0].split()
                if not tokens:
                    continue
                protein = tokens[0]
                location = f"{path}:{line_no}"
                contig = target_contig(protein, hq, location)
                if contig is None:
                    continue
                try:
                    start, end, strand = map(int, fields[1:4])
                except ValueError as error:
                    raise ValueError(f"{location}: invalid Prodigal coordinates for {protein}") from error
                start, end = sorted((start, end))
                if not 1 <= start <= end <= hq[contig][1] or strand not in (-1, 1):
                    raise ValueError(f"{location}: CDS coordinate out of range for {protein}")
                coordinate = (start, end, strand)
                if protein in coords:
                    if coords[protein] != coordinate:
                        raise ValueError(f"{location}: conflicting coordinates for {protein}")
                    continue
                coords[protein] = coordinate
                category = CATEGORY_INDEX.get(annotations.get(protein, ""), -1)
                genes[contig].append((start, end, protein, category))
    missing = annotations.keys() - coords.keys()
    if missing:
        raise ValueError(f"{len(missing):,} annotated HQ proteins lack FAA coordinates. "
                         f"Examples: {', '.join(sorted(missing)[:5])}. Check --faa-dirs and --extra-faa.")
    missing_contigs = hq.keys() - genes.keys()
    if missing_contigs:
        raise ValueError(f"{len(missing_contigs):,} HQ contigs have no FAA genes. "
                         f"Examples: {', '.join(sorted(missing_contigs)[:5])}")
    return genes


def selection_for_shifts(starts, ends, intervals, length, shifts):
    """Find directly overlapping genes and up to three flanking genes on each side."""
    shifts = np.asarray(shifts, dtype=np.int64)
    batch, gene_count = len(shifts), len(starts)
    direct = np.zeros((batch, gene_count), dtype=bool)
    indices = np.arange(gene_count)[None, :]
    flank = np.zeros_like(direct)
    active_all = np.ones(batch, dtype=bool)
    ones = np.ones(batch, dtype=np.int64)

    def add_piece(left, right, active):
        overlap = (starts[None, :] <= right[:, None]) & (ends[None, :] >= left[:, None])
        direct[:] |= overlap & active[:, None]
        has_overlap = overlap.any(axis=1)
        insertion = np.searchsorted(starts, left, side="left")
        first = np.where(has_overlap, overlap.argmax(axis=1), insertion)[:, None]
        last = np.where(has_overlap, gene_count - 1 - overlap[:, ::-1].argmax(axis=1),
                        insertion - 1)[:, None]
        before = (indices >= first - 3) & (indices < first)
        after = (indices > last) & (indices <= last + 3)
        flank[:] |= (overlap | before | after) & active[:, None]

    for start, end in intervals:
        left = (start - 1 + shifts) % length + 1
        right = left + (end - start)
        add_piece(left, np.minimum(right, length), active_all)
        wrap = right > length
        if wrap.any():
            add_piece(ones, right - length, wrap)
    return direct, flank


def stable_seed(seed, key):
    digest = hashlib.sha256(f"{seed}|{key}".encode()).digest()
    return int.from_bytes(digest[:8], "little")


def bh_adjust(pvalues):
    valid = sorted((p, i) for i, p in enumerate(pvalues) if p is not None and math.isfinite(p))
    result = [None] * len(pvalues)
    running, n = 1.0, len(valid)
    for rank in range(n, 0, -1):
        p, i = valid[rank - 1]
        running = min(running, p * n / rank)
        result[i] = running
    return result


def summarize_null(observed, null):
    result = dict.fromkeys((
        "observed_pct", "random_mean_pct", "random_sd_pct", "random_q025_pct", "random_q975_pct",
        "difference_pp", "observed_over_random", "p_enrichment", "p_depletion", "p_two_sided",
    ))
    if not math.isfinite(observed):
        return result
    mean = float(null.mean())
    lo, hi = np.quantile(null, (0.025, 0.975))
    tolerance = 1e-9
    upper = (1 + int(np.count_nonzero(null >= observed - tolerance))) / (len(null) + 1)
    lower = (1 + int(np.count_nonzero(null <= observed + tolerance))) / (len(null) + 1)
    result.update(observed_pct=observed, random_mean_pct=mean, random_sd_pct=float(null.std()),
                  random_q025_pct=float(lo), random_q975_pct=float(hi),
                  difference_pp=observed - mean, observed_over_random=observed / mean if mean else None,
                  p_enrichment=upper, p_depletion=lower, p_two_sided=min(1.0, 2 * min(upper, lower)))
    return result


def simulate_votu(job):
    """One RNG stream per vOTU, with a shared fractional shift per replicate."""
    votu, members, n_random, seed, batch_size = job
    k = len(FUNCTION_KEYS)
    sum_null = np.zeros((n_random, len(MODES), k))
    sum_observed = np.zeros(sum_null.shape[1:])
    valid_members = np.zeros(k, dtype=int)
    fraction = np.random.default_rng(stable_seed(seed, votu)).random(n_random)
    member_scores = {}
    for contig, length, starts, ends, categories, intervals in members:
        onehot = (categories[:, None] == np.arange(k)).astype(np.float64)
        denominators = onehot.sum(axis=0)
        weights = 100 / np.maximum(denominators, 1)
        direct, flank = selection_for_shifts(starts, ends, intervals, length, [0])
        observed = (direct @ onehot)[0] * weights
        sum_observed[0] += observed
        sum_observed[1] += (flank @ onehot)[0] * weights
        valid_members += denominators > 0
        shifts = (fraction * length).astype(np.int64)
        direct_total = np.zeros(k)
        for first in range(0, n_random, batch_size):
            last = min(first + batch_size, n_random)
            direct, flank = selection_for_shifts(starts, ends, intervals, length, shifts[first:last])
            scores = (direct @ onehot) * weights
            direct_total += scores.sum(axis=0)
            sum_null[first:last, 0] += scores
            sum_null[first:last, 1] += (flank @ onehot) * weights
        member_scores[contig] = (denominators.astype(int), direct_total / n_random, observed)
    return votu, len(members), sum_observed, sum_null, valid_members, member_scores


def random_groups(jobs, workers):
    if workers == 1:
        yield from map(simulate_votu, jobs)
        return
    with ProcessPoolExecutor(max_workers=workers, mp_context=mp.get_context("spawn")) as pool:
        tasks = iter(jobs)
        pending = deque(pool.submit(simulate_votu, task) for task in itertools.islice(tasks, workers * 2))
        try:
            while pending:
                result = pending.popleft().result()
                next_task = next(tasks, None)
                if next_task is not None:
                    pending.append(pool.submit(simulate_votu, next_task))
                yield result
        finally:
            for future in pending:
                future.cancel()


def build_random_jobs(contigs, genes, args):
    grouped = defaultdict(list)
    for contig, (votu, length, regions) in sorted(contigs.items()):
        ordered = sorted(genes[contig])
        starts = np.array([row[0] for row in ordered], dtype=np.int64)
        ends = np.array([row[1] for row in ordered], dtype=np.int64)
        categories = np.array([row[3] for row in ordered], dtype=np.int16)
        grouped[votu].append((contig, length, starts, ends, categories, regions))
    return [(v, members, args.n_random, args.seed, args.batch_size)
            for v, members in sorted(grouped.items())]


def recompute_random(hq, genes, args):
    jobs = build_random_jobs(hq, genes, args)
    shape = (len(LEVELS), len(MODES), len(FUNCTION_KEYS))
    obs_sum = np.zeros(shape)
    null_sum = np.zeros((len(LEVELS), args.n_random, *shape[1:]))
    n_units = np.zeros((len(LEVELS), len(FUNCTION_KEYS)), dtype=int)
    box_rows = []
    contig_scores = {}
    done = 0
    started = next_log = time.monotonic()
    workers = min(args.workers, len(jobs))
    print(f"D: {args.n_random:,} shifts, seed={args.seed}, workers={workers}, "
          f"{len(hq):,} contigs and {len(jobs):,} vOTUs", flush=True)
    for votu, size, observed, null, counts, member_scores in random_groups(jobs, workers):
        contig_scores.update(member_scores)
        obs_sum[0] += observed
        null_sum[0] += null
        n_units[0] += counts
        scale = np.maximum(counts, 1)
        votu_observed = observed / scale
        votu_null = null / scale
        obs_sum[1] += votu_observed
        null_sum[1] += votu_null
        n_units[1] += counts > 0
        null_means = votu_null[:, 0, :].mean(axis=0)
        for fi, function in enumerate(FUNCTION_KEYS):
            if counts[fi] > 0:
                random, actual = np.clip([null_means[fi], votu_observed[0, fi]], 0, 100)
                box_rows.extend(((votu, FUNCTIONS[function], "Random", float(random)),
                                 (votu, FUNCTIONS[function], "OVIPN", float(actual))))
        done += size
        now = time.monotonic()
        if done == len(hq) or now >= next_log:
            eta = (len(hq) - done) * (now - started) / done
            print(f"D: {done:,}/{len(hq):,} HQ contigs, ETA {eta:.0f}s", flush=True)
            next_log = now + 5

    tests, distributions = [], {}
    for li, level in enumerate(LEVELS):
        for mi, mode in enumerate(MODES):
            for fi, function in enumerate(FUNCTION_KEYS):
                n = int(n_units[li, fi])
                observed = float(obs_sum[li, mi, fi] / n) if n else float("nan")
                null = null_sum[li, :, mi, fi] / n if n else np.full(args.n_random, np.nan)
                tests.append({"level": level, "mode": mode, "function": function, "n_units": n,
                              "n_random": args.n_random, **summarize_null(observed, null)})
                distributions[f"{level}_{mode}_{fi}"] = null
    for row, q in zip(tests, bh_adjust([row["p_two_sided"] for row in tests])):
        row["q_BH_global"] = q
    return box_rows, tests, distributions, contig_scores


def panel_d(hq, genes, args, out):
    box_rows, tests, nulls, contig_scores = recompute_random(hq, genes, args)
    figure_tests = [row for row in tests if row["level"] == "vOTU" and row["mode"] == "overlap_only"]
    stats = []
    for row in figure_tests:
        q, p, delta = row["q_BH_global"], row["p_two_sided"], row["difference_pp"]
        if q is None:
            raise ValueError(f"No HQ vOTU has a classified protein of function {row['function']!r}")
        stats.append((FUNCTIONS[row["function"]], row["n_units"], p, q, delta,
                      row["n_random"]))
    write_tsv(out / "D_function_boxplot.tsv", ("vOTU", "function", "group", "capture_pct"), box_rows)
    write_tsv(out / "D_function_significance.tsv",
              ("function", "n", "p_raw", "p_bh", "difference_pp", "n_random"), stats)
    columns = tuple(tests[0])
    write_tsv(out / "D_global_random_tests.tsv", columns,
              (["NA" if row[key] is None else row[key] for key in columns] for row in tests))
    np.savez_compressed(out / "D_group_null_distributions.npz", **nulls)
    n_genes = sum(len(genes[c]) for c in hq)
    n_unclassified = sum(row[3] < 0 for c in hq for row in genes[c])
    n_tests = sum(row["p_two_sided"] is not None for row in tests)
    write_tsv(out / "D_random_summary.tsv",
              ("analysis", "hq_contigs", "hq_votus", "genes", "unclassified_genes", "n_random",
               "seed", "tested_hypotheses", "pairs", "annotations", "n_union_segments"),
              [("recomputed", len(hq), len({v for v, _, _ in hq.values()}), n_genes, n_unclassified,
                args.n_random, args.seed, n_tests, args.pairs.name, args.annotations.name,
                sum(len(regions) for _, _, regions in hq.values()))])
    print(f"D: recomputed {n_tests} tests and exported {len(box_rows):,} boxplot values", flush=True)
    return contig_scores


def read_genes(contigs, annotation_path, faa_dirs, extra_faa, require_annotated_contigs=False):
    if not contigs:
        return {}
    annotations = read_hq_annotations(annotation_path, contigs)
    if require_annotated_contigs:
        annotated = {protein.rpartition("_")[0] for protein in annotations}
        missing = sorted(contigs.keys() - annotated)
        if missing:
            raise ValueError(f"{len(missing):,} HQ contigs have no protein annotations in {annotation_path}. "
                             f"Examples: {', '.join(missing[:5])}. Use --annotations with the updated table.")
    return read_hq_genes(contigs, annotations, find_faa_files(faa_dirs, extra_faa))


def default_workers():
    slurm = os.environ.get("SLURM_CPUS_PER_TASK", "")
    if slurm.isdigit() and int(slurm) > 0:
        return min(8, int(slurm))
    count = len(os.sched_getaffinity(0)) if hasattr(os, "sched_getaffinity") else (os.cpu_count() or 1)
    return min(8, count)



def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pairs", type=Path, default=Path("pair_union_le10pct.tsv"))
    parser.add_argument("--retained", type=Path, default=Path("retained_contig_votu.tsv"))
    parser.add_argument("--hq-contigs", type=Path, default=Path("matched_HQ_comp90_contigs.txt"))
    parser.add_argument("--annotations", type=Path, default=Path("retained_contig_votu_phrokka.txt"))
    parser.add_argument("--faa-dirs", type=Path, nargs="+", default=[], help="search these directories for .faa files")
    parser.add_argument("--extra-faa", type=Path, nargs="+", default=[], help="explicit Prodigal .faa files")
    parser.add_argument("--n-random", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261111)
    parser.add_argument("--workers", type=int, default=default_workers())
    parser.add_argument("--batch-size", type=int, default=128)
    parser.add_argument("--output-dir", type=Path, default=Path("functional_enrichment_results"))
    args = parser.parse_args()
    if not (args.faa_dirs or args.extra_faa):
        parser.error("provide --faa-dirs, --extra-faa, or both")
    if min(args.n_random, args.workers, args.batch_size) <= 0:
        parser.error("--n-random, --workers, and --batch-size must be positive")

    mapping = read_retained(args.retained)
    pairs, lengths = read_qualified_pairs(args.pairs, mapping)
    hq = read_hq_intervals(args.hq_contigs, build_contig_unions(pairs), lengths, mapping)
    genes = read_genes(hq, args.annotations, args.faa_dirs, args.extra_faa, require_annotated_contigs=True)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    panel_d(hq, genes, args, args.output_dir)


if __name__ == "__main__":
    main()
