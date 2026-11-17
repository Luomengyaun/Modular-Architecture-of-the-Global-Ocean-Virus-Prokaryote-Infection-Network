"""Validate saved root text exports without loading MAT or recomputing a partition."""
import csv
import math
from pathlib import Path


def _labels(path):
    lines = path.read_text(encoding="utf-8-sig").splitlines()
    if len(lines) != 1:
        raise ValueError(f"{path.name} must contain exactly one tab-separated label line.")
    labels = lines[0].split("\t")
    if (any(not x or any(c.isspace() for c in x) for x in labels)
            or len(set(labels)) != len(labels)):
        raise ValueError(f"{path.name} contains empty, whitespace-bearing or duplicate labels.")
    return labels


def _pair(edge):
    if (not isinstance(edge, (tuple, list)) or len(edge) != 2
            or not all(isinstance(x, str) and x for x in edge)):
        raise ValueError("An exported edge must be a pair of nonempty string labels.")
    return tuple(edge)


def validate_exports(root_folder, observed, edges, groups=None):
    """Return ordered labels/modules and a read-only scientific validation report.

    ``edges`` is the normalized list already read from root/network.tsv.
    Optional ``groups`` is the validated Module.txt mapping, including -1.
    Recomputed Qb/Qr validate the exports only; saved observed scalars stay intact.
    """
    root = Path(root_folder)
    for key in ("n_rows", "n_cols", "n_edges", "n_modules"):
        if type(observed.get(key)) is not int or observed[key] < 1:
            raise ValueError(f"Invalid saved observed {key}.")
    for key in ("Qb", "Qr", "NODF"):
        value = observed.get(key)
        if key == "NODF" and value is None:
            continue
        if type(value) not in (int, float) or not math.isfinite(value):
            raise ValueError(f"Invalid saved observed {key}.")
    rows, cols = _labels(root / "label_row.txt"), _labels(root / "label_col.txt")
    if (len(rows), len(cols)) != (observed["n_rows"], observed["n_cols"]):
        raise ValueError("Label dimensions differ from the saved root result.")
    with (root / "node_modules.tsv").open(encoding="utf-8-sig", newline="") as handle:
        records = list(csv.reader(handle, delimiter="\t"))
    expected = [("pOTU", x) for x in rows] + [("vOTU", x) for x in cols]
    if (not records or records[0] != ["node_type", "node_id", "module"]
            or len(records) != len(expected) + 1):
        raise ValueError("Invalid node_modules.tsv header or node count.")
    modules, total_modules = [], observed["n_modules"]
    for number, (record, identity) in enumerate(zip(records[1:], expected), 2):
        if len(record) != 3 or tuple(record[:2]) != identity:
            raise ValueError(f"Partition labels/order differ from label exports at line {number}.")
        try:
            module = int(record[2])
        except ValueError as exc:
            raise ValueError(f"Invalid partition module at line {number}.") from exc
        if not 1 <= module <= total_modules:
            raise ValueError(f"Partition module outside 1..N at line {number}.")
        modules.append(module)
    if len(set(modules)) != total_modules:
        raise ValueError("Saved module count differs from the exported partition.")
    row_modules, col_modules = modules[:len(rows)], modules[len(rows):]
    row_map, col_map = dict(zip(rows, row_modules)), dict(zip(cols, col_modules))
    row_totals, col_totals = [0] * (total_modules + 1), [0] * (total_modules + 1)
    seen, used_rows, used_cols, inside = set(), set(), set(), 0
    for raw_edge in edges:
        edge = _pair(raw_edge)
        if edge in seen:
            raise ValueError(f"Duplicate exported edge: {edge}.")
        if edge[0] not in row_map or edge[1] not in col_map:
            raise ValueError(f"Exported edge has an unknown label: {edge}.")
        seen.add(edge)
        used_rows.add(edge[0]); used_cols.add(edge[1])
        rm, cm = row_map[edge[0]], col_map[edge[1]]
        row_totals[rm] += 1; col_totals[cm] += 1
        inside += rm == cm
    if len(seen) != observed["n_edges"] or len(used_rows) != len(rows) or len(used_cols) != len(cols):
        raise ValueError("Exported edges do not cover the saved edge/node counts.")
    if groups is not None:
        grouped = set()
        for module, members in groups.items():
            if type(module) is not int or (module != -1 and not 1 <= module <= total_modules):
                raise ValueError("Invalid Module.txt group ID.")
            for raw_edge in members:
                edge = _pair(raw_edge)
                if edge not in seen or edge in grouped:
                    raise ValueError("Module.txt groups contain an unknown or repeated edge.")
                rm, cm = row_map[edge[0]], col_map[edge[1]]
                if module != (rm if rm == cm else -1):
                    raise ValueError(f"Module.txt disagrees with node_modules.tsv for {edge}.")
                grouped.add(edge)
        if grouped != seen:
            raise ValueError("Module.txt groups do not cover the source edges.")
    count = len(seen)
    product = sum(k * d for k, d in zip(row_totals, col_totals))
    qb, qr = inside / count - product / (count * count), 2 * inside / count - 1
    for key, value in (("Qb", qb), ("Qr", qr)):
        if not math.isclose(value, observed[key], rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError(f"Recomputed exported-partition {key} differs from the saved observation.")
    summary = {key: observed[key] for key in ("n_rows", "n_cols", "n_edges", "n_modules")}
    summary.update(internal_edges=inside, degree_product_sum=product,
                   Qb_recomputed=qb, Qr_recomputed=qr,
                   Qb_difference=qb - observed["Qb"], Qr_difference=qr - observed["Qr"],
                   Qb_observed=observed["Qb"], Qr_observed=observed["Qr"],
                   NODF_observed=observed.get("NODF"), tolerance=1e-12)
    return {"row_labels": rows, "col_labels": cols, "row_modules": row_modules,
            "col_modules": col_modules, "summary": summary}
