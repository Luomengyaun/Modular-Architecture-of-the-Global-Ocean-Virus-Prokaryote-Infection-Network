"""Pure selection and validation rules copied unchanged from the original workflow.

This module contains no scheduler, task submission or automatic recursion loop.
It validates exported results and selects eligible child networks. The original
function names, arguments, bodies and numerical rules are retained verbatim.
"""
from collections import defaultdict
import csv
import json
import math
from pathlib import Path

class WorkflowError(RuntimeError):
    pass


def read_json(path):
    with Path(path).open(encoding="utf-8-sig") as handle:
        return json.load(handle)


def read_network(path):
    """Return sorted unique string-ID pairs and input statistics; header required."""
    pairs = set()
    rows = blanks = 0
    header_seen = False
    with Path(path).open(encoding="utf-8-sig") as handle:
        for number, line in enumerate(handle, 1):
            fields = line.split()
            if not fields:
                blanks += 1
                continue
            if not header_seen:
                if [value.lower() for value in fields] != ["potu", "votu"]:
                    raise WorkflowError(f"{path}:{number}: expected header pOTU<TAB>vOTU.")
                header_seen = True
                continue
            if len(fields) != 2:
                raise WorkflowError(f"{path}:{number}: expected exactly two columns.")
            if [value.lower() for value in fields] == ["potu", "votu"]:
                raise WorkflowError(f"{path}:{number}: unexpected repeated header.")
            pairs.add(tuple(fields))
            rows += 1
    if not pairs:
        raise WorkflowError(f"No edges in {path}.")
    return sorted(pairs), {
        "input_rows": rows, "blank_rows": blanks,
        "duplicate_rows_removed": rows - len(pairs), "n_edges": len(pairs),
        "n_rows": len({p for p, v in pairs}), "n_cols": len({v for p, v in pairs}),
    }


def write_network(path, edges):
    with Path(path).open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["pOTU", "vOTU"])
        writer.writerows(sorted(edges))


def read_module_edges(path, expected_edges):
    """Validate the complete edge partition, including the -1 cross-module edges."""
    groups = defaultdict(list)
    stated_counts = {}
    seen = set()
    row_groups, col_groups = {}, {}
    with Path(path).open(encoding="utf-8-sig", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not {"row", "col", "module", "module_num"}.issubset(reader.fieldnames or []):
            raise WorkflowError(f"Invalid Module.txt header: {path}")
        for number, row in enumerate(reader, 2):
            try:
                edge = (row["row"], row["col"])
                module = int(row["module"])
                count = int(row["module_num"])
            except (ValueError, TypeError, KeyError) as exc:
                raise WorkflowError(f"Invalid module row {number} in {path}") from exc
            if module != -1 and module < 1:
                raise WorkflowError(f"Invalid module ID {module} in {path}")
            if edge in seen:
                raise WorkflowError(f"Duplicate edge in {path}: {edge}")
            if module in stated_counts and stated_counts[module] != count:
                raise WorkflowError(f"Inconsistent module_num in {path}")
            if module > 0:
                if edge[0] in row_groups and row_groups[edge[0]] != module:
                    raise WorkflowError(f"pOTU assigned to multiple internal modules: {edge[0]}")
                if edge[1] in col_groups and col_groups[edge[1]] != module:
                    raise WorkflowError(f"vOTU assigned to multiple internal modules: {edge[1]}")
                row_groups[edge[0]] = module
                col_groups[edge[1]] = module
            seen.add(edge)
            stated_counts[module] = count
            groups[module].append(edge)
    if seen != set(expected_edges):
        raise WorkflowError(f"Module.txt edge set differs from input: {path}")
    if any(len(edges) != stated_counts[group] for group, edges in groups.items()):
        raise WorkflowError(f"Incorrect module_num in {path}")
    return dict(groups)


def significance_decision(result, config, required=True):
    """Return (valid significant votes, unknown votes, sig/nsig/NA).

Legacy complete records need no new fields. Failed tests must be explicit;
an unlabelled missing/non-finite p is corrupt data, not a non-significant vote.
"""
    if not required:
        if result.get("null_status") != "skipped":
            raise WorkflowError("Expected explicitly skipped root null tests.")
        return None, 0, "not_tested"
    if result.get("null_status") not in ("complete", "partial"):
        raise WorkflowError("Null tests missing or incomplete; this is not a non-significant leaf.")
    models = result.get("null_models")
    if not isinstance(models, list) or len(models) != 5:
        raise WorkflowError("Expected five completed null models.")
    names = [item.get("name") for item in models]
    if len(set(names)) != 5 or set(names) != set(config["null_models"]):
        raise WorkflowError("Null models are missing, duplicated, or unexpected.")
    count = unknown = 0
    for item in models:
        observed = item.get("Qb_observed")
        if isinstance(observed, bool) or not isinstance(observed, (int, float)) or not math.isfinite(observed) or not math.isclose(observed, result["Qb"], rel_tol=1e-12, abs_tol=1e-12):
            raise WorkflowError("Null tests use a different observed Qb than the exported partition.")
        if item.get("trials") != config["optimization_trials"]:
            raise WorkflowError("Null optimization-trial count differs from configuration.")
        p = item.get("p_upper")
        status = item.get("status", "complete")
        if status == "not_applicable":
            if (item["name"] != "FIXED" or not (result.get("n_rows") == 1 or result.get("n_cols") == 1)
                    or p is not None or item.get("replicates") != 0 or not item.get("error")):
                raise WorkflowError("Only an explicitly explained single-dimension FIXED may be not_applicable.")
            continue
        if status == "failed":
            done = item.get("replicates")
            if (p is not None or type(done) is not int or not 0 <= done < config["null_replicates"]
                    or item.get("requested_replicates") != config["null_replicates"] or not item.get("error")):
                raise WorkflowError("A failed null test must have a missing p, counts and an error.")
            unknown += 1
            continue
        if status != "complete":
            raise WorkflowError("Unknown null-model status: " + str(status))
        if isinstance(p, bool) or not isinstance(p, (int, float)) or not math.isfinite(p) or not 0 <= p <= 1:
            raise WorkflowError(f"Invalid upper-tail p value for {item['name']}.")
        if item.get("replicates") != config["null_replicates"] or item.get("trials") != config["optimization_trials"]:
            raise WorkflowError("Null replicate or optimization-trial count differs from configuration.")
        count += p < config["alpha"]
    if unknown and result.get("null_status") != "partial":
        raise WorkflowError("Failed models require null_status=partial.")
    needed = config["min_significant_models"]
    decision = "sig" if count >= needed else "nsig" if count + unknown < needed else "NA"
    return count, unknown, decision


def significant_model_count(result, config, required=True):
    return significance_decision(result, config, required)[0]


def validate_node_result(folder, config):
    folder = Path(folder)
    job = read_json(folder / "job.json")
    result = read_json(folder / "result.json")
    edges, stats = read_network(folder / "network.tsv")
    for key in ("node_id", "level", "optimization_trials"):
        if result.get(key) != job[key]:
            raise WorkflowError(f"Result {key} differs from job settings: {folder}")
    for key in ("n_edges", "n_rows", "n_cols"):
        if result.get(key) != stats[key]:
            raise WorkflowError(f"Result {key} differs from input: {folder}")
    qb = result.get("Qb")
    if isinstance(qb, bool) or not isinstance(qb, (int, float)) or not math.isfinite(qb):
        raise WorkflowError(f"Missing or invalid Qb: {folder}")
    if type(result.get("n_modules")) is not int or result["n_modules"] < 1:
        raise WorkflowError(f"Invalid module count: {folder}")
    if not math.isclose(result.get("connectance", -1), stats["n_edges"] / (stats["n_rows"] * stats["n_cols"]), rel_tol=1e-10):
        raise WorkflowError(f"Invalid connectance: {folder}")
    votes, unknown, decision = significance_decision(result, config, required=job["run_null"])
    result.update(significant_models=votes, unknown_models=unknown, significance_status=decision)
    result["null_failed_models"] = [item["name"] for item in result.get("null_models", []) if item.get("status") == "failed"]
    groups = read_module_edges(folder / "Module.txt", edges)
    result["cross_module_edges"] = len(groups.get(-1, []))
    if len([k for k in groups if k > 0]) > result["n_modules"]:
        raise WorkflowError(f"Module IDs exceed reported module count: {folder}")
    return result, edges, groups


def plan_children(node, result, edges, groups, config):
    """Pure recursion rule; all detected modules are reported, eligible ones queued."""
    level = node["level"]
    threshold = math.sqrt(2 * len(edges))
    sig_count, _, decision = significance_decision(result, config, required=bool(level > 0 or config["run_root_null"]))
    parent_ok = decision in ("sig", "not_tested")
    positive_groups = {key: value for key, value in groups.items() if key > 0}
    genuine_split = result["n_modules"] >= 2 and len(positive_groups) >= 2
    records, selected = [], {}
    for module, child_edges in sorted(positive_groups.items()):
        child_id = str(module) if node["node_id"] == "root" else f"{node['node_id']}_{module}"
        rows = len({p for p, v in child_edges})
        cols = len({v for p, v in child_edges})
        reason = ""
        if decision == "NA":
            reason = "parent_significance_unknown"
        elif not parent_ok:
            reason = "parent_not_significant"
        elif not genuine_split or len(child_edges) >= len(edges):
            reason = "no_further_split"
        elif len(child_edges) < threshold:
            reason = "below_resolution_threshold"
        elif config["max_depth"] is not None and level + 1 > config["max_depth"]:
            reason = "depth_limit"
        record = {
            "node_id": child_id, "parent_id": node["node_id"], "level": level + 1,
            "local_module": module, "n_edges": len(child_edges), "n_rows": rows, "n_cols": cols,
            "parent_edges": len(edges), "resolution_threshold": threshold,
            "parent_significant_models": sig_count,
            "genuine_split": genuine_split, "selected_for_analysis": not reason, "stop_reason": reason,
        }
        records.append(record)
        if not reason:
            selected[child_id] = sorted(child_edges)
    if decision == "NA":
        stop_reason = "significance_unknown"
    elif not parent_ok:
        stop_reason = "not_significant"
    elif not genuine_split:
        stop_reason = "no_further_split"
    elif not selected:
        stop_reason = "depth_limit" if any(r["stop_reason"] == "depth_limit" for r in records) else "no_eligible_children"
    else:
        stop_reason = "expanded"
    return records, selected, stop_reason
