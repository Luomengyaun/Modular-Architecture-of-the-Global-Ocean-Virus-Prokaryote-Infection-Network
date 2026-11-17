"""Prepare direct MATLAB job inputs; never calculate statistics or submit jobs.

Paths supplied as arguments are resolved against the public package directory.
Generated job files stay in example_output/, excluded from version control.
"""
import argparse
import hashlib
import json
from pathlib import Path
import sys

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from core_rules import read_json, read_network, write_network, read_module_edges
from root_null_exports import validate_exports

LIBRARY_FILES = (
    "main/Bipartite.m", "main/Options.m", "base/MatrixFunctions.m",
    "community/AdaptiveBrim.m", "community/BipartiteModularity.m",
    "statistics/NullModels.m", "nestedness/NestednessNODF.m",
    "nestedness/NODF.m", "nestedness/Nestedness.m",
)
MODELS = ["EQUIPROBABLE", "AVERAGE", "AVERAGE_COLS", "AVERAGE_ROWS", "FIXED"]


def path(value):
    value = Path(value).expanduser()
    return (value if value.is_absolute() else ROOT / value).resolve()


def digest(file):
    return hashlib.sha256(file.read_bytes()).hexdigest()


def save_job(folder, job):
    target = folder / "job.json"
    with target.open("w", encoding="utf-8", newline="\n") as handle:
        handle.write(json.dumps(job, indent=2, ensure_ascii=False, allow_nan=False) + "\n")
    print(target)
    return target


def new_output(value):
    folder = path(value)
    if folder.exists() and (not folder.is_dir() or any(folder.iterdir())):
        raise ValueError("Use a new or empty output directory: " + str(folder))
    return folder


def prepare_module(args):
    config = read_json(path(args.config))
    if len(config["null_models"]) != 5 or set(config["null_models"]) != set(MODELS):
        raise ValueError("The module example requires the original five distinct models.")
    if args.level < 0 or (args.node_id == "root") != (args.level == 0):
        raise ValueError("Use root/level 0, or a child node ID with a positive level.")
    if args.level > 0 and not args.root_network:
        raise ValueError("Provide the original normalized --root-network for a child seed.")
    library = path(args.bimat_path)
    if not (library / "main/Bipartite.m").is_file():
        raise ValueError("BiMat directory is missing main/Bipartite.m.")
    edges, _ = read_network(path(args.network))
    root_bytes = None
    if args.root_network:
        read_network(path(args.root_network))
        root_bytes = path(args.root_network).read_bytes()
    folder = new_output(args.output)
    folder.mkdir(parents=True, exist_ok=True)
    network = folder / "network.tsv"
    write_network(network, edges)
    root_hash = hashlib.sha256(root_bytes).hexdigest() if root_bytes is not None else digest(network)
    seed_text = "{}:{}:{}".format(config["seed"], root_hash, args.node_id)
    seed = int(hashlib.sha256(seed_text.encode()).hexdigest()[:8], 16)
    job = dict(config, node_id=args.node_id, level=args.level, seed=seed,
               network_file=str(network), output_dir=str(folder), bimat_path=str(library),
               cpus=args.cpus, run_null=bool(args.level > 0 or config["run_root_null"]),
               run_ntc=bool(config["run_ntc"] and args.level > 0))
    return save_job(folder, job)


def prepare_root(args):
    source = path(args.observed)
    config = read_json(path(args.config))
    observed = read_json(source / "result.json")
    saved_job = read_json(source / "job.json")
    allowed = saved_job["null_models"]
    if len(allowed) != 5 or set(allowed) != set(MODELS):
        raise ValueError("The saved observation must record the five original models.")
    if observed["node_id"] != "root" or saved_job["node_id"] != "root":
        raise ValueError("The observation must be the saved root result.")
    if observed["seed"] != saved_job["seed"]:
        raise ValueError("The saved root result and job have different seeds.")
    if observed["optimization_trials"] != saved_job["optimization_trials"]:
        raise ValueError("The saved root result and job have different trial counts.")
    if observed["optimization_trials"] != config["optimization_trials"]:
        raise ValueError("Null trials must match the saved observed trials.")
    names = args.models
    if not names or len(set(names)) != len(names) or any(name not in allowed for name in names):
        raise ValueError("Select distinct supported model names.")
    repeats = args.replicates if args.replicates is not None else config["null_replicates"]
    if repeats < 1:
        raise ValueError("Replicate count must be positive.")
    edges, _ = read_network(source / "network.tsv")
    groups = read_module_edges(source / "Module.txt", edges)
    validate_exports(source, observed, edges, groups)
    library = path(args.bimat_path)
    folder = new_output(args.output)
    job = dict(config, node_id="root", level=0, input_mode="edge_partition",
               network_file=str(source / "network.tsv"),
               label_row_file=str(source / "label_row.txt"),
               label_col_file=str(source / "label_col.txt"),
               partition_file=str(source / "node_modules.tsv"),
               output_dir=str(folder), bimat_path=str(library),
               seed=observed["seed"], cpus=args.cpus, null_replicates=repeats,
               null_models=allowed, global_model_index=[allowed.index(name)+1 for name in allowed if name in names])
    job.update({"observed_" + key: observed.get(key) for key in
                ("Qb", "Qr", "NODF", "n_rows", "n_cols", "n_edges", "n_modules")})
    definition = dict(job)
    for key in ("output_dir", "cpus", "global_model_index"):
        definition.pop(key)
    definition["source_sha256"] = {name: digest(source / name) for name in
        ("network.tsv", "job.json", "result.json", "Module.txt", "node_modules.tsv", "label_row.txt", "label_col.txt")}
    definition["code_sha256"] = {name: digest(ROOT / name) for name in
        ("root_null_worker.m", "NullAdaptiveBrim.m", "bimat_null_modularity.m",
         "bimat_nestedness_fast.m", "bimat_validate_fast_null.m")}
    definition["library_sha256"] = {name: digest(library / name) for name in LIBRARY_FILES}
    job["fingerprint"] = hashlib.sha256(json.dumps(definition, sort_keys=True, allow_nan=False).encode()).hexdigest()
    folder.mkdir(parents=True, exist_ok=True)
    return save_job(folder, job)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="action", required=True)
    module = commands.add_parser("module", help="Prepare one observed network/module calculation.")
    module.add_argument("--network", default="examples/network.tsv")
    module.add_argument("--node-id", default="root")
    module.add_argument("--level", type=int, default=0)
    module.add_argument("--root-network", help="The original normalized root/network.tsv, byte for byte.")
    module.add_argument("--output", default="example_output/observed")
    root = commands.add_parser("root", help="Prepare selected root null models from saved text outputs.")
    root.add_argument("--observed", default="example_output/observed")
    root.add_argument("--models", nargs="+", default=["FIXED"])
    root.add_argument("--replicates", type=int)
    root.add_argument("--output", default="example_output/root_null")
    for command in (module, root):
        command.add_argument("--config", default="analysis_config.json")
        command.add_argument("--bimat-path", required=True)
        command.add_argument("--cpus", type=int, default=1, help="MATLAB thread cap, not a batch allocation.")
    args = parser.parse_args(argv)
    if args.cpus < 1:
        parser.error("--cpus must be positive")
    return prepare_module(args) if args.action == "module" else prepare_root(args)


if __name__ == "__main__":
    main()
