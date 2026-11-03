"""Exact small-population validation, or nested genotyping at fixed pedigree.

All compute stages run in separate R processes. No internal time limit is used;
the scheduler's limit still applies. Run from Eddie through followup_job.sh.
"""
import argparse
import csv
from datetime import datetime
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time

ROOT = Path(__file__).resolve().parents[1]
LARGE = ROOT / "validation" / "large_population.R"
WORKER = ROOT / "validation" / "followup_worker.R"


def positive(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("must be positive")
    return number


def measured(command, log):
    """Wait without a timeout, recording the same process-RSS convention as before."""
    start = time.monotonic()
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1",
               VECLIB_MAXIMUM_THREADS="1")
    with log.open("w") as stream:
        process = subprocess.Popen(command, cwd=ROOT, env=env, stdout=stream,
                                   stderr=stream, start_new_session=True)
        try:
            _, status, usage = os.wait4(process.pid, 0)
        except BaseException:
            try:
                os.killpg(process.pid, signal.SIGKILL)
            except ProcessLookupError:
                pass
            _, status, _ = os.wait4(process.pid, 0)
            process.returncode = os.waitstatus_to_exitcode(status)
            raise
        process.returncode = os.waitstatus_to_exitcode(status)
    return dict(status="ok" if process.returncode == 0 else "error",
                process_wall_sec=round(time.monotonic() - start, 3),
                peak_process_rss_mib=usage.ru_maxrss /
                    (1024 ** 2 if sys.platform == "darwin" else 1024),
                log=str(log))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("experiment", choices=["exact", "genotyped"])
    parser.add_argument("--package-root", type=Path, default=ROOT.parent / "randPedGenPCA")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--sizes", type=positive, nargs="+", default=[200, 500, 1000])
    parser.add_argument("--pedigree-size", type=positive, default=20000)
    parser.add_argument("--genotyped", type=positive, nargs="+", default=[500, 1000, 2000, 4000])
    parser.add_argument("--markers", type=positive, default=1000)
    parser.add_argument("--generations", type=positive, default=10)
    parser.add_argument("--geno-generations", type=positive, default=2)
    parser.add_argument("--rho", type=positive, default=5)
    parser.add_argument("--q", type=positive, default=4)
    parser.add_argument("--sketch", type=positive, default=20)
    parser.add_argument("--trace-vectors", type=positive, default=50)
    parser.add_argument("--replicates", type=positive)
    parser.add_argument("--seed", type=positive, default=123)
    parser.add_argument("--dry-run", action="store_true")
    a = parser.parse_args()
    a.package_root = a.package_root.resolve()
    a.replicates = a.replicates or (3 if a.experiment == "exact" else 1)
    a.sizes = sorted(set(a.sizes))
    a.genotyped = sorted(set(a.genotyped))
    if a.sketch < a.rho or a.rho < 2 or a.trace_vectors < 2:
        parser.error("require sketch >= rho >= 2 and trace-vectors >= 2")
    if a.geno_generations > a.generations:
        parser.error("geno-generations must not exceed generations")
    relevant_sizes = a.sizes if a.experiment == "exact" else [a.pedigree_size]
    if min(relevant_sizes) < max(4 * a.generations, 2 * a.sketch + 1):
        parser.error("populations must have >=4 per generation and exceed twice sketch")
    if a.experiment == "exact" and max(a.sizes) > 2000:
        parser.error("exact validation is deliberately limited to <=2000 individuals")
    generation_sizes = [a.pedigree_size // a.generations +
                        (i < a.pedigree_size % a.generations)
                        for i in range(a.generations)]
    eligible = sum(generation_sizes[-a.geno_generations:])
    if a.experiment == "genotyped" and max(a.genotyped) > eligible:
        parser.error(f"only {eligible} individuals eligible in last generations")
    if a.seed + max(relevant_sizes) + a.replicates > 2147483647:
        parser.error("seed exceeds R's integer range")
    sources = [a.package_root / "R" / f for f in ("utils.R", "rhpca.R", "hutchH.R")]
    if not all(p.is_file() for p in sources + [LARGE, WORKER]):
        parser.error("missing worker or package sources; check --package-root")
    configs = [(a.q, a.sketch)]
    if a.experiment == "exact":
        configs.append((a.q + 2, 2 * a.sketch))
    print(f"Experiment={a.experiment}; configs={configs}; H and CHC; "
          f"replicates={a.replicates}; markers={a.markers}", flush=True)
    print(f"Sizes={a.sizes}" if a.experiment == "exact" else
          f"Fixed pedigree={a.pedigree_size}; nested genotyped={a.genotyped}", flush=True)
    if a.dry_run:
        return
    if not shutil.which("Rscript") or not hasattr(os, "wait4"):
        parser.error("requires Rscript and macOS/Linux os.wait4")
    output = (a.output or ROOT / "output" /
              ("followup_" + a.experiment + "_" + datetime.now().strftime("%Y%m%d_%H%M%S_%f"))).resolve()
    output.mkdir(parents=True, exist_ok=False)
    manifest = {k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()}
    manifest.update(output=str(output), internal_timeout=None,
                    source_sha256={str(p): hashlib.sha256(p.read_bytes()).hexdigest()
                                   for p in sources + [LARGE, WORKER, Path(__file__).resolve()]})
    (output / "config.json").write_text(json.dumps(manifest, indent=2))
    rows = []
    print(f"Output: {output}", flush=True)

    def stage(name, key, script, args, csv_file=None, **metadata):
        usage = measured(["Rscript", str(script)] + [str(v) for v in args],
                         output / (key + ".log"))
        row = dict(stage=name, **metadata, **usage)
        if usage["status"] == "ok" and csv_file is not None:
            with csv_file.open() as stream:
                row.update(next(csv.DictReader(stream)))
        rows.append(row)
        fields = list(dict.fromkeys(k for item in rows for k in item))
        pending = output / "results.pending.csv"
        with pending.open("w") as stream:
            writer = csv.DictWriter(stream, fieldnames=fields)
            writer.writeheader()
            writer.writerows(rows)
        pending.replace(output / "results.csv")
        print(f"{name} {metadata} {usage['status']}: "
              f"{usage['process_wall_sec']}s; RSS={usage['peak_process_rss_mib']:.1f} MiB", flush=True)
        if usage["status"] != "ok":
            raise SystemExit(f"Stopped; see {usage['log']}")

    stage("session", "session", WORKER, ["session", output / "session_info.txt"])
    datasets = []
    if a.experiment == "exact":
        for n in a.sizes:
            dataset = output / f"data_{n}.rds"
            stage("prepare", f"prepare_{n}", LARGE,
                  ["prepare", dataset, n, a.markers, a.generations,
                   a.geno_generations, 0.5, a.seed + n], n=n)
            datasets.append((n, None, dataset))
    else:
        n = a.pedigree_size
        master = output / "master.rds"
        stage("prepare_master", "prepare_master", LARGE,
              ["prepare", master, n, a.markers, a.generations,
               a.geno_generations, 1, a.seed + n], n=n, n_genotyped=eligible)
        for count in a.genotyped:
            dataset = output / f"data_g{count}.rds"
            stage("subset", f"subset_g{count}", WORKER,
                  ["subset", master, dataset, count, a.seed], n=n, n_genotyped=count)
            datasets.append((n, count, dataset))

    for n, count, dataset in datasets:
        tag = f"n{n}" if count is None else f"n{n}_g{count}"
        for center in ("FALSE", "TRUE"):
            reference = output / f"{tag}_c{center}_exact.rds"
            if a.experiment == "exact":
                target = output / f"{tag}_c{center}_exact.csv"
                stage("exact", f"{tag}_c{center}_exact", WORKER,
                      ["reference", dataset, center, a.rho, reference, target],
                      csv_file=target, n=n, centredH=center)
            for rep in range(a.replicates):
                for q, sketch in configs:
                    key = f"{tag}_c{center}_q{q}_l{sketch}_rep{rep+1}"
                    fit = output / (key + ".rds")
                    target = output / (key + ".csv")
                    meta = dict(n=n, centredH=center, q=q, l=sketch, seed=a.seed + rep,
                                fit_file=str(fit))
                    if count is not None:
                        meta["n_genotyped"] = count
                    stage("pca", key, LARGE,
                          ["pca", a.package_root, dataset, target, fit,
                           q, sketch, a.rho, center, a.trace_vectors, a.seed + rep],
                          csv_file=target, **meta)
                    target = output / (key + "_diagnostics.csv")
                    if a.experiment == "exact":
                        stage("compare_exact", key + "_diagnostics", WORKER,
                              ["compare", fit, reference, target], csv_file=target, **meta)
                    else:
                        stage("diagnose", key + "_diagnostics", LARGE,
                              ["diagnose", a.package_root, dataset, fit, target, "none"],
                              csv_file=target, **meta)
    print(f"Finished. Results: {output / 'results.csv'}", flush=True)


if __name__ == "__main__":
    main()
