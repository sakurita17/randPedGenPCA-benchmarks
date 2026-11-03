"""Large matrix-free exercise. Run --help for editable sizes and PCA budgets."""
import argparse
import csv
import hashlib
import json
import os
from pathlib import Path
import shutil
import signal
import subprocess
import sys
import time
from datetime import datetime

ROOT = Path(__file__).resolve().parents[1]
WORKER = ROOT / 'validation' / 'large_population.R'
FIELDS = ['stage', 'status', 'n', 'n_genotyped', 'markers', 'rho', 'q', 'l',
          'centredH', 'seed', 'trace_vectors', 'simulation_sec', 'pca_elapsed_sec',
          'diagnostic_elapsed_sec', 'process_wall_sec', 'peak_process_rss_mib',
          'symmetry_error', 'sketch_rank', 'orth_fro_error', 'score_scaling_max_error', 'score_scaling',
          'centred_scores_max_mean', 'trace_estimate', 'explained_sum',
          'max_relative_eigenpair_residual', 'eigenvalue_relative_change_vs_base',
          'normalised_subspace_change_vs_base', 'warnings', 'log', 'fit_file']


def positive(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError('must be positive')
    return number


def nonnegative(value):
    number = int(value)
    if number < 0:
        raise argparse.ArgumentTypeError('must be non-negative')
    return number


def run_process(command, log, timeout):
    start = time.monotonic()
    env = os.environ.copy()
    env.update(OMP_NUM_THREADS='1', OPENBLAS_NUM_THREADS='1', VECLIB_MAXIMUM_THREADS='1')
    with log.open('w') as stream:
        process = subprocess.Popen(command, cwd=ROOT, env=env, stdout=stream,
                                   stderr=stream, start_new_session=True)
        timed_out = False
        try:
            while True:
                pid, status, usage = os.wait4(process.pid, os.WNOHANG)
                if pid:
                    break
                if time.monotonic() - start > timeout:
                    timed_out = True
                    os.killpg(process.pid, signal.SIGKILL)
                    _, status, usage = os.wait4(process.pid, 0)
                    break
                time.sleep(0.2)
        except BaseException:
            os.killpg(process.pid, signal.SIGKILL)
            _, status, _ = os.wait4(process.pid, 0)
            process.returncode = os.waitstatus_to_exitcode(status)
            raise
        process.returncode = os.waitstatus_to_exitcode(status)
    return dict(status='timeout' if timed_out else ('ok' if process.returncode == 0 else 'error'),
                process_wall_sec=round(time.monotonic() - start, 3),
                peak_process_rss_mib=usage.ru_maxrss / (1024**2 if sys.platform == 'darwin' else 1024),
                log=str(log))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--sizes', type=positive, nargs='+', default=[10000, 50000, 100000])
    parser.add_argument('--markers', type=positive, default=1000)
    parser.add_argument('--generations', type=positive, default=10, help='includes founder generation')
    parser.add_argument('--geno-generations', type=positive, default=2)
    parser.add_argument('--geno-fraction', type=float, default=0.5, help='fraction in each eligible generation')
    parser.add_argument('--rho', type=positive, default=5)
    parser.add_argument('--q', type=nonnegative, default=4)
    parser.add_argument('--sketch', type=positive, default=20)
    parser.add_argument('--trace-vectors', type=nonnegative, default=50, help='0 disables Hutch++')
    parser.add_argument('--centering', choices=['both', 'false', 'true'], default='both')
    parser.add_argument('--compare', action='store_true', help='also run q+2 and twice the sketch size')
    parser.add_argument('--replicates', type=positive, default=1)
    parser.add_argument('--seed', type=positive, default=123)
    parser.add_argument('--timeout', type=positive, default=1800, help='seconds per process, including startup')
    parser.add_argument('--package-root', type=Path, default=ROOT.parent / 'RandPedGenPCA' / 'randPedGenPCA')
    parser.add_argument('--output', type=Path, default=None)
    parser.add_argument('--keep-datasets', action='store_true')
    parser.add_argument('--dry-run', action='store_true')
    a = parser.parse_args()
    if not hasattr(os, 'wait4'):
        parser.error('peak RSS measurement requires macOS or Linux (os.wait4)')
    if not 0 < a.geno_fraction <= 1 or a.geno_generations > a.generations:
        parser.error('require 0 < geno-fraction <= 1 and geno-generations <= generations')
    if min(a.sizes) < 4 * a.generations or a.sketch < a.rho:
        parser.error('require >= 4 individuals per generation and sketch >= rho')
    if a.rho >= min(a.sizes) or a.trace_vectors == 1:
        parser.error('require rho < smallest population and trace-vectors = 0 or >= 2')
    if a.seed + max(a.sizes) + a.replicates > 2147483647:
        parser.error('seed plus sizes/replicates exceeds R integer range')
    a.package_root = a.package_root.resolve()
    sources = [a.package_root / 'R' / name for name in ('utils.R', 'rhpca.R', 'hutchH.R')]
    if not all(p.is_file() for p in sources):
        parser.error('package R sources missing; set --package-root')
    centers = ['FALSE', 'TRUE'] if a.centering == 'both' else [a.centering.upper()]
    configs = [(a.q, a.sketch)] + ([(a.q + 2, a.sketch * 2)] if a.compare else [])
    sizes = sorted(set(a.sizes))
    print(f'Sizes={sizes}; markers={a.markers}; rho={a.rho}; configs={configs}; centered={centers}', flush=True)
    print(f'Eligible generations={a.geno_generations}/{a.generations}; sampled fraction={a.geno_fraction}; '
          f'trace vectors={a.trace_vectors}; replicates={a.replicates}', flush=True)
    if a.dry_run:
        return
    if not shutil.which('Rscript'):
        parser.error('Rscript is not on PATH')
    output = (a.output or ROOT / 'output' / ('large_' + datetime.now().strftime('%Y%m%d_%H%M%S_%f'))).resolve()
    output.mkdir(parents=True, exist_ok=False)
    manifest = {k: str(v) if isinstance(v, Path) else v for k, v in vars(a).items()}
    manifest['output'] = str(output)
    manifest['source_sha256'] = {str(p): hashlib.sha256(p.read_bytes()).hexdigest() for p in sources}
    (output / 'config.json').write_text(json.dumps(manifest, indent=2))
    print(f'Output: {output}', flush=True)
    rows = []

    def save(row):
        rows.append(row)
        pending = output / 'results.pending.csv'
        with pending.open('w') as stream:
            writer = csv.DictWriter(stream, fieldnames=FIELDS)
            writer.writeheader(); writer.writerows(rows)
        pending.replace(output / 'results.csv')
        print(f'{row["stage"]} n={row["n"]} {row["status"]}: '
              f'{row.get("pca_elapsed_sec", row["process_wall_sec"])}s; '
              f'peak process RSS={row["peak_process_rss_mib"]:.1f} MiB', flush=True)

    # Separate processes prevent generation and diagnostic peaks contaminating PCA RSS.
    for n in sizes:
        dataset = output / f'data_{n}.rds'
        result = run_process(['Rscript', str(WORKER), 'prepare', str(dataset), str(n),
                   str(a.markers), str(a.generations), str(a.geno_generations),
                   str(a.geno_fraction), str(a.seed + n)], output / f'prepare_{n}.log', a.timeout)
        save(dict(stage='prepare', n=n, markers=a.markers, seed=a.seed+n, **result))
        if result['status'] != 'ok':
            raise SystemExit(f'Stopped after preparation failure. See {result["log"]}')
        for center in centers:
            for rep in range(a.replicates):
                base_fit = None
                for q, sketch in configs:
                    key = f'n{n}_c{center}_q{q}_l{sketch}_rep{rep+1}'
                    csv_file = output / f'{key}.csv'; fit_file = output / f'{key}.rds'
                    result = run_process(['Rscript', str(WORKER), 'pca', str(a.package_root), str(dataset),
                        str(csv_file), str(fit_file), str(q), str(sketch), str(a.rho), center,
                        str(a.trace_vectors), str(a.seed + rep)], output / f'{key}.log', a.timeout)
                    row = dict(stage='pca', n=n, q=q, l=sketch, centredH=center, seed=a.seed+rep,
                               fit_file=str(fit_file), **result)
                    if result['status'] == 'ok':
                        with csv_file.open() as stream:
                            row.update(next(csv.DictReader(stream)))
                    save(row)
                    if result['status'] != 'ok':
                        raise SystemExit(f'Stopped after PCA failure. See {result["log"]}')
                    diagnostic = output / f'{key}_diagnostics.csv'
                    result = run_process(['Rscript', str(WORKER), 'diagnose', str(a.package_root), str(dataset),
                        str(fit_file), str(diagnostic), str(base_fit) if base_fit else 'none'],
                        output / f'{key}_diagnostics.log', a.timeout)
                    row = dict(stage='diagnose', n=n, q=q, l=sketch, centredH=center, seed=a.seed+rep, **result)
                    if result['status'] == 'ok':
                        with diagnostic.open() as stream:
                            row.update(next(csv.DictReader(stream)))
                    save(row)
                    if result['status'] != 'ok':
                        raise SystemExit(f'Stopped after diagnostic failure. See {result["log"]}')
                    if base_fit is None:
                        base_fit = fit_file
        if not a.keep_datasets:
            dataset.unlink()
    print(f'Finished. Results: {output / "results.csv"}', flush=True)


if __name__ == '__main__':
    main()
