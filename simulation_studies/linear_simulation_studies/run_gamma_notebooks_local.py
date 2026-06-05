"""
Run the Gamma benchmark Colab notebooks (folder33..folder48) LOCALLY, in
parallel, without GitHub. For each notebook we:
  * drop the `cell-install` cell (no devtools::install_github) and instead load
    the locally-installed package with library(zicbcf);
  * setwd() into a per-folder output directory so the relative OUT_CSV lands
    there -- this is required because the ZI-sensitivity notebooks (folders
    37-40, 45-48) reuse the same N=500 OUT_CSV filename as the standard
    notebooks (33, 41); isolating by folder prevents clobbering.

Each notebook is executed as a plain Rscript job; up to --jobs run concurrently
(default 6 = two folders' worth of A/B/C notebooks at a time).

Usage:
    python3 run_gamma_notebooks_local.py                 # all of 33..48, 6 at a time
    python3 run_gamma_notebooks_local.py --jobs 6
    python3 run_gamma_notebooks_local.py --folders 33 41 # just these folders (smoke test)

Outputs:
    CSVs   -> simulation_studies/results_gamma/folderNN/results_*.csv
    R + log-> <tmp>/gamma_local_runs/folderNN/<notebook>.R(.log)
"""
import argparse
import json
import os
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
NB_BASE = os.path.join(HERE, "colab_notebooks")
OUT_BASE = os.path.join(HERE, "results_gamma")
RUN_BASE = os.path.join(tempfile.gettempdir(), "gamma_local_runs")


def extract_r(nb_path, out_dir):
    """Notebook -> runnable R: setwd + library(zicbcf) + every code cell except install."""
    nb = json.load(open(nb_path, encoding="utf-8"))
    parts = [f"setwd({json.dumps(out_dir)})", "library(zicbcf)", ""]
    for c in nb["cells"]:
        if c["cell_type"] != "code" or c.get("id") == "cell-install":
            continue
        parts.append("".join(c["source"]))
    return "\n\n".join(parts)


def build_jobs(folders):
    jobs = []
    for fld in folders:
        nbdir = os.path.join(NB_BASE, fld)
        if not os.path.isdir(nbdir):
            print(f"  ! skipping missing {fld}", file=sys.stderr)
            continue
        out_dir = os.path.join(OUT_BASE, fld)
        run_dir = os.path.join(RUN_BASE, fld)
        os.makedirs(out_dir, exist_ok=True)
        os.makedirs(run_dir, exist_ok=True)
        for fn in sorted(os.listdir(nbdir)):
            if not fn.endswith(".ipynb"):
                continue
            r_path = os.path.join(run_dir, fn[:-6] + ".R")
            with open(r_path, "w", encoding="utf-8") as f:
                f.write(extract_r(os.path.join(nbdir, fn), out_dir))
            jobs.append((fld, fn, r_path))
    return jobs


def run_job(job):
    fld, fn, r_path = job
    log_path = r_path + ".log"
    with open(log_path, "w") as lf:
        rc = subprocess.run(["Rscript", r_path], stdout=lf,
                            stderr=subprocess.STDOUT).returncode
    tag = "OK " if rc == 0 else "ERR"
    print(f"  [{tag}] {fld}/{fn}  (Rscript exit {rc}; log: {log_path})", flush=True)
    return (fld, fn, rc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=6, help="concurrent notebooks (default 6)")
    ap.add_argument("--folders", nargs="*", type=int,
                    default=list(range(33, 49)), help="folder numbers (default 33..48)")
    args = ap.parse_args()

    folders = [f"folder{n}" for n in args.folders]
    jobs = build_jobs(folders)
    print(f"Running {len(jobs)} notebooks from {len(folders)} folders, "
          f"{args.jobs} at a time...\n")

    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        results = list(ex.map(run_job, jobs))

    bad = [f"{f}/{n}" for f, n, rc in results if rc != 0]
    print(f"\nDone: {len(results) - len(bad)}/{len(results)} notebooks exited cleanly.")
    if bad:
        print("Non-zero exits (check their .log):")
        for b in bad:
            print("  -", b)
    print(f"CSVs under: {OUT_BASE}/folderNN/")


if __name__ == "__main__":
    main()
