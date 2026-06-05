"""
Run the nonlinear-DGP Colab notebooks (colab_notebooks/folder1..40) LOCALLY,
without GitHub. For each notebook we:
  * drop the `cell-install` cell (no devtools::install_github) and load the
    locally-installed package with library(zicbcf);
  * setwd() into a single flat results directory so every relative OUT_CSV lands
    together. This is safe because the notebooks now use unique filenames --
    ZI-sensitivity files carry a `_lvl_<n>` suffix (e.g.
    results_bcf_dgp_a_N500_lvl_1.csv), so they never clobber the standard N=500
    file. (The R script + log for each notebook still live in a per-folder temp.)

Mirrors simulation_studies/run_gamma_notebooks_local.py.

Usage:
    python3 run_nonlinear_notebooks_local.py                  # all of 1..40, 6 at a time
    python3 run_nonlinear_notebooks_local.py --jobs 4
    python3 run_nonlinear_notebooks_local.py --folders 9 25   # just these folders
    python3 run_nonlinear_notebooks_local.py --folders 9 --smoke   # fast end-to-end check
                                                                   # (N_SIM=1, NBURN=NSIM=100)

Outputs:
    CSVs   -> nonlinear_simulation_studies/results_notebooks/results_*.csv  (flat)
    R+log  -> <tmp>/nonlinear_local_runs/folderNN/<notebook>.R(.log)
"""
import argparse
import json
import os
import re
import subprocess
import sys
import tempfile
from concurrent.futures import ThreadPoolExecutor

HERE = os.path.dirname(os.path.abspath(__file__))
NB_BASE = os.path.join(HERE, "colab_notebooks")
OUT_BASE = os.path.join(HERE, "results_notebooks")
RUN_BASE = os.path.join(tempfile.gettempdir(), "nonlinear_local_runs")


def extract_r(nb_path, out_dir, smoke=False):
    nb = json.load(open(nb_path, encoding="utf-8"))
    parts = [f"setwd({json.dumps(out_dir)})", "library(zicbcf)", ""]
    for c in nb["cells"]:
        if c["cell_type"] != "code" or c.get("id") == "cell-install":
            continue
        src = "".join(c["source"])
        if smoke:
            src = re.sub(r"N_SIM\s*<-\s*\d+L", "N_SIM   <- 1L", src)
            src = re.sub(r"NBURN\s*<-\s*\d+L", "NBURN   <- 100L", src)
            src = re.sub(r"NSIM\s*<-\s*\d+L", "NSIM    <- 100L", src)
        parts.append(src)
    return "\n\n".join(parts)


def build_jobs(folders, smoke):
    jobs = []
    os.makedirs(OUT_BASE, exist_ok=True)  # single flat results dir (unique filenames)
    for fld in folders:
        nbdir = os.path.join(NB_BASE, fld)
        if not os.path.isdir(nbdir):
            print(f"  ! skipping missing {fld}", file=sys.stderr)
            continue
        run_dir = os.path.join(RUN_BASE, fld)
        os.makedirs(run_dir, exist_ok=True)
        for fn in sorted(os.listdir(nbdir)):
            if not fn.endswith(".ipynb"):
                continue
            r_path = os.path.join(run_dir, fn[:-6] + ".R")
            with open(r_path, "w", encoding="utf-8") as f:
                f.write(extract_r(os.path.join(nbdir, fn), OUT_BASE, smoke=smoke))
            jobs.append((fld, fn, r_path))
    return jobs


def run_job(job):
    fld, fn, r_path = job
    log_path = r_path + ".log"
    with open(log_path, "w") as lf:
        rc = subprocess.run(["Rscript", r_path], stdout=lf, stderr=subprocess.STDOUT).returncode
    tag = "OK " if rc == 0 else "ERR"
    print(f"  [{tag}] {fld}/{fn}  (Rscript exit {rc}; log: {log_path})", flush=True)
    return (fld, fn, rc)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=6, help="concurrent notebooks (default 6)")
    ap.add_argument("--folders", nargs="*", type=int, default=list(range(1, 41)),
                    help="folder numbers (default 1..40)")
    ap.add_argument("--smoke", action="store_true",
                    help="fast end-to-end check: N_SIM=1, NBURN=NSIM=100")
    args = ap.parse_args()

    folders = [f"folder{n}" for n in args.folders]
    jobs = build_jobs(folders, args.smoke)
    print(f"Running {len(jobs)} notebooks from {len(folders)} folders, {args.jobs} at a time"
          f"{' [SMOKE]' if args.smoke else ''}...\n")

    with ThreadPoolExecutor(max_workers=args.jobs) as ex:
        results = list(ex.map(run_job, jobs))

    bad = [f"{f}/{n}" for f, n, rc in results if rc != 0]
    print(f"\nDone: {len(results) - len(bad)}/{len(results)} notebooks exited cleanly.")
    if bad:
        print("Non-zero exits (check their .log):")
        for b in bad:
            print("  -", b)
    print(f"CSVs under: {OUT_BASE}/")


if __name__ == "__main__":
    main()
