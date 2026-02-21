#!/usr/bin/env python3
"""
amplification.py

Post-process insilico.py output to simulate stochastic thermocycling and produce
a cycled products folder (resampled FASTA files) suitable for downstream scripts
that expect a directory of FASTAs.

Inputs:
  - amplification_report.tsv (from insilico.py)
  - products directory containing product_XXXXXX.fasta files

Outputs:
  - products_cycled/ directory with N output FASTA files (resampled)
  - cycled_counts.tsv: counts per source product_id
  - cycled_index.tsv: mapping cycled file -> source product_id
  - cycled_params.json: run parameters

Design notes:
  - We do NOT generate new sequences. PCR amplifies the same product sequences.
  - We simulate relative abundances after C cycles, then resample N output reads/files
    proportional to those abundances.
  - Stochasticity is applied in log-space to avoid astronomical molecule counts.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import csv
import json
import math
import random
from collections import Counter
import matplotlib.pyplot as plt


FASTA_EXTENSIONS = {".fa", ".fasta", ".fna"}

def binomial_draw(n: int, p: float, rng: random.Random) -> int:
    """Fast Binomial(n,p) without numpy: exact for small n, normal approx for large n."""
    if n <= 0 or p <= 0.0:
        return 0
    if p >= 1.0:
        return n

    # Exact for small n (keeps early-cycle stochasticity accurate)
    if n < 200:
        return sum(1 for _ in range(n) if rng.random() < p)

    mu = n * p
    var = n * p * (1.0 - p)
    if var < 1e-9:
        # essentially deterministic
        k = int(round(mu))
        return max(0, min(n, k))

    # Normal approximation
    k = int(round(rng.gauss(mu, math.sqrt(var))))
    return max(0, min(n, k))


def poisson_draw(lam: float, rng: random.Random) -> int:
    """Poisson sampler (Knuth). Good for small lam (<~10), which is what you want here."""
    if lam <= 0:
        return 0
    L = math.exp(-lam)
    k = 0
    p = 1.0
    while p > L:
        k += 1
        p *= rng.random()
    return k - 1


def simulate_curves(p_list, start_lambda: float, eff_decay: float, ct_threshold: float, cycles, replicates, seed,
                    baseline_mu=50.0, baseline_sigma=5.0, K=1e8):
    rng = random.Random(seed)  # seed can be None for random
    curves = []
    cts = []

    for _ in range(replicates):
        # Sparse starting templates: most targets start at 0 copies, some at 1–few.
        Ns = [poisson_draw(start_lambda, rng) for _ in p_list]

        F = []
        ct = None

        for t in range(cycles + 1):
            # total before this cycle's amplification step
            total_N = sum(Ns)

            # Amplify (skip at t = 0)
            if t > 0:
                sat = max(0.0, 1.0 - total_N / K)
                for i, p in enumerate(p_list):
                    n = Ns[i]
                    # Late-cycle drag / inhibition: efficiency decays with cycle number
                    cycle_drag = math.exp(-eff_decay * (t - 1))

                    effective_p = p * sat * cycle_drag

                    new = binomial_draw(n, effective_p, rng)
                    Ns[i] = n + new

                total_N = sum(Ns)  # recompute after amplification

            # Observe fluorescence (baseline + noise)
            baseline = rng.gauss(baseline_mu, baseline_sigma)
            fluorescence = max(0.0, total_N + baseline)
            F.append(fluorescence)

            # Ct call based on fluorescence
            if ct is None and fluorescence >= ct_threshold:
                ct = t

        curves.append(F)
        cts.append(ct if ct is not None else -1)

    return curves, cts


def write_curves_csv(out_path, curves):
    with open(out_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["replicate", "cycle", "fluorescence"])
        for r, F in enumerate(curves):
            for t, val in enumerate(F):
                w.writerow([r, t, val])

def plot_curves_png(out_path, curves):
    plt.figure()
    for F in curves[:20]:  # plot first 20 to avoid spaghetti
        plt.plot(range(len(F)), F)
    plt.yscale("log")
    plt.xlabel("Cycle")
    plt.ylabel("Fluorescence (arb, log scale)")
    plt.title("Synthetic qPCR amplification curves")
    plt.tight_layout()
    plt.savefig(out_path, dpi=200)
    plt.close()

    


def read_single_fasta(path: Path) -> tuple[str, str]:
    """Read a single-record FASTA (as written by insilico.py)."""
    header = None
    seq_parts = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                header = line[1:].strip()
            else:
                seq_parts.append(line.strip())
    if header is None:
        raise ValueError(f"No FASTA header found in {path}")
    return header, "".join(seq_parts).upper()

def write_single_fasta(path: Path, header: str, sequence: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i+80] + "\n")

def parse_report(report_path: Path) -> list[dict]:
    """Read amplification_report.tsv into a list of dict rows."""
    rows = []
    with report_path.open() as f:
        reader = csv.DictReader(f, delimiter="\t")
        for row in reader:
            # Minimal sanity checks
            if not row.get("product_id"):
                continue
            rows.append(row)
    return rows

def clamp(x: float, lo: float, hi: float) -> float:
    return max(lo, min(hi, x))

def compute_p_i(
    total_mismatches: int,
    product_length: int,
    p_max: float,
    mm_factor: float,
    len_ref: int,
    len_penalty: float,
) -> float:
    """
    Heuristic per-product "success probability" per cycle.

    - Start at p_max for a perfect match near len_ref
    - Multiply down by mm_factor^mismatch_count
    - Penalize length deviations above len_ref using exp(-len_penalty * (L-len_ref)/1000)
    """
    p = p_max * (mm_factor ** total_mismatches)
    if product_length > len_ref:
        p *= math.exp(-len_penalty * (product_length - len_ref) / 1000.0)
    # Keep within [0, 1]
    return clamp(p, 0.0, 1.0)

def weighted_sample_indices(weights: list[float], k: int, rng: random.Random) -> list[int]:
    """
    Sample k indices with replacement according to weights.
    Uses cumulative distribution for simplicity.
    """
    total = sum(weights)
    if total <= 0:
        raise ValueError("All weights are zero; cannot sample.")
    cdf = []
    acc = 0.0
    for w in weights:
        acc += w / total
        cdf.append(acc)

    out = []
    for _ in range(k):
        r = rng.random()
        # linear scan is fine for typical product counts (hundreds to thousands)
        # could switch to bisect if you want
        for i, c in enumerate(cdf):
            if r <= c:
                out.append(i)
                break
    return out

def main():
    ap = argparse.ArgumentParser(
        description="Simulate stochastic thermocycling on insilico PCR products and generate products_cycled/."
    )
    ap.add_argument("--report", default="products/amplification_report.tsv",
                    help="Path to amplification_report.tsv from insilico.py.")
    ap.add_argument("--products-dir", default="products",
                    help="Directory containing product_*.fasta created by insilico.py.")
    ap.add_argument("--out-dir", default="products_cycled",
                    help="Output directory for cycled FASTAs + reports.")
    ap.add_argument("--cycles", type=int, default=35, help="Number of PCR cycles to simulate.")
    ap.add_argument("--n-out", type=int, default=10000,
                    help="How many cycled FASTA files to generate (resampled proportional to abundance).")

    # Efficiency model parameters
    ap.add_argument("--p-max", type=float, default=0.92,
                    help="Max per-cycle success probability for a perfect-match product.")
    ap.add_argument("--mm-factor", type=float, default=0.80,
                    help="Per-mismatch multiplicative penalty (p *= mm_factor^mismatches).")
    ap.add_argument("--len-ref", type=int, default=500,
                    help="Reference amplicon length (bp) where length penalty starts.")
    ap.add_argument("--len-penalty", type=float, default=1.2,
                    help="Penalty strength for longer products; higher = harsher penalty.")

    # Stochasticity knobs (log-space)
    ap.add_argument("--sigma-per-cycle", type=float, default=0.06,
                    help="Log-space noise per cycle (std). Overall noise ~ sigma_per_cycle*sqrt(cycles).")
    ap.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility.")

    # ct graph
    ap.add_argument("--replicates", type=int, default=50,
                help="How many synthetic qPCR replicates to simulate for Ct curves.")
    

    ap.add_argument("--start_lambda", type=float, default=0.8,
                        help="Mean starting copies per target (Poisson). Use <1 for sparse templates.")
    ap.add_argument("--eff_decay", type=float, default=0.12,
                        help="Per-cycle exponential decay of efficiency. Higher = slower late-cycle growth.")
    ap.add_argument("--ct_threshold", type=float, default=1e6,
                        help="Fluorescence threshold for Ct call (raise to make Ct later / absent).")



    args = ap.parse_args()

    report_path = Path(args.report)
    products_dir = Path(args.products_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    rows = parse_report(report_path)
    if not rows:
        raise SystemExit(f"No rows found in report: {report_path}")

    # Build product metadata list
    products = []
    for r in rows:
        pid = r["product_id"]
        try:
            total_mm = int(r.get("total_mismatches", "0"))
        except ValueError:
            total_mm = 0
        try:
            plen = int(r.get("product_length", "0"))
        except ValueError:
            plen = 0

        fasta_path = products_dir / f"{pid}.fasta"
        if not fasta_path.exists():
            # Some users may rename/move; skip missing
            continue

        p_i = compute_p_i(
            total_mismatches=total_mm,
            product_length=plen,
            p_max=args.p_max,
            mm_factor=args.mm_factor,
            len_ref=args.len_ref,
            len_penalty=args.len_penalty,
        )

        products.append({
            "product_id": pid,
            "product_length": plen,
            "total_mismatches": total_mm,
            "p_i": p_i,
            "fasta_path": str(fasta_path),
        })

    if not products:
        raise SystemExit("No product FASTAs found matching report product_id values.")

    rng = random.Random(args.seed)

    # Compute simulated post-cycling weights in log space:
    # weight_i ~ exp( cycles * log(1+p_i) + Normal(0, sigma*sqrt(cycles)) )
    weights = []
    for pr in products:
        p_i = pr["p_i"]
        # Guard against log(1+0)
        growth = math.log1p(p_i)
        noise = rng.gauss(0.0, args.sigma_per_cycle * math.sqrt(max(args.cycles, 1)))
        log_w = args.cycles * growth + noise
        w = math.exp(log_w)
        weights.append(w)

    # Resample N output FASTAs according to weights
    chosen_idx = weighted_sample_indices(weights, args.n_out, rng)
    counts = Counter(chosen_idx)

    # Write outputs
    # 1) cycled_counts.tsv
    counts_path = out_dir / "cycled_counts.tsv"
    with counts_path.open("w") as f:
        f.write("\t".join([
            "source_product_id", "cycled_count",
            "p_i", "product_length", "total_mismatches",
        ]) + "\n")
        for i, pr in enumerate(products):
            c = counts.get(i, 0)
            if c == 0:
                continue
            f.write("\t".join([
                pr["product_id"], str(c),
                f"{pr['p_i']:.6f}",
                str(pr["product_length"]),
                str(pr["total_mismatches"]),
            ]) + "\n")

    # 2) cycled_index.tsv + FASTAs
    index_path = out_dir / "cycled_index.tsv"
    with index_path.open("w") as idxf:
        idxf.write("\t".join(["cycled_file", "source_product_id"]) + "\n")

        out_i = 1
        for i in chosen_idx:
            pr = products[i]
            src_pid = pr["product_id"]
            src_fa = Path(pr["fasta_path"])
            _, seq = read_single_fasta(src_fa)

            cycled_pid = f"cycled_{out_i:07d}"
            out_fa = out_dir / f"{cycled_pid}.fasta"
            header = f"{cycled_pid} source={src_pid} cycles={args.cycles} p_i={pr['p_i']:.4f}"
            write_single_fasta(out_fa, header, seq)

            idxf.write(f"{out_fa.name}\t{src_pid}\n")
            out_i += 1

    # 3) params json
    params_path = out_dir / "cycled_params.json"
    with params_path.open("w") as f:
        json.dump({
            "report": str(report_path),
            "products_dir": str(products_dir),
            "out_dir": str(out_dir),
            "cycles": args.cycles,
            "n_out": args.n_out,
            "p_max": args.p_max,
            "mm_factor": args.mm_factor,
            "len_ref": args.len_ref,
            "len_penalty": args.len_penalty,
            "sigma_per_cycle": args.sigma_per_cycle,
            "seed": args.seed,
            "n_input_products_used": len(products),
        }, f, indent=2)

    print(f"Done.")
    print(f"Input products used: {len(products)}")
    print(f"Cycled FASTAs written: {args.n_out}")
    print(f"Output folder: {out_dir.resolve()}")
    print(f"Counts: {counts_path.resolve()}")
    print(f"Index:  {index_path.resolve()}")
    print(f"Params: {params_path.resolve()}")

    p_list = [pr["p_i"] for pr in products]
    curves, cts = simulate_curves(
        p_list=p_list,
        start_lambda=args.start_lambda,
        eff_decay=args.eff_decay,
        ct_threshold=args.ct_threshold,
        cycles=args.cycles,
        replicates=args.replicates,
        seed=args.seed,
        baseline_mu=50.0,
        baseline_sigma=5.0,
        K=1e8,
    )


    write_curves_csv(out_dir / "ct_curves.csv", curves)
    with (out_dir / "ct_values.tsv").open("w") as f:
        f.write("replicate\tct\n")
        for i, ct in enumerate(cts):
            f.write(f"{i}\t{ct}\n")

    plot_curves_png(out_dir / "ct_curves.png", curves)



if __name__ == "__main__":
    main()
