#!/usr/bin/env python3
from pathlib import Path
import subprocess
import sys
import re
import argparse
from collections import defaultdict

def run_cmd(cmd, **kwargs):
    print(" ".join(cmd))
    subprocess.run(cmd, check=True, **kwargs)

def concat_products(products_dir: Path):
    prod_files = sorted(products_dir.glob("*.fasta"))
    if not prod_files:
        sys.exit(f"ERROR: No .fasta files found in {products_dir}")
    all_fa = products_dir / "all_products.fasta"
    with all_fa.open("w") as out:
        for p in prod_files:
            text = p.read_text()
            out.write(text)
            if not text.endswith("\n"):
                out.write("\n")
    return all_fa

def build_index(minimap2_bin: str, refs: Path, mmi: Path, rebuild: bool):
    if rebuild and mmi.exists():
        mmi.unlink()
    if not mmi.exists():
        run_cmd([minimap2_bin, "-d", str(mmi), str(refs)])

def align(minimap2_bin: str, samtools_bin: str, mmi: Path, all_fa: Path, out_bam: Path):
    p1 = subprocess.Popen(
        [minimap2_bin, "-a", "-x", "sr", "--secondary=no", str(mmi), str(all_fa)],
        stdout=subprocess.PIPE,
        text=True
    )
    p2 = subprocess.Popen(
        [samtools_bin, "sort", "-o", str(out_bam)],
        stdin=p1.stdout,
        text=True
    )
    p1.stdout.close()
    p2.communicate()
    if p2.returncode != 0:
        sys.exit("ERROR: samtools sort failed.")
    run_cmd([samtools_bin, "index", str(out_bam)])

def compute_coverage(samtools_bin: str, bam: Path, out_cov: Path):
    with out_cov.open("w") as f:
        run_cmd([samtools_bin, "coverage", str(bam)], stdout=f)

def rank(refs: Path, bam: Path, cov_tsv: Path, out_ranked: Path, out_top: Path, samtools_bin: str, mapq_min: int):
    # Parse samtools coverage
    cov = {}
    with cov_tsv.open() as f:
        header = f.readline().strip().split("\t")
        idx = {name: i for i, name in enumerate(header)}
        for line in f:
            parts = line.rstrip("\n").split("\t")
            rname = parts[idx["#rname"]] if "#rname" in idx else parts[idx["rname"]]
            numreads = int(parts[idx["numreads"]])
            covbases = int(parts[idx["covbases"]])
            breadth = float(parts[idx["coverage"]])  # percent
            cov[rname] = {"numreads": numreads, "covbases": covbases, "breadth": breadth}

    # Identity proxy from NM tags for primary alignments passing MAPQ
    nm_sum = defaultdict(int)
    aln_sum = defaultdict(int)

    cig_re = re.compile(r"(\d+)([MIDNSHP=X])")
    nm_re = re.compile(r"\tNM:i:(\d+)")

    with subprocess.Popen([samtools_bin, "view", str(bam)], stdout=subprocess.PIPE, text=True) as p:
        for line in p.stdout:
            fields = line.rstrip("\n").split("\t")
            rname = fields[2]
            if rname == "*" or rname not in cov:
                continue
            mapq = int(fields[4])
            if mapq < mapq_min:
                continue

            cigar = fields[5]
            aln_len = 0
            for n, op in cig_re.findall(cigar):
                n = int(n)
                if op in ("M", "=", "X", "D"):
                    aln_len += n

            m = nm_re.search(line)
            nm = int(m.group(1)) if m else 0

            aln_sum[rname] += aln_len
            nm_sum[rname] += nm

    rows = []
    for rname, d in cov.items():
        aln_len = aln_sum[rname]
        nm = nm_sum[rname]
        identity = (1.0 - (nm / aln_len)) * 100.0 if aln_len > 0 else 0.0
        rows.append({
            "serotype": rname,
            "breadth_pct": d["breadth"],
            "aligned_bases": d["covbases"],
            "numreads": d["numreads"],
            "identity_proxy_pct": identity,
            "aln_len_for_identity": aln_len
        })

    rows.sort(key=lambda x: (x["breadth_pct"], x["aligned_bases"], x["identity_proxy_pct"]), reverse=True)

    with out_ranked.open("w") as out:
        out.write("\t".join([
            "serotype","breadth_pct","aligned_bases","numreads","identity_proxy_pct","aln_len_for_identity"
        ]) + "\n")
        for r in rows:
            out.write("\t".join([
                r["serotype"],
                f"{r['breadth_pct']:.2f}",
                str(r["aligned_bases"]),
                str(r["numreads"]),
                f"{r['identity_proxy_pct']:.2f}",
                str(r["aln_len_for_identity"]),
            ]) + "\n")

    with out_top.open("w") as out:
        out.write("Top serotype candidates (ranked):\n")
        for r in rows[:10]:
            out.write(
                f"- {r['serotype']}: breadth={r['breadth_pct']:.2f}%, "
                f"aligned_bases={r['aligned_bases']}, reads={r['numreads']}, "
                f"identity~={r['identity_proxy_pct']:.2f}%\n"
            )

def parse_args():
    ap = argparse.ArgumentParser(description="Serotype ranking from in-silico PCR products using minimap2+samtools.")
    ap.add_argument("--products", default="products_cycled", help="Folder containing product_*.fasta files.")
    ap.add_argument("--refs", default="cps_refs.fasta", help="cps reference FASTA (one locus per serotype).")
    ap.add_argument("--mmi", default="cps_refs.mmi", help="minimap2 index path for refs.")
    ap.add_argument("--rebuild-index", action="store_true", help="Force rebuild of the minimap2 index.")
    ap.add_argument("--mapq-min", type=int, default=20, help="Minimum MAPQ to count alignments toward identity proxy.")
    ap.add_argument("--minimap2", default="minimap2", help="Path/name of minimap2 binary.")
    ap.add_argument("--samtools", default="samtools", help="Path/name of samtools binary.")
    
    return ap.parse_args()



def main():
    args = parse_args()
    products_dir = Path(args.products)
    refs = Path(args.refs)
    mmi = Path(args.mmi)

    if not products_dir.exists():
        sys.exit(f"ERROR: Products folder not found: {products_dir}")
    if not refs.exists():
        sys.exit(f"ERROR: Reference FASTA not found: {refs}")

    all_fa = concat_products(products_dir)

    out_bam = products_dir / "products_vs_cps.bam"
    out_cov = products_dir / "coverage.tsv"
    out_ranked = products_dir / "serotype_ranked.tsv"
    out_top = products_dir / "serotype_top.txt"

    build_index(args.minimap2, refs, mmi, rebuild=True)
    align(args.minimap2, args.samtools, mmi, all_fa, out_bam)
    compute_coverage(args.samtools, out_bam, out_cov)
    rank(refs, out_bam, out_cov, out_ranked, out_top, args.samtools, args.mapq_min)

    print("\nDone.")
    print(f"Ranked results: {out_ranked}")
    print(f"Top summary:    {out_top}")
    print(f"BAM:            {out_bam}")

if __name__ == "__main__":
    main()
