#!/usr/bin/env python3
from pathlib import Path
from collections import defaultdict
import argparse
import re

FASTA_EXTENSIONS = {".fa", ".fasta", ".fna"}

_rc_map = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def revcomp(seq: str) -> str:
    return seq.translate(_rc_map)[::-1]

def read_fasta(path: Path):
    header = None
    seq = []
    with path.open() as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(seq).upper()
                header = line[1:].strip()
                seq = []
            else:
                seq.append(line.strip())
        if header is not None:
            yield header, "".join(seq).upper()

def write_single_fasta(path: Path, header: str, sequence: str):
    with path.open("w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), 80):
            f.write(sequence[i:i+80] + "\n")

def find_with_mismatches(seq: str, primer: str, max_mismatches: int, three_prime_lock: int):
    """
    Return list of (start_pos, mismatches) where primer matches seq[start:start+len(primer)]
    with <= max_mismatches total mismatches AND exact match in last `three_prime_lock` bases.
    """
    hits = []
    plen = len(primer)

    if three_prime_lock < 0:
        three_prime_lock = 0
    if three_prime_lock > plen:
        three_prime_lock = plen

    for i in range(0, len(seq) - plen + 1):
        window = seq[i:i+plen]

        # Enforce exact match at 3′ end (last bases of primer)
        if three_prime_lock and window[-three_prime_lock:] != primer[-three_prime_lock:]:
            continue

        mismatches = 0
        for a, b in zip(window, primer):
            if a != b:
                mismatches += 1
                if mismatches > max_mismatches:
                    break

        if mismatches <= max_mismatches:
            hits.append((i, mismatches))

    return hits

bed_line_re = re.compile(
    r"^\s*(?!#)(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+([+-])\s+([ACGTNacgtn]+)"
)

def load_primers(bed_path: Path):
    primers_by_amp = defaultdict(lambda: {"LEFT": [], "RIGHT": []})
    with bed_path.open() as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            m = bed_line_re.match(line)
            if not m:
                continue
            chrom, start, end, name, amp_s, strand, primer_seq = m.groups()
            amp = int(amp_s)
            primer_seq = primer_seq.upper()

            if "_LEFT" in name:
                side = "LEFT"
            elif "_RIGHT" in name:
                side = "RIGHT"
            else:
                side = "LEFT" if strand == "+" else "RIGHT"

            # Search sequence on the forward-oriented template string
            # If primer is defined on '-' strand, its binding site appears as revcomp(primer) in the forward template.
            search_seq = primer_seq if strand == "+" else revcomp(primer_seq)

            primers_by_amp[amp][side].append({
                "name": name,
                "strand": strand,
                "primer_seq": primer_seq,
                "search_seq": search_seq,
                "len": len(primer_seq),
            })
    return primers_by_amp

def best_product_for_amplicon(
    template_seq: str,
    left_primers,
    right_primers,
    min_len: int,
    max_len: int,
    max_mismatches: int,
    three_prime_lock: int,
):
    # Each hit is (pos, primer_dict, mismatches)
    left_hits = []
    right_hits = []

    for p in left_primers:
        for pos, mm in find_with_mismatches(template_seq, p["search_seq"], max_mismatches, three_prime_lock):
            left_hits.append((pos, p, mm))

    for p in right_primers:
        for pos, mm in find_with_mismatches(template_seq, p["search_seq"], max_mismatches, three_prime_lock):
            right_hits.append((pos, p, mm))

    if not left_hits or not right_hits:
        return None, len(left_hits), len(right_hits)

    # Best candidate chosen by:
    # 1) fewest total mismatches (L+R)
    # 2) shortest product length
    # 3) earliest left position (tie-break)
    best = None  # (total_mm, product_len, lpos, rpos, lp, rp, lp_mm, rp_mm, product_seq)

    for lpos, lp, lp_mm in left_hits:
        for rpos, rp, rp_mm in right_hits:
            if rpos <= lpos:
                continue

            product_end = rpos + rp["len"]
            product_len = product_end - lpos

            if product_len < min_len or product_len > max_len:
                continue

            total_mm = lp_mm + rp_mm
            product_seq = template_seq[lpos:product_end]
            # Key for ordering: only sortable types
            key = (total_mm, product_len, lpos, rpos, lp_mm, rp_mm, lp["name"], rp["name"])

            # Store both the key and the payload
            cand = (key, lp, rp, product_seq)

            if best is None or key < best[0]:
                best = cand


    return best, len(left_hits), len(right_hits)

def run(
    templates_dir: Path,
    bed_path: Path,
    out_dir: Path,
    min_product_len: int = 1,
    max_product_len: int = 20000,
    max_mismatches: int = 0,
    three_prime_lock: int = 8,
):
    out_dir.mkdir(exist_ok=True)

    primers_by_amp = load_primers(bed_path)
    amplicons = sorted(primers_by_amp.keys())

    template_files = sorted(
        [p for p in templates_dir.iterdir() if p.suffix.lower() in FASTA_EXTENSIONS],
        key=lambda p: p.name.lower()
    )

    report_path = out_dir / "amplification_report.tsv"
    product_index = 1

    with report_path.open("w") as rep:
        rep.write("\t".join([
            "product_id", "template_file", "template_record", "amplicon",
            "status", "product_length",
            "left_primer_name", "right_primer_name",
            "left_pos_0based", "right_pos_0based",
            "left_hits", "right_hits",
            "left_mismatches", "right_mismatches", "total_mismatches",
        ]) + "\n")

        for tf in template_files:
            for record_header, seq in read_fasta(tf):
                for amp in amplicons:
                    L = primers_by_amp[amp]["LEFT"]
                    R = primers_by_amp[amp]["RIGHT"]

                    best, left_hits, right_hits = best_product_for_amplicon(
                        seq, L, R,
                        min_product_len, max_product_len,
                        max_mismatches=max_mismatches,
                        three_prime_lock=three_prime_lock
                    )

                    if best is None:
                        continue

                    key, lp, rp, product_seq = best
                    total_mm, product_len, lpos, rpos, lp_mm, rp_mm, lp_name, rp_name = key


                    pid = f"product_{product_index:06d}"
                    write_single_fasta(out_dir / f"{pid}.fasta", pid, product_seq)

                    rep.write("\t".join([
                        pid, tf.name, record_header, str(amp),
                        "OK", str(product_len),
                        lp["name"], rp["name"],
                        str(lpos), str(rpos),
                        str(left_hits), str(right_hits),
                        str(lp_mm), str(rp_mm), str(total_mm),
                    ]) + "\n")

                    product_index += 1

    print(f"Done. Products: {product_index-1}")
    print(f"Products folder: {out_dir.resolve()}")
    print(f"Report: {report_path.resolve()}")

def parse_args():
    ap = argparse.ArgumentParser(description="In-silico PCR using ARTIC-style primer.bed against template FASTAs.")
    ap.add_argument("--templates", default="holocene", help="Folder containing template FASTA files.")
    ap.add_argument("--bed", default="primal_resources/primer.bed", help="Primer BED file.")
    ap.add_argument("--out", default="products", help="Output folder for products and report.")
    ap.add_argument("--min-len", type=int, default=450, help="Minimum allowed product length (bp).")
    ap.add_argument("--max-len", type=int, default=650, help="Maximum allowed product length (bp).")
    ap.add_argument("--max-mismatches", type=int, default=2, help="Max mismatches per primer.")
    ap.add_argument("--three-prime-lock", type=int, default=8, help="Last N bases at 3′ end must match exactly.")
    return ap.parse_args()

if __name__ == "__main__":
    args = parse_args()
    run(
        Path(args.templates),
        Path(args.bed),
        Path(args.out),
        min_product_len=args.min_len,
        max_product_len=args.max_len,
        max_mismatches=args.max_mismatches,
        three_prime_lock=args.three_prime_lock,
    )
