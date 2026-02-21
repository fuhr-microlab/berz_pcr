#!/usr/bin/env python3
from __future__ import annotations

from pathlib import Path
from collections import defaultdict
import argparse
import gzip
import re
import sys

FASTA_EXTENSIONS = {".fa", ".fasta", ".fna"}

_rc_map = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def revcomp(seq: str) -> str:
    return seq.translate(_rc_map)[::-1]

def open_fastq(path: Path):
    if str(path).endswith(".gz"):
        return gzip.open(path, "rt")
    return path.open("r")

def iter_fastq_seqs(handle):
    """Yield (read_index_1based, seq_string_upper)."""
    i = 0
    while True:
        h = handle.readline()
        if not h:
            return
        seq = handle.readline()
        plus = handle.readline()
        qual = handle.readline()
        if not qual:
            return
        i += 1
        yield i, seq.strip().upper()

bed_line_re = re.compile(
    r"^\s*(?!#)(\S+)\s+(\d+)\s+(\d+)\s+(\S+)\s+(\d+)\s+([+-])\s+([ACGTNacgtn]+)"
)

def load_primers(bed_path: Path):
    """
    Returns list of primer dicts:
      {name, amp, strand, side, seq, len}
    """
    primers = []
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

            primers.append({
                "name": name,
                "amplicon": amp,
                "strand": strand,
                "side": side,
                "seq": primer_seq,
                "len": len(primer_seq),
            })
    return primers

def mismatches_at(seq: str, primer: str, start: int, max_mm: int) -> int | None:
    """Return mismatch count if <= max_mm else None."""
    mm = 0
    for a, b in zip(seq[start:start+len(primer)], primer):
        if a != b:
            mm += 1
            if mm > max_mm:
                return None
    return mm

def build_seed_index(primers, seed_len: int):
    """
    Build mapping: seed -> list of (primer_idx, offset_in_primer)
    where seed == primer[offset:offset+seed_len]
    We include multiple seeds per primer to increase recall.
    """
    idx = defaultdict(list)
    for pi, p in enumerate(primers):
        s = p["seq"]
        L = p["len"]
        if L < seed_len:
            continue

        # Use several seeds (start, end, middle) for recall
        offsets = {0, L - seed_len, (L - seed_len)//2}
        for off in sorted(offsets):
            seed = s[off:off+seed_len]
            idx[seed].append((pi, off))

    return idx

def scan_one_orientation(seq: str, primers, seed_index, seed_len: int, max_mm: int, stats, seen_in_read):
    """
    Update stats for primer hits found in seq.
    stats[primer_idx] = dict counters
    seen_in_read is a set of primer_idx already counted for reads_with_hit.
    """
    n = len(seq)
    if n < seed_len:
        return

    # For speed, scan seeds in read, use seed_index to propose candidates
    for pos in range(0, n - seed_len + 1):
        seed = seq[pos:pos+seed_len]
        cand = seed_index.get(seed)
        if not cand:
            continue

        for primer_idx, off in cand:
            p = primers[primer_idx]
            start = pos - off
            if start < 0:
                continue
            end = start + p["len"]
            if end > n:
                continue

            mm = mismatches_at(seq, p["seq"], start, max_mm)
            if mm is None:
                continue

            # Count hit
            st = stats[primer_idx]
            st["total_hits"] += 1
            st["sum_mismatches"] += mm

            if primer_idx not in seen_in_read:
                st["reads_with_hit"] += 1
                seen_in_read.add(primer_idx)

def run(fastq_path: Path, bed_path: Path, out_dir: Path, max_mismatches: int = 2, seed_len: int = 8, limit_reads: int = 0):
    out_dir.mkdir(parents=True, exist_ok=True)

    primers = load_primers(bed_path)
    if not primers:
        raise SystemExit(f"No primers parsed from BED: {bed_path}")

    seed_index = build_seed_index(primers, seed_len=seed_len)
    if not seed_index:
        raise SystemExit(f"Seed index empty (seed_len={seed_len}). Are primers shorter than seed_len?")

    # Stats per primer index
    stats = []
    for _ in primers:
        stats.append({
            "reads_with_hit": 0,
            "total_hits": 0,
            "sum_mismatches": 0,
        })

    total_reads = 0

    with open_fastq(fastq_path) as fh:
        for ridx, seq in iter_fastq_seqs(fh):
            total_reads += 1
            if limit_reads and total_reads > limit_reads:
                break

            # Track primers already seen in this read (across both orientations)
            seen_in_read = set()

            # Forward
            scan_one_orientation(seq, primers, seed_index, seed_len, max_mismatches, stats, seen_in_read)

            # Reverse-complement orientation
            rc = revcomp(seq)
            scan_one_orientation(rc, primers, seed_index, seed_len, max_mismatches, stats, seen_in_read)

            if total_reads % 1_000_000 == 0:
                print(f"Processed {total_reads:,} reads...", file=sys.stderr)

    report_path = out_dir / "primer_hits_report.tsv"
    with report_path.open("w") as out:
        out.write("\t".join([
            "primer_name",
            "amplicon",
            "side",
            "strand",
            "primer_len",
            "primer_seq",
            "reads_processed",
            "reads_with_hit",
            "total_hits",
            "hit_rate_reads",
            "avg_mismatches_per_hit",
        ]) + "\n")

        for i, p in enumerate(primers):
            st = stats[i]
            reads_with_hit = st["reads_with_hit"]
            total_hits = st["total_hits"]
            avg_mm = (st["sum_mismatches"] / total_hits) if total_hits else 0.0
            hit_rate = (reads_with_hit / total_reads) if total_reads else 0.0

            out.write("\t".join([
                p["name"],
                str(p["amplicon"]),
                p["side"],
                p["strand"],
                str(p["len"]),
                p["seq"],
                str(total_reads),
                str(reads_with_hit),
                str(total_hits),
                f"{hit_rate:.6g}",
                f"{avg_mm:.4f}",
            ]) + "\n")

    print(f"Done.")
    print(f"Reads processed: {total_reads:,}")
    print(f"Report written: {report_path.resolve()}")

def parse_args():
    ap = argparse.ArgumentParser(description="Search FASTQ reads for BED primers allowing mismatches, output summary TSV.")
    ap.add_argument("--fastq", default="metagenomic_generator/meta.fastq.gz", help="Input FASTQ(.gz)")
    ap.add_argument("--bed", default="primal_resources_6/primer.bed", help="Primer BED file")
    ap.add_argument("--out", default="metagenomic_generator", help="Output folder")
    ap.add_argument("--max-mismatches", type=int, default=2, help="Allow up to N mismatches (default 2)")
    ap.add_argument("--seed-len", type=int, default=8, help="Seed length for indexing (default 8)")
    ap.add_argument("--limit-reads", type=int, default=0, help="Debug: stop after this many reads (0 = all)")
    return ap.parse_args()

if __name__ == "__main__":
    args = parse_args()
    run(
        fastq_path=Path(args.fastq),
        bed_path=Path(args.bed),
        out_dir=Path(args.out),
        max_mismatches=args.max_mismatches,
        seed_len=args.seed_len,
        limit_reads=args.limit_reads,
    )
