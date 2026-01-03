#!/usr/bin/env python3
from pathlib import Path
import random
import argparse

FASTA_EXTENSIONS = {".fa", ".fasta", ".fna"}

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
                seq.append(line)
        if header is not None:
            yield header, "".join(seq).upper()

def write_single_fasta(path: Path, header: str, seq: str):
    with path.open("w") as f:
        f.write(f">{header}\n")
        for i in range(0, len(seq), 80):
            f.write(seq[i:i+80] + "\n")

def run(
    input_dir: Path,
    output_dir: Path,
    chunk_size: int = 1000,
    target_reads: int = 1000,
    seed=None,
    avoid_duplicates: bool = True,
    max_attempts_multiplier: int = 50,
    weighted: bool = True,
):

    output_dir.mkdir(exist_ok=True)

    rng = random.Random(seed) if seed is not None else random.Random()

    fasta_files = sorted(
        [p for p in input_dir.iterdir() if p.suffix.lower() in FASTA_EXTENSIONS],
        key=lambda p: p.name.lower(),
    )

    # if you want non weighted genomes use lines below, dont forget other line as well
    
    templates = []
    weights = []

    for f in fasta_files:
        for rec_i, (_hdr, seq) in enumerate(read_fasta(f), start=1):
            if len(seq) >= chunk_size:
                template_id = f"{f.name}::rec{rec_i}"
                templates.append((template_id, seq))
                if weighted:
                    weights.append(len(seq) - chunk_size + 1)

    if not templates:
        raise SystemExit(f"No sequences of length >= {chunk_size} found in {input_dir}")


    max_possible = sum(len(seq) - chunk_size + 1 for _, seq in templates)

    written = 0
    seg_idx = 1
    seen = set()

    max_attempts = target_reads * max_attempts_multiplier
    attempts = 0

    while written < target_reads and attempts < max_attempts:
        attempts += 1
        # if you want non weighted genomes use
        if weighted:
            template_id, seq = rng.choices(templates, weights=weights, k=1)[0]
        else:
            template_id, seq = rng.choice(templates)

        max_start = len(seq) - chunk_size
        start = rng.randint(0, max_start)
        end = start + chunk_size

        if avoid_duplicates:
            key = (template_id, start)
            if key in seen:
                continue
            seen.add(key)

        chunk = seq[start:end]
        seg_name = f"segment_{seg_idx:06d}"
        out_path = output_dir / f"{seg_name}.fasta"
        write_single_fasta(out_path, seg_name, chunk)

        seg_idx += 1
        written += 1

    print("Done.")
    print(f"Input: {input_dir.resolve()}")
    print(f"Output: {output_dir.resolve()}")
    print(f"Chunk size: {chunk_size}")
    print(f"Wrote: {written} reads (target {target_reads})")
    print(f"Attempts: {attempts}")
    if avoid_duplicates:
        print(f"Unique windows used: {len(seen)}")
    print(f"Approx max possible unique windows (upper bound): {max_possible}")
    if written < target_reads:
        print("WARNING: Could not reach target read count (not enough unique windows or too strict settings).")

def parse_args():
    ap = argparse.ArgumentParser(description="Randomly shear cps loci into fixed-length FASTA reads.")
    ap.add_argument("--input", default="resources", help="Input folder containing FASTA loci.")
    ap.add_argument("--output", default="holocene", help="Output folder for sheared reads.")
    ap.add_argument("--chunk", type=int, default=1000, help="Chunk/read length in bp.")
    ap.add_argument("--target", type=int, default=1000, help="Target number of reads to write.")
    ap.add_argument("--seed", default="1337", help="Random seed, or 'none' for non-reproducible.")
    ap.add_argument("--avoid-dups", type=int, default=1, help="1 to avoid duplicate windows, 0 to allow.")
    ap.add_argument("--max-attempts-mult", type=int, default=50, help="Safety multiplier for sampling attempts.")
    ap.add_argument("--weighted", type=int, default=1,
                help="1 = weight by possible window count (length-chunk+1), 0 = uniform per record.")

    
    return ap.parse_args()

if __name__ == "__main__":
    args = parse_args()
    seed = None if str(args.seed).strip().lower() in ("none", "") else int(args.seed)
    run(
        input_dir=Path(args.input),
        output_dir=Path(args.output),
        chunk_size=args.chunk,
        target_reads=args.target,
        seed=seed,
        avoid_duplicates=bool(args.avoid_dups),
        max_attempts_multiplier=args.max_attempts_mult,
        weighted=bool(args.weighted),

    )
