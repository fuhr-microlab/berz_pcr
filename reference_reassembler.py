#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path
from collections import defaultdict


def check_program(program):
    try:
        subprocess.run(
            [program, "--version"],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
            check=False,
        )
    except FileNotFoundError:
        sys.exit(f"ERROR: Required program not found in PATH: {program}")


def run_cmd(cmd, description):
    print(f"\n[RUNNING] {description}")
    print(" ".join(map(str, cmd)))

    result = subprocess.run(cmd)
    if result.returncode != 0:
        sys.exit(f"\nERROR: Step failed: {description}")


def read_fasta(fasta_path):
    records = {}
    current_header = None
    current_seq = []

    with open(fasta_path, "r") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue

            if line.startswith(">"):
                if current_header is not None:
                    records[current_header] = "".join(current_seq).upper()

                current_header = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)

        if current_header is not None:
            records[current_header] = "".join(current_seq).upper()

    return records


def write_wrapped_fasta(records, output_path, width=80):
    with open(output_path, "w") as out:
        for header, seq in records.items():
            out.write(f">{header}\n")
            for i in range(0, len(seq), width):
                out.write(seq[i:i + width] + "\n")


def parse_depth_file(depth_path):
    """
    Returns:
        coverage[ref][1-based-position] = depth
    """
    coverage = defaultdict(dict)

    with open(depth_path, "r") as f:
        for line in f:
            ref, pos, depth = line.rstrip("\n").split("\t")
            coverage[ref][int(pos)] = int(depth)

    return coverage


def make_masked_references(reference_records, coverage, min_depth):
    masked = {}
    stats = []

    for ref, seq in reference_records.items():
        seq_list = list(seq)
        length = len(seq_list)
        covered = 0
        total_depth = 0

        for i in range(length):
            pos = i + 1
            depth = coverage.get(ref, {}).get(pos, 0)
            total_depth += depth

            if depth >= min_depth:
                covered += 1
            else:
                seq_list[i] = "N"

        breadth = covered / length if length > 0 else 0
        mean_depth = total_depth / length if length > 0 else 0

        masked[ref] = "".join(seq_list)

        stats.append({
            "reference": ref,
            "length_bp": length,
            "covered_bp": covered,
            "breadth": breadth,
            "mean_depth": mean_depth,
        })

    return masked, stats


def write_stats(stats, output_path):
    with open(output_path, "w") as out:
        out.write("reference\tlength_bp\tcovered_bp\tbreadth_percent\tmean_depth\n")

        for row in stats:
            out.write(
                f"{row['reference']}\t"
                f"{row['length_bp']}\t"
                f"{row['covered_bp']}\t"
                f"{row['breadth'] * 100:.2f}\t"
                f"{row['mean_depth']:.2f}\n"
            )


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Map filtered FASTQ.gz reads back to pneumo_control.fasta and produce "
            "a masked reference FASTA where uncovered positions are N."
        )
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Filtered input FASTQ.gz",
    )

    parser.add_argument(
        "-r",
        "--reference",
        default="pneumo_control.fasta",
        help="Reference FASTA with serotype sequences [default: pneumo_control.fasta]",
    )

    parser.add_argument(
        "-o",
        "--outprefix",
        default=None,
        help="Output prefix [default: input filename without .fastq.gz]",
    )

    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="Threads [default: 4]",
    )

    parser.add_argument(
        "--preset",
        choices=["sr", "map-ont"],
        default="sr",
        help="minimap2 preset: sr for short reads, map-ont for Nanopore [default: sr]",
    )

    parser.add_argument(
        "--min-mapq",
        type=int,
        default=10,
        help="Minimum mapping quality to keep alignments [default: 10]",
    )

    parser.add_argument(
        "--min-depth",
        type=int,
        default=1,
        help="Minimum read depth for a base to remain unmasked [default: 1]",
    )

    args = parser.parse_args()

    input_fastq = Path(args.input)
    reference = Path(args.reference)

    if not input_fastq.exists():
        sys.exit(f"ERROR: Input FASTQ not found: {input_fastq}")

    if not reference.exists():
        sys.exit(f"ERROR: Reference FASTA not found: {reference}")

    if args.outprefix is None:
        name = input_fastq.name
        for suffix in [".fastq.gz", ".fq.gz", ".gz"]:
            if name.endswith(suffix):
                name = name[: -len(suffix)]
                break
        outprefix = name
    else:
        outprefix = args.outprefix

    bam = Path(f"{outprefix}.mapped.sorted.bam")
    depth_file = Path(f"{outprefix}.depth.tsv")
    masked_fasta = Path(f"{outprefix}.masked_reference_reassembly.fasta")
    stats_file = Path(f"{outprefix}.coverage_summary.tsv")

    check_program("minimap2")
    check_program("samtools")

    map_cmd = [
        "bash",
        "-c",
        (
            f"minimap2 -t {args.threads} -ax {args.preset} {reference} {input_fastq} "
            f"| samtools view -@ {args.threads} -b -F 4 -q {args.min_mapq} "
            f"| samtools sort -@ {args.threads} -o {bam}"
        ),
    ]

    run_cmd(map_cmd, "Mapping reads and creating sorted BAM")

    run_cmd(
        ["samtools", "index", str(bam)],
        "Indexing BAM",
    )

    depth_cmd = [
        "samtools",
        "depth",
        "-aa",
        str(bam),
    ]

    print(f"\n[RUNNING] Calculating per-base depth")
    print(" ".join(depth_cmd) + f" > {depth_file}")

    with open(depth_file, "w") as out:
        result = subprocess.run(depth_cmd, stdout=out)

    if result.returncode != 0:
        sys.exit("\nERROR: samtools depth failed")

    print("\n[RUNNING] Creating masked FASTA")

    reference_records = read_fasta(reference)
    coverage = parse_depth_file(depth_file)

    masked_records, stats = make_masked_references(
        reference_records=reference_records,
        coverage=coverage,
        min_depth=args.min_depth,
    )

    write_wrapped_fasta(masked_records, masked_fasta)
    write_stats(stats, stats_file)

    print("\nDone.")
    print(f"Sorted BAM:              {bam}")
    print(f"BAM index:               {bam}.bai")
    print(f"Depth table:             {depth_file}")
    print(f"Masked reassembly FASTA: {masked_fasta}")
    print(f"Coverage summary:        {stats_file}")


if __name__ == "__main__":
    main()