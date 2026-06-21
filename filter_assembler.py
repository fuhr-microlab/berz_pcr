#!/usr/bin/env python3

import argparse
import subprocess
import sys
from pathlib import Path


def check_program(program):
    """Check that an external command exists."""
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


def main():
    parser = argparse.ArgumentParser(
        description=(
            "Filter a gzipped metagenomic FASTQ against pneumo_control.fasta "
            "and output matching reads as a new FASTQ.gz."
        )
    )

    parser.add_argument(
        "-i",
        "--input",
        required=True,
        help="Input metagenomic FASTQ.gz file",
    )

    parser.add_argument(
        "-r",
        "--reference",
        default="pneumo_control.fasta",
        help="Reference FASTA containing pneumococcal/serotype sequences "
             "[default: pneumo_control.fasta]",
    )

    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output filtered FASTQ.gz file [default: <input>.pneumo_filtered.fastq.gz]",
    )

    parser.add_argument(
        "-t",
        "--threads",
        type=int,
        default=4,
        help="Number of threads to use [default: 4]",
    )

    parser.add_argument(
        "--preset",
        default="map-ont",
        choices=["map-ont", "sr"],
        help="minimap2 preset: map-ont for Nanopore/long reads, sr for short reads [default: map-ont]",
    )

    parser.add_argument(
        "-q",
        "--mapq",
        type=int,
        default=10,
        help="Minimum mapping quality to keep read [default: 10]",
    )

    args = parser.parse_args()

    input_fastq = Path(args.input)
    reference = Path(args.reference)

    if not input_fastq.exists():
        sys.exit(f"ERROR: Input file not found: {input_fastq}")

    if not reference.exists():
        sys.exit(f"ERROR: Reference file not found: {reference}")

    if args.output is None:
        name = input_fastq.name
        if name.endswith(".fastq.gz"):
            name = name.replace(".fastq.gz", "")
        elif name.endswith(".fq.gz"):
            name = name.replace(".fq.gz", "")
        elif name.endswith(".gz"):
            name = name.replace(".gz", "")
        output_fastq = Path(f"{name}.pneumo_filtered.fastq.gz")
    else:
        output_fastq = Path(args.output)

    bam_file = Path("tmp_pneumo_filtered.bam")
    sorted_bam = Path("tmp_pneumo_filtered.sorted.bam")

    check_program("minimap2")
    check_program("samtools")

    # Map reads to pneumo reference and keep mapped reads above MAPQ threshold.
    # -F 4 removes unmapped reads.
    map_cmd = [
        "bash",
        "-c",
        (
            f"minimap2 -t {args.threads} -ax {args.preset} {reference} {input_fastq} "
            f"| samtools view -@ {args.threads} -b -F 4 -q {args.mapq} -o {bam_file}"
        ),
    ]

    run_cmd(map_cmd, "Mapping and filtering reads")

    sort_cmd = [
        "samtools",
        "sort",
        "-@",
        str(args.threads),
        "-o",
        str(sorted_bam),
        str(bam_file),
    ]

    run_cmd(sort_cmd, "Sorting filtered BAM")

    fastq_cmd = [
        "bash",
        "-c",
        f"samtools fastq -@ {args.threads} {sorted_bam} | gzip > {output_fastq}",
    ]

    run_cmd(fastq_cmd, "Writing filtered FASTQ.gz")

    # Cleanup temporary files
    for tmp in [bam_file, sorted_bam]:
        if tmp.exists():
            tmp.unlink()

    print("\nDone.")
    print(f"Filtered reads written to: {output_fastq}")


if __name__ == "__main__":
    main()