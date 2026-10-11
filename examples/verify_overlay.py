#!/usr/bin/env python3
"""Check an overlay against an independent reading of the same bigWig.

Builds the overlay the way the package does, then, for a sample of copies,
reads the track directly with pyBigWig, reverses it by hand where the copy is
on the minus strand, and compares base by base. Any disagreement is a bug in
the package or a mismatch between the alignment and the coordinates, and either
way you want to know before the figure.

    python examples/verify_overlay.py \\
        --coordinates copies.bed --genome hg38.fa \\
        --bigwig hg38.phyloP447way.bw --sample 10

The coordinates file is a six-column BED (chrom, start, end, name, score,
strand), one row per fragment, as many rows per copy as it has fragments.
"""

import argparse
import sys

import numpy as np
import pandas as pd

from ma_mapper import mapper, sequence_alignment


def build_alignment(coordinates, genome, path):
    """Fetch the sequences and pad them to a common width.

    Padding is not a real alignment, but it is a valid one, and it keeps this
    check independent of which aligner you use.
    """
    records = sequence_alignment.sequence_io(coordinates, genome)
    width = max(len(record.seq) for record in records)
    with open(path, "w") as handle:
        for record in records:
            handle.write(f">{record.id}\n{str(record.seq).ljust(width, '-')}\n")
    return path, len(records), width


def independent_read(bigwig, group, strand):
    """The same values, read straight out of the file."""
    pieces = []
    for row in group.sort_values("start").itertuples():
        values = bigwig.values(row.chrom, int(row.start), int(row.end), numpy=True)
        pieces.append(values[::-1] if strand == "-" else values)
    if strand == "-":
        pieces.reverse()
    return np.concatenate(pieces)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--coordinates", required=True, help="six-column BED of the copies")
    parser.add_argument("--genome", required=True, help="genome FASTA, indexed if possible")
    parser.add_argument("--bigwig", required=True, help="the track to check")
    parser.add_argument("--alignment", help="write the padded alignment here (default: alongside --coordinates)")
    parser.add_argument("--sample", type=int, default=10, help="how many copies to verify")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args(argv)

    import pyBigWig

    columns = ["chrom", "start", "end", "name", "score", "strand"]
    coordinates = pd.read_csv(args.coordinates, sep="\t", header=None, usecols=range(6))
    coordinates.columns = columns
    print(f"{len(coordinates)} rows, {coordinates.name.nunique()} copies, "
          f"{(coordinates.strand == '-').sum()} minus-strand rows")

    alignment_path = args.alignment or (args.coordinates + ".padded.fasta")
    alignment_path, n_records, width = build_alignment(coordinates, args.genome, alignment_path)
    print(f"alignment: {n_records} records x {width} columns -> {alignment_path}")

    matrix, surviving, alignment = mapper.map_and_overlay(
        alignment=alignment_path,
        data_file=args.bigwig,
        data_format="bigwig",
        filter=False,
        return_details=True,
    )
    print(f"overlay {matrix.shape}, {int((alignment > 0).sum())} non-gap cells")

    row_of = {name: index for index, name in enumerate(surviving["name"].drop_duplicates())}
    names = list(row_of)
    rng = np.random.default_rng(args.seed)
    picks = [names[i] for i in rng.choice(len(names), min(args.sample, len(names)), replace=False)]
    # make sure the awkward cases are in the sample
    fragmented = coordinates.name.value_counts()
    if (fragmented > 1).any():
        picks.append(fragmented[fragmented > 1].index[0])
    minus = coordinates[coordinates.strand == "-"]
    if len(minus):
        picks.append(minus.name.iloc[0])

    bigwig = pyBigWig.open(args.bigwig)
    mismatched = []
    try:
        for name in dict.fromkeys(picks):
            group = coordinates[coordinates.name == name]
            strand = group.strand.iloc[0]
            expected = independent_read(bigwig, group, strand)
            row = row_of[name]
            placed = matrix[row][alignment[row] > 0].astype(float)
            same_length = len(expected) == len(placed)
            agrees = same_length and np.allclose(
                np.nan_to_num(expected, nan=0.0), np.nan_to_num(placed, nan=0.0), atol=1e-4
            )
            print(f"  {name:20s} {strand} fragments={len(group):2d} bases={len(expected):5d} "
                  f"agrees={agrees}")
            if not agrees:
                mismatched.append(name)
    finally:
        bigwig.close()

    if mismatched:
        print(f"\n{len(mismatched)} copies disagree with the independent reading: {mismatched}")
        return 1
    print(f"\nall {len(dict.fromkeys(picks))} sampled copies agree base by base")
    return 0


if __name__ == "__main__":
    sys.exit(main())
