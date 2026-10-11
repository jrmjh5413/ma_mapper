#!/usr/bin/env python3
"""Build a small example from nothing and run the whole chain on it.

No downloads and no reference genome: the script writes a 400 bp chromosome, a
four-copy family on it, and a bigWig whose value at position ``p`` is ``p``,
then maps that track onto an alignment of the copies. Because the expected
answer is arithmetic, the script checks itself and says so.

    python examples/toy_end_to_end.py /tmp/ma_mapper_example

Run it after installing to confirm the installation works end to end, or read
it as the shortest complete example of the API.
"""

import os
import sys

import numpy as np
import pandas as pd

from ma_mapper import configure_logging, mapper, sequence_alignment

CHROM = "chrT"
LENGTH = 400

# name -> fragments, strand. copyC is split in two, copyB and copyD are on the
# minus strand, which is where orientation mistakes show up.
COPIES = {
    "copyA": ([(100, 120)], "+"),
    "copyB": ([(200, 220)], "-"),
    "copyC": ([(300, 310), (320, 330)], "+"),
    "copyD": ([(50, 65)], "-"),
}


def write_genome(path):
    sequence = ("ACGT" * (LENGTH // 4))[:LENGTH]
    with open(path, "w") as handle:
        handle.write(f">{CHROM}\n")
        for i in range(0, len(sequence), 60):
            handle.write(sequence[i : i + 60] + "\n")
    import pysam

    pysam.faidx(path)
    return path


def write_bigwig(path):
    """A track whose value at position p is p, so orientation is checkable."""
    import pyBigWig

    handle = pyBigWig.open(path, "w")
    handle.addHeader([(CHROM, LENGTH)])
    starts = list(range(LENGTH))
    handle.addEntries(
        [CHROM] * LENGTH,
        starts,
        ends=[p + 1 for p in starts],
        values=[float(p) for p in starts],
    )
    handle.close()
    return path


def coordinate_table():
    rows = []
    for name, (fragments, strand) in COPIES.items():
        for start, end in fragments:
            rows.append([CHROM, start, end, name, 10, strand])
    return pd.DataFrame(rows)


def expected(name):
    """What the track should read for this copy, in the element's orientation."""
    fragments, strand = COPIES[name]
    values = []
    for start, end in fragments:
        piece = list(range(start, end))
        values.append(piece[::-1] if strand == "-" else piece)
    if strand == "-":
        values.reverse()
    return np.array([v for piece in values for v in piece], dtype=float)


def main(workdir):
    os.makedirs(workdir, exist_ok=True)
    configure_logging()

    genome = write_genome(os.path.join(workdir, "genome.fa"))
    bigwig = write_bigwig(os.path.join(workdir, "positions.bw"))
    coordinates = coordinate_table()

    # 1. sequence. The header of each record carries its coordinates, which is
    #    how the overlay later knows where each row came from.
    records = sequence_alignment.sequence_io(coordinates, genome)
    print(f"fetched {len(records)} copies, lengths {[len(r.seq) for r in records]}")

    # 2. "align" them. A real run would call sequence_alignment.mafft_align
    #    here; padding to a common width keeps this script free of external
    #    binaries while producing a valid alignment.
    width = max(len(r.seq) for r in records)
    alignment_path = os.path.join(workdir, "family.aligned")
    with open(alignment_path, "w") as handle:
        for record in records:
            handle.write(f">{record.id}\n{str(record.seq).ljust(width, '-')}\n")

    # 3. map the track onto the alignment.
    matrix, surviving, alignment = mapper.map_and_overlay(
        alignment=alignment_path,
        data_file=bigwig,
        data_format="bigwig",
        filter=False,
        return_details=True,
    )
    print(f"overlay matrix {matrix.shape}, {int((alignment > 0).sum())} non-gap cells")

    # 4. check it: every value must sit on its own base.
    failures = []
    for row, name in enumerate(surviving["name"].drop_duplicates()):
        placed = matrix[row][alignment[row] > 0]
        if not np.allclose(placed, expected(name)):
            failures.append(name)
    print("every value in the column of its own base:", not failures)
    if failures:
        print("  MISMATCH for:", failures)
        return 1

    # 5. one number per column, over the copies that have a base there.
    profile = mapper.normalise(alignment, matrix, method="average")
    print(f"per-column mean over {matrix.shape[1]} columns, "
          f"first five: {np.round(profile[:5], 1)}")

    out = os.path.join(workdir, "matrix.npy")
    np.save(out, matrix)
    surviving.to_csv(os.path.join(workdir, "matrix.rows.tsv"), sep="\t", index=False)
    print(f"wrote {out} and matrix.rows.tsv in {workdir}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1] if len(sys.argv) > 1 else "ma_mapper_example"))
