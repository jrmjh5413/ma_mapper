#!/usr/bin/env python3
"""Measure how each stage scales with the size of a family.

One row per family: copies, alignment columns, and the wall time and peak
memory of parsing, filtering, extracting and mapping. Every family runs in its
own subprocess, so the memory figure is that family's own high-water mark
rather than the maximum over everything measured so far.

    # one family, in this process, printing JSON
    python examples/benchmark_scaling.py --worker --family MER11A \\
        --matrix-dir <dir> --bigwig <track.bw>

    # the sweep, one subprocess per family, writing a TSV
    python examples/benchmark_scaling.py --families MER11A L1PA2 AluJb \\
        --matrix-dir <dir> --bigwig <track.bw> --out scaling.tsv

``--matrix-dir`` holds one directory per family with ``alignment.npy`` (the
parsed alignment, as integer codes) and ``metadata.tsv`` (its coordinates, one
row per fragment, in alignment row order). The alignment is written back out as
FASTA first, so the parse stage is measured on a real alignment at its real
dimensions rather than on a reconstruction.
"""

import argparse
import json
import os
import resource
import subprocess
import sys
import time

import numpy as np
import pandas as pd

BED_COLUMNS = ["chrom", "start", "end", "name", "score", "strand"]


class Stopwatch:
    """Wall time and the process's peak RSS, in seconds and MB."""

    def __init__(self):
        self.timings = {}

    def __call__(self, label):
        self.label = label
        return self

    def __enter__(self):
        self.started = time.perf_counter()
        return self

    def __exit__(self, *exc):
        elapsed = time.perf_counter() - self.started
        peak_mb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024
        self.timings[self.label] = {"seconds": round(elapsed, 3), "peak_rss_mb": round(peak_mb, 1)}
        return False


def write_fasta(alignment, metadata, path):
    """Turn a parsed alignment plus its coordinates back into an aligned FASTA."""
    from ma_mapper import sequence_alignment

    records = sequence_alignment.parse_alignment_to_fasta(alignment)
    names = metadata["name"].drop_duplicates().tolist()
    if len(names) != alignment.shape[0]:
        raise SystemExit(
            f"{alignment.shape[0]} alignment rows but {len(names)} names in the metadata"
        )
    grouped = {name: group for name, group in metadata.groupby("name", sort=False)}
    with open(path, "w") as handle:
        for name, record in zip(names, records):
            group = grouped[name].sort_values("start")
            coords = ",".join(f"{int(a)}-{int(b)}" for a, b in zip(group.start, group.end))
            header = f"{name}::{group.chrom.iloc[0]}:{coords}({group.strand.iloc[0]})"
            handle.write(f">{header}\n{str(record.seq)}\n")
    return path


def run_one(family, matrix_dir, bigwig, workdir, keep_fasta=False, measure_parse=False):
    from ma_mapper import extract_bigwig, mapper

    family_dir = os.path.join(matrix_dir, family)
    alignment_npy = os.path.join(family_dir, "alignment.npy")
    metadata_tsv = os.path.join(family_dir, "metadata.tsv")
    for path in (alignment_npy, metadata_tsv):
        if not os.path.isfile(path):
            raise SystemExit(f"missing {path}")

    stored = np.load(alignment_npy).astype(np.uint8)
    metadata = pd.read_csv(metadata_tsv, sep="\t")
    metadata = metadata[BED_COLUMNS]

    watch = Stopwatch()
    # Parsing is measured on a real alignment at its real dimensions, which
    # means writing the stored matrix back out as FASTA first. That write is
    # scaffolding, not a cost of the tool, and on the largest families it takes
    # longer than everything else put together -- so it is opt-in, and the
    # default path measures filtering, extraction and mapping straight from the
    # stored arrays.
    if measure_parse:
        fasta = os.path.join(workdir, f"{family}.aligned")
        with watch("write_fasta_scaffolding"):
            write_fasta(stored, metadata, fasta)
        with watch("parse_alignment"):
            alignment = mapper.parse_alignment(fasta)
        with watch("extract_coordinates"):
            coordinates = mapper.extract_coordinate_from_alignment(fasta)
        if not keep_fasta:
            os.remove(fasta)
    else:
        alignment = stored
        coordinates = metadata

    with watch("create_filter"):
        row_filter, col_filter = mapper.create_filter(
            alignment, col_threshold=0.05, col_content_threshold=0.05, row_threshold=0.10
        )
    with watch("extract_bigwig"):
        extracted = extract_bigwig.bigwig_io(coordinates, bigwig)
    with watch("map_data"):
        matrix = mapper.map_data(extracted, alignment, filter=[row_filter, col_filter])

    return {
        "family": family,
        "measured_parse": bool(measure_parse),
        "copies": int(alignment.shape[0]),
        "columns": int(alignment.shape[1]),
        "cells": int(alignment.shape[0]) * int(alignment.shape[1]),
        "copies_kept": int(row_filter.sum()),
        "columns_kept": int(col_filter.sum()),
        "matrix_mb": round(matrix.nbytes / 1024 / 1024, 1),
        "stages": watch.timings,
        "total_seconds": round(sum(s["seconds"] for s in watch.timings.values()), 2),
        "peak_rss_mb": max(s["peak_rss_mb"] for s in watch.timings.values()),
    }


def main(argv=None):
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument("--family", help="one family (with --worker)")
    parser.add_argument("--families", nargs="+", help="the families to sweep")
    parser.add_argument("--matrix-dir", required=True, help="one directory per family")
    parser.add_argument("--bigwig", required=True, help="a bigWig to extract")
    parser.add_argument("--workdir", default=".", help="where the temporary FASTAs go")
    parser.add_argument("--out", help="write the table here (TSV)")
    parser.add_argument("--worker", action="store_true", help="measure one family and print JSON")
    parser.add_argument("--keep-fasta", action="store_true")
    parser.add_argument(
        "--measure-parse",
        action="store_true",
        help="also measure parse_alignment, which needs the alignment written out as "
        "FASTA first (slow on the largest families, and the write is scaffolding)",
    )
    args = parser.parse_args(argv)

    os.makedirs(args.workdir, exist_ok=True)

    if args.worker:
        if not args.family:
            raise SystemExit("--worker needs --family")
        print(json.dumps(run_one(args.family, args.matrix_dir, args.bigwig,
                                 args.workdir, args.keep_fasta, args.measure_parse)))
        return 0

    if not args.families:
        raise SystemExit("give --families, or --worker with --family")

    rows = []
    for family in args.families:
        command = [
            sys.executable, os.path.abspath(__file__), "--worker",
            "--family", family,
            "--matrix-dir", args.matrix_dir,
            "--bigwig", args.bigwig,
            "--workdir", args.workdir,
        ]
        if args.measure_parse:
            command.append("--measure-parse")
        started = time.perf_counter()
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # wait4 gives this child's own rusage, not the maximum over all children.
        _, status, usage = os.wait4(process.pid, 0)
        stdout = process.stdout.read().decode()
        stderr = process.stderr.read().decode()
        if status != 0:
            print(f"{family}: FAILED\n{stderr[-800:]}", file=sys.stderr)
            continue
        record = json.loads(stdout.strip().splitlines()[-1])
        record["child_peak_rss_mb"] = round(usage.ru_maxrss / 1024, 1)
        record["child_wall_seconds"] = round(time.perf_counter() - started, 2)
        record["child_cpu_seconds"] = round(usage.ru_utime + usage.ru_stime, 2)
        rows.append(record)
        print(
            f"{family:16s} {record['copies']:>7,} copies x {record['columns']:>6,} columns  "
            f"{record['child_wall_seconds']:>7.1f} s  {record['child_peak_rss_mb']:>8.1f} MB",
            flush=True,
        )

    if not rows:
        return 1

    flat = []
    for record in rows:
        line = {k: v for k, v in record.items() if k != "stages"}
        for stage, values in record["stages"].items():
            line[f"{stage}_s"] = values["seconds"]
        flat.append(line)
    table = pd.DataFrame(flat).sort_values("copies")
    if args.out:
        table.to_csv(args.out, sep="\t", index=False)
        print(f"\nwrote {args.out}")
    else:
        print(table.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())
