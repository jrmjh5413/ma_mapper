"""Command line for ma_mapper.

Five steps, each its own subcommand, so a whole overlay can be built without
writing any Python:

``ma-mapper sequence``
    coordinates -> one FASTA record per copy, with the coordinates in the header
``ma-mapper align``
    that FASTA -> a multiple alignment (MAFFT)
``ma-mapper overlay``
    alignment + dataset -> a copies x columns matrix, with its row names
``ma-mapper plot``
    matrix -> a heatmap or a per-column profile
``ma-mapper index-maf``
    build a MAF index once, before any parallel extraction

``ma-mapper info`` prints the versions of everything the package depends on,
which is the first thing to paste into a bug report.
"""

import argparse
import os
import sys

from . import __version__, configure_logging, logger


def _add_filter_arguments(parser):
    group = parser.add_argument_group(
        "alignment filtering",
        "Which columns and copies survive into the matrix. The defaults are the "
        "package defaults, which are stricter than most analyses want: see "
        "docs/package_docs.md.",
    )
    group.add_argument(
        "--filter",
        choices=("one_pass", "iterative", "none"),
        default="one_pass",
        help="one_pass applies columns -> rows -> columns once (default); "
        "iterative repeats until the masks stop changing; none keeps everything.",
    )
    group.add_argument(
        "--col-threshold",
        type=float,
        default=0.50,
        help="keep a column non-gap in more than this fraction of the copies that "
        "span it (default: 0.50)",
    )
    group.add_argument(
        "--col-content-threshold",
        type=float,
        default=0.10,
        help="and non-gap in at least this fraction of all copies (default: 0.10)",
    )
    group.add_argument(
        "--row-threshold",
        type=float,
        default=0.50,
        help="keep a copy whose non-gap bases are at least this dense within its own "
        "span (default: 0.50). This is a density, not a length: it removes most "
        "copies of an ancient, deletion-riddled family, so stratified analyses "
        "usually want 0.10 or none.",
    )
    group.add_argument(
        "--legacy-second-pass",
        action="store_true",
        help="reproduce the pre-0.9.0 column arithmetic, for rebuilding old matrices "
        "bit-for-bit",
    )


def _filter_kwargs(args):
    kwargs = {
        "col_threshold": args.col_threshold,
        "col_content_threshold": args.col_content_threshold,
        "row_threshold": args.row_threshold,
        "legacy_second_pass": args.legacy_second_pass,
    }
    return kwargs


def cmd_sequence(args):
    from . import sequence_alignment

    if args.repeatmasker:
        if not args.subfamily:
            raise SystemExit("--repeatmasker needs --subfamily")
        coordinates = sequence_alignment.extract_coord_from_repeatmasker_table(
            subfamily=args.subfamily,
            repeatmasker_table=args.repeatmasker,
            internal_id_table=args.internal_id_table,
        )
        logger.info("%d copies of %s", len(coordinates), args.subfamily)
    else:
        if not args.coordinates:
            raise SystemExit("give either --coordinates or --repeatmasker with --subfamily")
        coordinates = args.coordinates

    sequence_alignment.sequence_io(
        coordinate_table=coordinates,
        source_fasta=args.genome,
        save_to_file=args.output,
        generate_new_id=args.generate_new_id,
    )
    print(args.output)


def cmd_align(args):
    from . import sequence_alignment

    output = sequence_alignment.mafft_align(
        args.input,
        nthread=args.threads,
        output_filepath=args.output,
        mafft_arg=args.mafft_arg,
    )
    print(output)


def cmd_overlay(args):
    import numpy as np

    from . import mapper

    filter_mode = None if args.filter == "none" else args.filter
    matrix, coordinate_table, alignment_matrix = mapper.map_and_overlay(
        alignment=args.alignment,
        data_file=args.data,
        data_format=args.format,
        filter=filter_mode is not None,
        filter_mode=filter_mode or "one_pass",
        coordinate_file=args.coordinates,
        return_details=True,
        **_filter_kwargs(args),
    )
    out = args.output
    stem = out[:-4] if out.endswith(".npy") else out
    np.save(out, matrix)
    # Row names and the alignment itself, so the matrix stays interpretable:
    # which copy each row is, and which of its cells are gaps.
    rows_path = stem + ".rows.tsv"
    coordinate_table.to_csv(rows_path, sep="\t", index=False)
    alignment_path = stem + ".alignment.npy"
    np.save(alignment_path, alignment_matrix)
    logger.info("matrix %s -> %s", matrix.shape, out)
    print(out)
    print(rows_path)
    print(alignment_path)


def cmd_plot(args):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    from . import mapper, plots

    matrix = np.load(args.matrix, allow_pickle=True)
    if args.kind == "heatmap":
        plots.plot_heatmap(
            data=[matrix],
            save_to_file=args.output,
            image_res=args.dpi,
            plot_title=args.title,
        )
    else:
        # The profile divides by the number of copies with a base in each
        # column, so it needs the alignment *as it was filtered for this
        # matrix*. `overlay` saves exactly that next to the matrix.
        stem = args.matrix[:-4] if args.matrix.endswith(".npy") else args.matrix
        sidecar = stem + ".alignment.npy"
        if args.alignment_matrix:
            alignment_matrix = np.load(args.alignment_matrix)
        elif os.path.isfile(sidecar):
            alignment_matrix = np.load(sidecar)
        elif args.alignment:
            alignment_matrix = mapper.parse_alignment(args.alignment)
        else:
            raise SystemExit(
                "--kind profile needs the alignment: pass --alignment-matrix, or keep "
                f"the {os.path.basename(sidecar)} written by `overlay` beside the matrix"
            )
        if alignment_matrix.shape != matrix.shape:
            raise SystemExit(
                f"the alignment is {alignment_matrix.shape} but the matrix is "
                f"{matrix.shape}: they were built with different filters"
            )
        profile = mapper.normalise(alignment_matrix, matrix, method=args.method)
        fig, ax = plt.subplots(figsize=(8, 2))
        ax.plot(np.arange(profile.size), profile, linewidth=0.8)
        ax.set_xlabel("alignment column")
        ax.set_ylabel(args.method)
        if args.title:
            ax.set_title(args.title)
        ax.margins(x=0)
        fig.tight_layout()
        fig.savefig(args.output, dpi=args.dpi)
        plt.close(fig)
    print(args.output)


def cmd_index_maf(args):
    from . import extract_maf

    path = extract_maf.build_maf_index(
        args.maf, target_species=args.target_species, chrom=args.chrom
    )
    print(path)


def cmd_info(args):
    import importlib

    print(f"ma_mapper {__version__}")
    print(f"python {sys.version.split()[0]} ({sys.executable})")
    for name in (
        "numpy",
        "pandas",
        "Bio",
        "pysam",
        "pyBigWig",
        "pybedtools",
        "cyvcf2",
        "h5py",
        "matplotlib",
        "logomaker",
        "compress_pickle",
        "scipy",
    ):
        try:
            module = importlib.import_module(name)
            print(f"  {name:16s} {getattr(module, '__version__', 'unknown')}")
        except Exception as exc:
            print(f"  {name:16s} MISSING ({exc.__class__.__name__})")
    import shutil

    for tool in ("mafft", "hmmalign", "bedtools", "tabix"):
        path = shutil.which(tool)
        print(f"  {tool:16s} {path or 'not on PATH'}")


def build_parser():
    parser = argparse.ArgumentParser(
        prog="ma-mapper",
        description="Map genome-wide data onto a multiple alignment of transposable elements.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--version", action="version", version=f"ma_mapper {__version__}")
    parser.add_argument(
        "--verbose", action="store_true", help="log progress to stderr"
    )
    parser.add_argument("--logfile", help="also append log records to this file")
    subparsers = parser.add_subparsers(dest="command", required=True)

    sequence = subparsers.add_parser(
        "sequence",
        help="fetch the genomic sequence of every copy",
        description="Write one FASTA record per copy, with the coordinates in the header.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    sequence.add_argument("--coordinates", help="BED-like table of copies (6+ columns)")
    sequence.add_argument("--repeatmasker", help="RepeatMasker table, instead of --coordinates")
    sequence.add_argument("--subfamily", nargs="+", help="subfamily name(s) to take from it")
    sequence.add_argument("--internal-id-table", dest="internal_id_table", help="optional id mapping")
    sequence.add_argument("--genome", required=True, help="genome FASTA (indexed if possible)")
    sequence.add_argument("--output", required=True, help="output FASTA")
    sequence.add_argument(
        "--generate-new-id",
        dest="generate_new_id",
        action="store_true",
        help="number the copies instead of grouping rows by name (does not reassemble "
        "fragmented elements)",
    )
    sequence.set_defaults(func=cmd_sequence)

    align = subparsers.add_parser(
        "align",
        help="align a FASTA with MAFFT",
        description="MAFFT's default algorithm. Above ~50,000 copies see "
        "docs/aligning_large_families.md instead.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    align.add_argument("--input", required=True, help="unaligned FASTA")
    align.add_argument("--output", help="aligned FASTA (default: <input>.aligned)")
    align.add_argument("--threads", type=int, help="MAFFT --thread")
    align.add_argument(
        "--mafft-arg",
        dest="mafft_arg",
        default="--treeout --reorder",
        help="extra MAFFT flags, as one quoted string",
    )
    align.set_defaults(func=cmd_align)

    overlay = subparsers.add_parser(
        "overlay",
        help="map a dataset onto an alignment",
        description="Write the copies x columns matrix as .npy, plus a .rows.tsv of "
        "the copies that survived filtering, in matrix row order.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    overlay.add_argument("--alignment", required=True, help="aligned FASTA")
    overlay.add_argument("--data", required=True, help="BAM, bigWig, BED, MAF or VCF")
    overlay.add_argument(
        "--format",
        required=True,
        help="bigwig, bed, maf, vcf, or a BAM readout: read_sum, read_max, read_min, "
        "read_forward, read_reverse, normal_max, normal_min, normal_forward, "
        "normal_reverse",
    )
    overlay.add_argument(
        "--coordinates",
        help="use these coordinates instead of the alignment headers, matched by name",
    )
    overlay.add_argument("--output", required=True, help="output .npy")
    _add_filter_arguments(overlay)
    overlay.set_defaults(func=cmd_overlay)

    plot = subparsers.add_parser(
        "plot",
        help="draw a matrix",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    plot.add_argument("--matrix", required=True, help=".npy written by `overlay`")
    plot.add_argument(
        "--alignment", help="aligned FASTA, for --kind profile if no sidecar is present"
    )
    plot.add_argument(
        "--alignment-matrix",
        dest="alignment_matrix",
        help="the .alignment.npy sidecar written by `overlay` (found automatically if "
        "it sits beside the matrix)",
    )
    plot.add_argument("--output", required=True, help="output PNG")
    plot.add_argument("--kind", choices=("heatmap", "profile"), default="heatmap")
    plot.add_argument(
        "--method",
        choices=("average", "perc_coverage", "median"),
        default="average",
        help="how --kind profile collapses each column",
    )
    plot.add_argument("--dpi", type=int, default=600)
    plot.add_argument("--title")
    plot.set_defaults(func=cmd_plot)

    index_maf = subparsers.add_parser(
        "index-maf",
        help="build a MAF index up front",
        description="Biopython builds this index on first use. Doing it once here "
        "avoids every worker process building the same one at the same time.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    index_maf.add_argument("--maf", required=True)
    index_maf.add_argument("--target-species", dest="target_species", default="Homo_sapiens")
    index_maf.add_argument("--chrom", help="default: the first dot-separated part of the filename")
    index_maf.set_defaults(func=cmd_index_maf)

    info = subparsers.add_parser("info", help="print versions of dependencies and external tools")
    info.set_defaults(func=cmd_info)

    return parser


def main(argv=None):
    parser = build_parser()
    args = parser.parse_args(argv)
    if args.verbose or args.logfile:
        configure_logging(logfile=args.logfile, stream=args.verbose)
    try:
        args.func(args)
    except (FileNotFoundError, ValueError, TypeError, KeyError) as exc:
        parser.exit(2, f"ma-mapper {args.command}: {exc}\n")
    return 0


if __name__ == "__main__":  # pragma: no cover
    sys.exit(main())
