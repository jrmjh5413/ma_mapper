"""Read a bigWig signal track per base for each copy of a family.

Bases the track does not cover come back as ``nan``, not 0.0 -- "not measured"
and "measured as zero" are different statements, and the distinction has to
survive into the overlay matrix. ``mapper.normalise`` is the place that decides
how to treat them.
"""

import os

import numpy as np
import pyBigWig

from . import logger
from ._io_common import group_coordinates, read_coordinate_table, save_or_return


def _open_bigwig(bigwig_file):
    """Return ``(handle, needs_closing)`` for a path or an open pyBigWig handle."""
    if isinstance(bigwig_file, str):
        return pyBigWig.open(bigwig_file), True
    if hasattr(bigwig_file, "values") and hasattr(bigwig_file, "chroms"):
        return bigwig_file, False
    raise TypeError(
        "bigwig_file must be a path or an open pyBigWig handle, got "
        f"{type(bigwig_file).__name__}"
    )


def extract_bigwig(bigwig_file, chrom, start_list, end_list, strand):
    """Per-base signal for one copy, in element orientation.

    Parameters
    ----------
    bigwig_file : str or pyBigWig handle
        A path is opened and closed here; an open handle is left alone.
    chrom : str
        Chromosome name as spelled in the bigWig. A name the file does not
        carry raises ``KeyError`` naming the mismatch, rather than the
        bare RuntimeError pyBigWig raises.
    start_list, end_list : list of int
        Fragment coordinates, 0-based half-open (BED), ascending.
    strand : {'+', '-'}
        ``'-'`` reverses each fragment and the fragment order.

    Returns
    -------
    numpy.ndarray
        float64, length ``sum(end - start)``; ``nan`` where the track has no data.
    """
    bigwig, close_after = _open_bigwig(bigwig_file)
    try:
        if chrom not in bigwig.chroms():
            example = next(iter(bigwig.chroms()), "?")
            raise KeyError(
                f"chromosome '{chrom}' is not in the bigWig "
                f"(it uses names like '{example}'): check for a chr-prefix mismatch"
            )
        chrom_length = bigwig.chroms()[chrom]

        bigwig_arrays = []
        for i in range(len(start_list)):
            start = int(start_list[i])
            end = int(end_list[i])
            if start < 0 or end > chrom_length:
                raise ValueError(
                    f"fragment {chrom}:{start}-{end} falls outside the bigWig's "
                    f"{chrom} (length {chrom_length})"
                )
            bigwig_array = bigwig.values(chrom, start, end, numpy=True)
            if strand == "-":
                bigwig_array = np.flip(bigwig_array)
            bigwig_arrays.append(bigwig_array)
        if strand == "-":
            bigwig_arrays.reverse()
        return np.concatenate(bigwig_arrays)
    finally:
        if close_after:
            bigwig.close()


def bigwig_io(coordinate_table, bigwig, save_to_file=False, generate_new_id=False):
    """Extract a bigWig track for every copy in a coordinate table.

    Parameters
    ----------
    coordinate_table : pandas.DataFrame or str
        BED-like table, one or more rows per copy.
    bigwig : str
        Path to the bigWig file.
    save_to_file : bool or str
        ``True`` writes ``bigwig_out.p`` beside the coordinate table; a string
        is used as the output path. Before 0.9.0 passing a path wrote nothing
        and reported nothing.
    generate_new_id : bool
        Treat every row as its own copy instead of grouping by column 3.

    Returns
    -------
    list of numpy.ndarray
        One array per copy, or the output path if the result was saved.
    """
    coordinate_local = read_coordinate_table(coordinate_table)
    meta_id, chrom_list, start_list, end_list, strand_list = group_coordinates(
        coordinate_local, generate_new_id=generate_new_id
    )

    if not os.path.isfile(bigwig):
        raise FileNotFoundError(f"bigwig file not found: {bigwig}")

    logger.info("extracting bigwig %s for %d copies", os.path.basename(bigwig), len(meta_id))
    bigwig_out = []
    with pyBigWig.open(bigwig) as bigwig_file:
        for i in range(len(chrom_list)):
            bigwig_out.append(
                extract_bigwig(
                    bigwig_file, chrom_list[i], start_list[i], end_list[i], strand_list[i]
                )
            )

    return save_or_return(bigwig_out, save_to_file, coordinate_table, "bigwig_out.lzma", "bigwig_out")
