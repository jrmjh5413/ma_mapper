"""Count reads, or smooth them into a profile, per base of each copy.

Two families of readout:

``read_*``
    the number of read 5' ends at each base, which is what a ChIP-exo or
    MNase footprint is read from;
``normal_*``
    the same ends smoothed with a Gaussian of half-width ``probe_length``,
    which is easier to look at but not a count of anything.

Each comes in ``forward``, ``reverse``, ``min``, ``max`` and (for reads)
``sum`` flavours, where forward and reverse are relative to the *copy's* strand,
not the genome's.
"""

import math
import os
import sys

import numpy as np
import pysam

from . import logger
from ._io_common import group_coordinates, read_coordinate_table, save_or_return

if sys.version_info >= (3, 8, 0):
    from typing import List, Literal
else:  # pragma: no cover - Python 3.7 and older
    from typing_extensions import List, Literal


def normal_array(width=1, sigma=1, odd=0):
    """An unnormalised Gaussian of half-width ``width``."""
    sigma2 = sigma**2

    def normal_func(x):
        return math.exp(-x * x / (2 * sigma2))

    # width is the half of the distribution
    values = list(map(normal_func, range(-width, width + odd)))
    return np.array(values)


_FORMAT = Literal[
    "read_max",
    "read_min",
    "read_forward",
    "read_reverse",
    "read_sum",
    "normal_max",
    "normal_min",
    "normal_forward",
    "normal_reverse",
]
BAM_FORMATS = (
    "read_max",
    "read_min",
    "read_forward",
    "read_reverse",
    "read_sum",
    "normal_max",
    "normal_min",
    "normal_forward",
    "normal_reverse",
)


def extract_bam(
    bam_file: str,
    chrom: str,
    start_list: List,
    end_list: List,
    strand: str,
    bam_format: _FORMAT,
    offset: int = 5,
    probe_length: int = 100,
    smoothing_length: int = 100,
    drop_dup=False,
    drop_unmap=False,
    drop_2=False,
    drop_multimap=False,
) -> np.ndarray:
    """Per-base read counts or smoothed profile for one copy, element orientation.

    Parameters
    ----------
    bam_file : str or pysam.AlignmentFile
        A path is opened and closed here; an open handle is left alone.
    chrom : str
        Chromosome name as spelled in the BAM.
    start_list, end_list : list of int
        Fragment coordinates, 0-based half-open (BED), ascending.
    strand : {'+', '-'}
        ``'-'`` reverses the output and swaps the forward and reverse
        profiles, so "forward" always means "along the element".
    bam_format : str
        One of ``BAM_FORMATS``. An unknown value raises ``ValueError``; before
        0.9.0 it raised ``UnboundLocalError`` from the next line.
    offset : int
        Bases to shift each read's 5' end by, 5 by default. **This is not
        zero**: the default silently moves every count 5 bp downstream, which
        suits a Tn5-style correction and does not suit anything else. Pass
        ``offset=0`` for raw 5' ends.
    probe_length : int
        Half-width of the Gaussian for the ``normal_*`` formats, and the
        padding applied around each fragment.
    smoothing_length : int
        Sigma of that Gaussian.
    drop_dup, drop_unmap, drop_2, drop_multimap : bool
        Skip reads flagged as duplicate, unmapped, secondary/supplementary, or
        carrying ``NH`` > 1. ``drop_multimap`` only works on aligners that
        write the ``NH`` tag (STAR, HISAT2 do; bwa and bowtie2 do not).

    Returns
    -------
    numpy.ndarray
        Length ``sum(end - start)``. ``read_*`` formats return uint32 counts,
        ``normal_*`` return float32.
    """
    if bam_format not in BAM_FORMATS:
        raise ValueError(f"bam_format must be one of {BAM_FORMATS}, got {bam_format!r}")

    close_after = False
    if isinstance(bam_file, str):
        bam_file = pysam.AlignmentFile(bam_file, "rb")
        close_after = True

    try:
        normal = normal_array(width=probe_length, sigma=smoothing_length, odd=1).astype(np.float32)

        result_array = []
        for i in range(len(start_list)):
            start = int(start_list[i])
            end = int(end_list[i])
            fragment_length = end - start
            padded = fragment_length + (2 * probe_length)
            # float32 and uint32, not float16 and uint8. A uint8 count wraps
            # silently to 0 at 256 reads, which a ChIP-exo summit reaches
            # routinely; float16 carries three decimal digits and saturates at
            # 65504, so a few hundred summed Gaussians lost precision.
            profile_normals_forward = np.zeros(padded, dtype=np.float32)
            profile_normals_reverse = np.zeros(padded, dtype=np.float32)
            profile_reads_forward = np.zeros(padded, dtype=np.uint32)
            profile_reads_reverse = np.zeros(padded, dtype=np.uint32)
            mapped_reads = bam_file.fetch(chrom, start, end)

            for mapped_read in mapped_reads:
                if drop_dup and mapped_read.is_duplicate:
                    continue
                if drop_unmap and mapped_read.is_unmapped:
                    continue
                if drop_2 and (mapped_read.is_secondary or mapped_read.is_supplementary):
                    continue
                if drop_multimap and mapped_read.has_tag("NH") and mapped_read.get_tag("NH") > 1:
                    continue
                if mapped_read.is_reverse is False:
                    read_position = mapped_read.reference_start - start + offset + probe_length
                else:
                    read_position = mapped_read.reference_end - start - offset - 1 + probe_length

                start_in_window = read_position - probe_length
                end_in_window = read_position + probe_length + 1

                if (start_in_window < 0) or (end_in_window > padded):
                    continue

                read_is_along_element = bool(mapped_read.is_reverse) == (strand == "-")
                if read_is_along_element:
                    np.add(
                        profile_normals_forward[start_in_window:end_in_window],
                        normal,
                        out=profile_normals_forward[start_in_window:end_in_window],
                    )
                    profile_reads_forward[read_position] += 1
                else:
                    np.add(
                        profile_normals_reverse[start_in_window:end_in_window],
                        normal,
                        out=profile_normals_reverse[start_in_window:end_in_window],
                    )
                    profile_reads_reverse[read_position] += 1

            if bam_format == "normal_min":
                output = np.minimum(profile_normals_forward, profile_normals_reverse)
            elif bam_format == "normal_max":
                output = np.maximum(profile_normals_forward, profile_normals_reverse)
            elif bam_format == "read_min":
                output = np.minimum(profile_reads_forward, profile_reads_reverse)
            elif bam_format == "read_max":
                output = np.maximum(profile_reads_forward, profile_reads_reverse)
            elif bam_format == "normal_forward":
                output = profile_normals_forward
            elif bam_format == "normal_reverse":
                output = profile_normals_reverse
            elif bam_format == "read_forward":
                output = profile_reads_forward
            elif bam_format == "read_reverse":
                output = profile_reads_reverse
            else:  # read_sum
                output = profile_reads_forward + profile_reads_reverse

            if strand == "-":
                output = np.flip(output)
            # The padding is symmetric, so the same slice applies before and
            # after the flip.
            output = output[probe_length : probe_length + fragment_length]
            result_array.append(output)

        if result_array:
            return np.concatenate(result_array)
        dtype = np.float32 if bam_format.startswith("normal") else np.uint32
        return np.zeros(int(np.sum(np.subtract(end_list, start_list))), dtype=dtype)
    finally:
        if close_after:
            bam_file.close()


def bam_io(
    coordinate_table,
    bam_file: str,
    bam_format: _FORMAT,
    save_to_file: bool = False,
    generate_new_id: bool = False,
    **kwargs,
) -> List:
    """Extract a BAM readout for every copy in a coordinate table.

    Parameters
    ----------
    coordinate_table : pandas.DataFrame or str
        BED-like table, one or more rows per copy.
    bam_file : str
        Path to an indexed BAM.
    bam_format : str
        One of ``BAM_FORMATS``.
    save_to_file : bool or str
        ``True`` writes ``bam_<format>.p`` beside the coordinate table; a
        string is used as the output path.
    generate_new_id : bool
        Treat every row as its own copy instead of grouping by column 3.
    **kwargs
        ``offset``, ``probe_length``, ``smoothing_length``, ``drop_dup``,
        ``drop_unmap``, ``drop_2`` and ``drop_multimap`` are passed through to
        :func:`extract_bam`. Note ``offset`` defaults to 5, not 0.

    Returns
    -------
    list of numpy.ndarray
        One array per copy, or the output path if the result was saved.
    """
    if bam_format not in BAM_FORMATS:
        raise ValueError(f"bam_format must be one of {BAM_FORMATS}, got {bam_format!r}")
    if not os.path.isfile(bam_file):
        raise FileNotFoundError(f"bam file not found: {bam_file}")

    coordinate_local = read_coordinate_table(coordinate_table)
    meta_id, chrom_list, start_list, end_list, strand_list = group_coordinates(
        coordinate_local, generate_new_id=generate_new_id
    )

    passthrough = (
        "offset",
        "probe_length",
        "smoothing_length",
        "drop_dup",
        "drop_unmap",
        "drop_2",
        "drop_multimap",
    )
    eb_kwargs = {key: value for key, value in kwargs.items() if key in passthrough}
    unknown = set(kwargs) - set(passthrough)
    if unknown:
        logger.warning("bam_io ignoring unknown arguments: %s", sorted(unknown))

    logger.info("extracting %s from %d copies", bam_format, len(meta_id))
    bam = []
    with pysam.AlignmentFile(bam_file, "rb") as handle:
        for i in range(len(chrom_list)):
            bam.append(
                extract_bam(
                    handle,
                    chrom_list[i],
                    start_list[i],
                    end_list[i],
                    strand_list[i],
                    bam_format,
                    **eb_kwargs,
                )
            )

    return save_or_return(
        bam, save_to_file, coordinate_table, f"bam_{bam_format}.lzma", f"bam {bam_format}"
    )
